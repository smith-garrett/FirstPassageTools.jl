# firstpassagedistributions.jl

using LinearAlgebra
using Distributions
import Distributions: pdf, logpdf, cdf, minimum, maximum, mean, var, @distr_support
using Roots  # needed for numerically finding quantiles
using ExponentialUtilities  # More numerically stable than built-in exp(::Matrix)

"""
    fpdistribution(T, A, p0)

Creates a new first-passage time distribution from the transient transition rate matrix T,
the absorbing transition rate matrix A, and the vector of initial probabilities of the
transient states p0.

Checks to make sure that the inputs meet the requirements for being a first passage time
distribution for a continuous-time, discrete-state Markov process.
"""
struct fpdistribution{T1, T2, T3} <: ContinuousUnivariateDistribution
    T::T1  # transient matrix
    A::T2  # absorbing matrix
    p0::T3  # initial condition

    # Internal constructor function
    fpdistribution(T::T1, A::T2, p0::T3) where {T1, T2, T3} = begin
        # At least one column of T sum to be less than zero
        @assert any(sum(T, dims=1) .< 0) "Transient T matrix incorrect:\n$T\n"
        # Checks to make sure that the total exit rate from T is equal to the total rate of
        # entry into A
        @assert isapprox(vec(-ones(1, size(T, 1)) * T), vec(A'*ones(size(A, 1)))) "Absorbing matrix A incorrect"
        # Size of p0 should correspond to the number of transient states
        @assert size(p0, 1) == size(T, 1) == size(A, 2) "Dimension mismatch with T, A, and/or p0"
        # p0 should be a probability distribution
        @assert isapprox(sum(p0), 1.0) "Initial condition p0 not a probability distribution"
        new{T1, T2, T3}(T, A, p0)
    end
end

# Helper functions
Distributions.@distr_support fpdistribution 0 +Inf
Distributions.minimum(d::fpdistribution) = 0.0
Distributions.maximum(d::fpdistribution) = +Inf
Distributions.mean(d::fpdistribution) = -sum(inv(d.T) * d.p0)
Distributions.var(d::fpdistribution) = 2*sum(d.T^(-2) * d.p0) - mean(d)^2

"""
    pdf(d::fpdistribution, t::Real)

Return the probability density function of d evaluated at the value t.
"""
Distributions.pdf(d::fpdistribution, t::T) where T<:Real = begin
    #insupport(d, t) ? sum(d.A * exp(t * d.T) * d.p0) : zero(t)
    insupport(d, t) ? sum(d.A * expv(t, d.T, d.p0)) : zero(t)
end

"""
    pdf(d::fpdistribution, t, dims)

Return the conditional probability density at `t` for the absorbing dimensions given in
`dims`.
"""
Distributions.pdf(d::fpdistribution, t::T, dims) where T<:Real = begin
    #insupport(d, t) ? getindex(d.A * exp(t * d.T) * d.p0 ./ splittingprobabilities(d), dims) : zero(t)
    insupport(d, t) ? getindex(d.A * expv(t, d.T, d.p0) ./ splittingprobabilities(d), dims) : zero(t)
end

"""
    logpdf(d::fpdistribution, t)

Returns the log probability density of `d` evaluated at `t`.
"""
Distributions.logpdf(d::fpdistribution, t::T) where T<:Real = begin
    log(pdf(d, t))
end


"""
    logpdf(d::fpdistribution, t::Real)

Returns the log probability density of `d` evaluated at `t` for the absorbing dimensions in
`dims`.
"""
Distributions.logpdf(d::fpdistribution, t::T, dims) where T<:Real = begin
    log.(pdf(d, t, dims))
end

"""
    cdf(d::fpdistribution, t::Real)

Returns the cumulative distribution function of `d` evaluated at `t`.
"""
Distributions.cdf(d::fpdistribution, t::T) where T<:Real = begin
    insupport(d, t) ? 1 - sum(expv(t, d.T, d.p0)) : zero(t)
end

Distributions.cdf(d::fpdistribution, t::T, dims) where T<:Real = begin
    error("cdf not implemented for multiple absorbing states.")
end

"""
    quantile(d::fpdistribution, p)

Return the ``p``-th quantile for the first-passage time distribution `d`.

The quantile function for these distributions has no closed-form solution, so this method
uses the `find_zero` function from the `Roots` package to find the quantile numerically.
"""
Distributions.quantile(d::fpdistribution, p) = begin
    if p <= zero(p)
        return zero(p)
    elseif p >= one(p)
        return Inf
    else
        # Use the mean of the distribution as the initial guess
        find_zero((x -> cdf(d, x) - p), mean(d))
    end
end

Distributions.quantile(d::fpdistribution, p, dims) = begin
    error("Quantile function not implemented for multiple absorbing states.")
end

Distributions.rand(d::fpdistribution, rng::AbstractVector{<:Real}) = begin
    quantile(d, rng)
end

"""
    splittingprobabilities(T, A, p0)

Compute the splitting probabilities, i.e., the probabilities of being absorbed into each of
the absorbing states.
"""
function splittingprobabilities(d::fpdistribution)
    return -d.A * inv(d.T) * d.p0
end

"""
    quasistationary(d::fpdistribution)

Return the quasi-stationary distribution of the distribution `d`, i.e., the probability
distribution over transient states conditional on not yet being absorbed.
"""
function quasistationary(d::fpdistribution)
    vals, vecs = eigen(d.T)
    if length(findall(vals .== maximum(vals))) != 1
        @warn "The largest eigenvalue was not unique. Quasi-stationary distribution might not make sense."
    end
    idx = argmax(vals)
    vec = vecs[:,idx]
    return vec / sum(vec)
end

