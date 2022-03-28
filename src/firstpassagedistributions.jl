# firstpassagedistributions.jl
# Tested with Julia 1.6 on macOS 11.6.4

using LinearAlgebra
using Distributions
using Roots  # needed for numerically finding quantiles

# Creating new methods for getting a first passage time distributions
struct fpdistribution <: ContinuousUnivariateDistribution
    T::Matrix  # transient matrix
    A::Matrix  # absorbing matrix
    p0::Vector  # initial condition

    # Internal constructor function
    fpdistribution(T, A, p0) = begin
        # At least one column of T sum to be less than zero
        @assert any(sum(T, dims=1) .< 0) "Transient T matrix incorrect:\n$T\n"
        # Definition of A: BUT ONLY FOR SYSTEMS WITH A SINGLE ABSORBING STATE!
        @assert isapprox(vec(-ones(1, size(T, 1)) * T), vec(A'*ones(size(A, 1)))) "Absorbing matrix A incorrect"
        # Size of p0 should correspond to the number of transient states
        @assert size(p0, 1) == size(T, 1) == size(A, 2) "Dimension mismatch with T, A, and/or p0"
        # p0 should be a probability distribution
        @assert isapprox(sum(p0), 1.0) "Initial condition p0 not a probability distribution"
        new(T, A, p0)
    end
end

# Helper functions
Distributions.@distr_support fpdistribution 0 +Inf
Distributions.minimum(d::fpdistribution) = 0.0
Distributions.maximum(d::fpdistribution) = +Inf
Distributions.mean(d::fpdistribution) = -sum(inv(d.T) * d.p0)
Distributions.var(d::fpdistribution) = 2*sum(d.T^(-2) * d.p0) - mean(d)^2

# Probability density and cumulative density functions
Distributions.pdf(d::fpdistribution, t::Real) = begin
    t >= 0 ? sum(d.A * exp(t * d.T) * d.p0) : zero(t)
end

Distributions.logpdf(d::fpdistribution, t::Real) = begin
	log(pdf(d, t))
end

Distributions.cdf(d::fpdistribution, t::Real) = begin
    t >= 0 ? 1 - sum(exp(t * d.T) * d.p0) : zero(t)
end

"""
    quantile(d::fpdistribution, p)

Return the ``p``-th quantile for the first-passage time distribution `d`.

The quantile function for these distributions has no closed-form solution, so this method
uses the `find_zero` function from the `Roots` package to find the quantile numerically.
"""
Distributions.quantile(d::fpdistribution, p) = begin
    if p <= 0
        return 0.0
    elseif p >= 1.0
        return Inf
    else
        find_zero((x -> cdf(d, x) - p), p)
    end
end

Distributions.rand(d::fpdistribution, rng::AbstractVector{<:Real}) = begin
    quantile(d, rng)
end

"""
    splittingprobabilities(T, A, p0)

Compute the splitting probabilities, i.e., the probabilities of being absorbed into each of
the absorbing states.

Currently not very useful, because only exit times (unconditional first-passage tmes) are
implemented
"""
function splittingprobabilities(d::fpdistribution)
    return -d.A * inv(d.T) * d.p0
end