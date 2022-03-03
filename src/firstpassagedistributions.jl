# firstpassagedistributions.jl
# Tested with Julia 1.6 on macOS 11.6.4

using LinearAlgebra
using Distributions
using Roots  # needed for numerically finding quantiles
using ExponentialAction  # needed for differentiable expv()

"""
    splittingprobabilities(T, A, p0)

Compute the splitting probabilities, i.e., the probabilities of being absorbed into each of
the absorbing states.

Currently not very useful, because only exit times (unconditional first-passage tmes) are
implemented
"""
function splittingprobabilities(T::Matrix{Real}, A::Matrix{Real}, p0::Vector{Real})
    return -A * inv(T) * p0
end

# Creating new methods for getting a first passage time distributions
struct fpdistribution{T1::Matrix{Real}, T2::Matrix{Real}, T3::Vector{Real}} <: ContinuousUnivariateDistribution
    T::T1  # transient matrix
    A::T2  # absorbing matrix
    p0::T3  # initial condition
end

# Helper functions
Distributions.@distr_support fpdistribution 0 +Inf
Distributions.minimum(d::fpdistribution) = 0.0
Distributions.maximum(d::fpdistribution) = +Inf
Distributions.mean(d::fpdistribution) = -sum(inv(d.T) * d.p0)
Distributions.var(d::fpdistribution) = 2*sum(d.T^(-2) * d.p0) - mean(d)^2

# Probability density and cumulative density functions
Distributions.pdf(d::fpdistribution, t::Float64) = begin
    t >= 0 ? d.A * expv(t, d.T, d.p0) : zero(t)
end

Distributions.logpdf(d::fpdistribution, t::Float64) = begin
	log(pdf(d, t))
end

Distributions.cdf(d::fpdistribution, t::Float64) = begin
    t >= 0 ? 1 - sum(expv(t, d.T, d.p0)) : zero(t)
end

"""
    quantile(d::fpdistribution, p)

Return the `$p$`-th quantile for the first-passage time distribution `d`.

The quantile function for these distributions has no closed-form solution, so this method
uses the `find_zero` function from the `Roots` package to find the quantile numerically.
"""
Distributions.quantile(d::fpdistribution, p) = begin
    find_zero((x -> cdf(d, x) - p), 0.1)
end

Distributions.rand(d::fpdistribution, rng::AbstractVector{<:Real}) = begin
    quantile(d, rng)
end