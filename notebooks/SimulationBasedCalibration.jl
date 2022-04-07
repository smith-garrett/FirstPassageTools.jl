#' ---
#' title: Simulation-based calibration of FirstPassageTools.jl
#' author: Garrett Smith
#' date: 07.04.2022
#' ---
#'
#' # Simulation-based Calibration 
#'
#' Simulation-based calibration [@talts2018validating] is used to ensure that a Bayesian
#' sampling algorithm is functioning properly. This script does a reasonably thorough job of
#' testing the NUTS sampler available in `Turing.jl` for a simple, one-parameter,
#' non-hierarchical first-passage time model.

using Distributions
using Turing
using Plots
using Pkg
Pkg.activate("../../FirstPassageTools.jl/")
using FirstPassageTools
using Zygote
Turing.setadbackend(:zygote)

#' # Defining things
#' 
#' Here, we define a simple model with two transient states and a single absorbing state.

T = [-3.0 0; 3.0 -3.0]
A = [0.0 3]
p0 = [1.0, 0]

#' Here, I define the prior on τ, the scaling parameter to be fit. If the sampling algorithm
#' works as it's supposed to, then the posterior should eventually approach the prior when
#' repeatedly sampling τs from the prior, simulating data, and then fitting the model to
#' that data.

prior = truncated(Normal(1.0, 0.5), lower=0.0)

#' Now, we define the model, a combination of the prior and a likelihood.

@model function mod(y)
    # Define prior
    τ ~ prior
    
    # Define likelihood
    y ~ filldist(fpdistribution(τ*T, τ*A, p0), length(y))
end

#' Now, we define the SBC function, which takes a single argument that is the number of
#' times to run the prior-generate-fit sequence, the number of data points to generate, and
#' the number of posterior samples to sample.

function sbc(n=1000)
    ranks = zeros(n)
    #ranks = Array{Int64, 1}
    Threads.@threads for i = 1:n
        # Sample a τ from the prior
        curr = rand(prior)
        # Generate data from the fpdistribution
        dat = rand(fpdistribution(curr*T, curr*A, p0), n)
        # Fit the fpdistribution to the simulated data
        posterior = sample(mod(dat), NUTS(0.65), n, progress=false)
        # Get the rank of curr in the posterior
        ranks[i] = count(x -> x < curr, posterior[:τ])
        #push!(ranks, count(x -> x < curr, posterior[:τ]))
    end
    return ranks
end

#' # Running the calibration
#'
#' Now we run the calibration:

n = 10
rks = sbc(n)
println(rks)

#' # Visualizing the results
#'
#' If the sampler is well-calibrated, the whole histogram should fall within the confidence
#' interval in the plot.

function plot_sbc(ranks, n)
    bn = Binomial(length(ranks), (n+1)^-1)
    fig = histogram(ranks, bins=0:n)
    hspan!(fig, quantile(bn, [0.005, 0.995]), color=:grey, alpha=0.25)
    savefig(fig, "SBCResults.pdf")
end

plot_sbc(rks, n)
