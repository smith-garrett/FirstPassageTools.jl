#'---
#'author: Garrett Smith
#'title: Parameter recovery with a hierarchical model
#'---
#'
#' # Parameter recovery test using a hierarchical model
#'
#' The goal of this script is to test the recovery of parameters from a hierarchical model
#with a `fpdistribution` likelihood.

using Distributions
using Turing
using Plots, StatsPlots
using Pkg
Pkg.activate("../../FirstPassageTools.jl/")
using FirstPassageTools

#' ## Generating the true data
#'
#' First, we set up the transition rate matrices for the first-passage time distribution we
#' want to fit.

T = 4*[-1 0 0; 1 -1 1; 0 1 -2]
A = 4*[0 0 1]
p0 = [1.0, 0, 0]

#' Scaling the transition rate matrices by `τ = 2.5` should give mean first-passage times of
#' around 400ms. Generating and fitting the paramters (τ and the separate τᵢ) will be done
#' on the log scale and then exponentiated in order to keep the transition rates positive.

nparticipants = 25
true_tau = log(2.5)
true_sd = 0.25
true_tau_i = rand(Normal(0, true_sd), nparticipants)

#' The data will be saved in wide format: Each participant's data corresponds to a row, and
#' each column is a data point.

ndata = 40
data = zeros(nparticipants, ndata)
param = exp.(true_tau .+ true_tau_i)
for i = 1:nparticipants
    data[i,:] = rand(fpdistribution(param[i]*T, param[i]*A, p0), ndata)
end

#' ## Specifying the model
#'
#' First, we specify the priors on the log scale:

pr_tau = Normal(1, 0.25)
pr_sd = Exponential(0.5)  # Prior on the SD of the τᵢ

#' Next, we can write the full model including the likelihood

@model function mod(y)
    np = size(y, 1)
    nd = size(y, 2)
    # Priors
    τ ~ pr_tau
    # Need to initialize as a TArray b/c PG sampler
    τᵢ = tzeros(Float64, np)
    sd ~ pr_sd
    τᵢ ~ filldist(Normal(0, sd), np)

    # Likelihood
    mult = exp.(τ .+ τᵢ)
    #y ~ filldist(arraydist([fpdistribution(mult[p]*T, mult[p]*A, p0) for p in 1:np]), nd)
    for d = 1:nd
        for p = 1:np
            y[p, d] ~ fpdistribution(mult[p]*T, mult[p]*A, p0)
        end
    end
end

#' ## Sampling
#'
#' Here, we'll use the particle Gibbs sampler with 50 particles to sample from the
#' posterior. We'll use four chains of 1000 samples each. Make sure to execute this script
#' with `julia -t 4 HierarchicalParameterRecovery.jl`.

#posterior = sample(mod(data), PG(50), MCMCThreads(), 500, 4)
posterior = sample(mod(data), SMC(200), MCMCThreads(), 500, 4)

#' ## Evaluating parameter recovery
#' 
#' First, we summarize the chains:

print(describe(posterior))

#' And plot them:

histogram(posterior[:τ][:], xlabel="τ")
vline!([true_tau], label="True value")
savefig("tau_posterior.pdf")

histogram(posterior[:sd][:], xlabel="Std. deviation")
vline!([true_sd], label="True value")
savefig("sd_posterior.pdf")

plot_vec = []
for p = 1:nparticipants
    curr = posterior.name_map.parameters[p]
    plt = histogram(posterior[curr][:])
    plt = vline!(plt, [true_tau_i[p]])
    push!(plot_vec, plt)
end
plot(plot_vec..., legend=false, link=:all)
savefig("tau_i_posterior.pdf")

#' If the posterior contains the true values of the parameters, we can say the parameters
#' were recovered successfully.
#'
#' Let's also look at the Gelman-Rubin statistic for the chains:

print(gelmandiag(posterior))

