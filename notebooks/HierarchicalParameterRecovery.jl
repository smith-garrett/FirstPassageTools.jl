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
using Zygote
Turing.setadbackend(:zygote)
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
true_tau = 2.5
true_sd = 0.25
true_tau_i = exp.(rand(Normal(0, true_sd), nparticipants))

#' The data will be saved in wide format: Each participant's data corresponds to a row, and
#' each column is a data point.

ndata = 10
data = zeros(nparticipants, ndata)
param = true_tau .* true_tau_i
for i = 1:nparticipants
    data[i,:] = rand(fpdistribution(param[i]*T, param[i]*A, p0), ndata)
end

#' ## Specifying the model
#'
#' Now, we can write the full model including the likelihood. Note that we're using a
#' non-centered parameterization. This is because pilot simulations suggested that
#' sampling was biased in the centered parameterization.

@model function mod(y)
    np, nd = size(y)
    # Priors
    # Using the non-centered parameterization for τ
    #τ ~ Normal()
    #τ̂ = 1 + 0.1*τ  # Corresponds to Normal(1, 0.1)
    τ ~ Normal(1, 0.1)
    sd ~ Exponential(0.25)
    τᵢ ~ filldist(Normal(), np)
    τ̂ᵢ = sd.*τᵢ  # Corresponds to MvNormal(0, sd)
    # Likelihood
    #mult = exp.(τ̂ .+ τ̂ᵢ)
    mult = exp.(τ .* τ̂ᵢ)
    y ~ filldist(arraydist([fpdistribution(mult[p]*T, mult[p]*A, p0) for p in 1:np]), nd)
    return τ̂, τ̂ᵢ
end

#' ## Sampling
#'
#' Here, we'll use the NUTS sampler with a burnin of 100 samples and an acceptance rate of 
#' 0.65 posterior. We'll use four chains of 1000 samples each. Make sure to execute this 
#' script with `julia -t 4 HierarchicalParameterRecovery.jl`.

posterior = sample(mod(data), NUTS(100, 0.65; init_ϵ=0.1), MCMCThreads(), 500, 4)

#' ## Evaluating parameter recovery
#' 
#' First, we summarize the chains:

posterior

#' And plot histograms of the parameters on the millisecond scale:
#+ echo=false
gen = generated_quantities(mod(data), posterior)
gen = vcat(gen...)
gen_tau = [x[1] for x in gen]
gen_tau_i = transpose(hcat([x[2] for x in gen]...))

#+ fig_cap="Posterior distribution of τ. Vertical lines show the true value."
histogram(exp.(gen_tau), xlabel="τ")
vline!([true_tau])

#' fig_cap="Posterior distribution of the SD. Vertical lines show the true value."
histogram(posterior[:sd][:], xlabel="Std. deviation of the τᵢ")
vline!([true_sd], label="True value")

#' fig_cap="Posterior distributions of the τᵢ. Vertical lines show the true value."
plot_vec = []
for p = 1:nparticipants
    plt = histogram(exp.(gen_tau_i[:,p]), xticks=false, yticks=false)
    plt = vline!(plt, [true_tau_i[p]])
    push!(plot_vec, plt)
end
plot(plot_vec..., legend=false)

#' If the posterior contains the true values of the parameters, we can say the parameters
#' were recovered successfully.
#'
#' Let's also look at the Gelman-Rubin statistic for the chains:

gelmandiag(posterior)

