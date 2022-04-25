#'---
#'author: Garrett Smith
#'title: Parameter recovery with a hierarchical model
#'---
#'
#' # Parameter recovery test using a hierarchical model
#'
#' The goal of this script is to test the recovery of parameters from a hierarchical model
#' with a `fpdistribution` likelihood.

#+ echo=false
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

#T = 4*[-1.0 0 0; 1 -1 1; 0 1 -2]
T = 4*[-1 1.0 0; 1 -2 1; 0 1 -2]
A = 4*[0 0 1.0]
p0 = [1.0, 0, 0]

#' Scaling the transition rate matrices by `τ = 2.5` should give mean first-passage times of
#' around 400ms. Generating and fitting the paramters (τ and the separate τᵢ) will be done
#' on the log scale and then exponentiated in order to keep the transition rates positive.

nparticipants = 50
true_tau = 2.5
true_sd = 0.2
true_tau_i = exp.(rand(Normal(0, true_sd), nparticipants));

#' The data will be saved in wide format: Each participant's data corresponds to a row, and
#' each column is a data point.

#+ echo=false
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

# Switch to param = exp(tau) + exp(tau_i). This will prevent really small params b/c multiplication. 
@model function mod(y)
    np, nd = size(y)
    # Priors
    # Using the non-centered parameterization for τ
    τ ~ Normal()
    τ̂ = 1 + 0.1*τ  # Corresponds to Normal(1, 0.1)
    sd ~ Exponential(0.25)
    τᵢ ~ filldist(Normal(), np)
    τ̂ᵢ = sd .* τᵢ  # Corresponds to MvNormal(0, sd)
    # Likelihood
    mult = exp.(τ̂ .+ τ̂ᵢ)
    y ~ filldist(arraydist([fpdistribution(mult[p]*T, mult[p]*A, p0) for p in 1:np]), nd)
    return τ̂, τ̂ᵢ, mult
end

#+ echo=false
#@model function mod_centered(y)
#    np, nd = size(y)
#    # Priors
#    # Using the centered parameterization for τ
#    τ ~ Normal(1, 0.1)
#    sd ~ Exponential(0.25)
#    τᵢ ~ filldist(Normal(0, sd), np)
#    # Likelihood
#    mult = exp.(τ .+ τᵢ)
#    y ~ filldist(arraydist([fpdistribution(mult[p]*T, mult[p]*A, p0) for p in 1:np]), nd)
#end
#
#@model function mod_noncentered_tau(y)
#    np, nd = size(y)
#    # Priors
#    # Using the centered parameterization for τ
#    τ ~ Normal()
#    τ̂ = 1.0 + 0.1*τ
#    sd ~ Exponential(0.25)
#    τᵢ ~ filldist(Normal(0, sd), np)
#    # Likelihood
#    mult = exp.(τ̂ .+ τᵢ)
#    y ~ filldist(arraydist([fpdistribution(mult[p]*T, mult[p]*A, p0) for p in 1:np]), nd)
#end
#
#@model function mod_noncentered_tau_i(y)
#    np, nd = size(y)
#    # Priors
#    # Using the centered parameterization for τ
#    τ ~ Normal(1, 0.1)
#    sd ~ Exponential(0.25)
#    τᵢ ~ filldist(Normal(), np)
#    τ̂ᵢ = sd .* τᵢ
#    # Likelihood
#    mult = exp.(τ .+ τ̂ᵢ)
#    y ~ filldist(arraydist([fpdistribution(mult[p]*T, mult[p]*A, p0) for p in 1:np]), nd)
#end
#' ## Sampling
#'
#' Here, we'll use the NUTS sampler with a burnin of 100 samples and an acceptance rate of 
#' 0.65 posterior. We'll use four chains of 1000 samples each. Make sure to execute this 
#' script with `julia -t 4 HierarchicalParameterRecovery.jl`.

#posterior = sample(mod(data), NUTS(250, 0.65), MCMCThreads(), 1000, 4)
posterior = sample(mod(data), NUTS(250, 0.7), 1000);
#posterior_centered = sample(mod_centered(data), NUTS(100, 0.65), 1000);
#posterior_noncentered_tau = sample(mod_noncentered_tau(data), NUTS(100, 0.65), 1000);
#posterior_noncentered_tau_i = sample(mod_noncentered_tau_i(data), NUTS(100, 0.65), 1000);

#' ## Evaluating parameter recovery
#' 
#' First, we summarize the chains:

posterior

#' Let's also look at the Gelman-Rubin statistic for the chains:
#gelmandiag(posterior)

# #' And the centered version:
#posterior_centered

# #' τ non-centered, τᵢ centered:
#posterior_noncentered_tau

# #' τ centered, τᵢ non-centered
#posterior_noncentered_tau_i

#' And plot histograms of the parameters on the millisecond scale:
#+ echo=false
gen = generated_quantities(mod(data), posterior);
gen = vcat(gen...);
gen_tau = [x[1] for x in gen];
gen_tau_i = transpose(hcat([x[2] for x in gen]...));
shrunken = transpose(hcat([x[3] for x in gen]...));

#+ echo=false; fig_cap="Posterior distribution of τ. Vertical line shows the true value"
histogram(exp.(gen_tau), xlabel="τ")
vline!([true_tau])

#+ echo=false; fig_cap="Posterior distribution of the SD. Vertical lines show the true value."
histogram(posterior[:sd][:], xlabel="Std. deviation of the τᵢ")
vline!([true_sd], label="True value")

#+ echo=false; fig_cap="Posterior distributions of the τᵢ. Vertical lines show the true value."
plot_vec = [];
for p = 1:nparticipants
    plt = histogram(exp.(gen_tau_i[:,p]), xticks=false, yticks=false)
    plt = vline!(plt, [true_tau_i[p]])
    push!(plot_vec, plt)
end
plot(plot_vec..., legend=false)

#+ echo=false; fig_cap="Posterior distributions of τ*τᵢ. Vertical lines show the true value."
plot_vec = [];
for p = 1:nparticipants
    plt = histogram(shrunken[:,p], xticks=false, yticks=false)
    plt = vline!(plt, [param[p]])
    push!(plot_vec, plt)
end
plot(plot_vec..., legend=false)

#' If the posterior contains the true values of the parameters, we can say the parameters
#' were recovered successfully.
#'
