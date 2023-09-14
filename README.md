# FirstPassageTools

[![Build Status](https://github.com/smith-garrett/FirstPassageTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/garrett-m-smith/FirstPassageTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/smith-garrett/FirstPassageTools.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/garrett-m-smith/FirstPassageTools.jl)
[![DOI](https://zenodo.org/badge/465749757.svg)](https://zenodo.org/badge/latestdoi/465749757)

`FirstPassageTools` is a Julia package for setting up first-passage time distributions for
continuous-time, discrete-state Markov processes. The first-passage time distributions can
then be fit to empirical first-passage time data using [Turing.jl](https://turing.ml).

# Installation

FirstPassageTools.jl can be installed from Juila's central package repository:

```julia
julia> ]
pkg> add FirstPassageTools
```

FirstPassageTools.jl can also be installed directly from Github in the Julia REPL:

```julia
julia> ]
pkg> add https://github.com/smith-garrett/FirstPassageTools.jl
```

# Usage

To set up a first-passage time distribution, one needs to provide two transition rate
matrices and a vector with the initial probability distribution over the transient states
(`p0` in the example below). The first matrix (`T` below) determines the transition rates
between transient states of the system. The second determines the transition rates from the
transient states to one or more absorbing states (`A` below). The rows of the transient
matrix T should sum to the values given in the absorbing matrix. For both matrices, the
$i,j$-th entry should provide the transition rate from state $j$ to state $i$.

Currently, the most complete functionality is available for the first-passage time to *any*
absorbing state by time t. When a system has more than one absorbing state, one might be
interested in the conditional first-passage time to reach a particular absorbing state
before all others. Some methods are implemented for this, but not everything.

For a single transient state and a single absorbing state, the first-passage time
distribution is equivalent to the exponential distribution:

```julia
julia> using FirstPassageTools
julia> T = [-1.0;;]  # transition rate matrices need to be 2-dimensional
julia> A = [1.0;;]
julia> p0 = [1.0]
julia> fp = fpdistribution(T, A, p0)
```

Available methods for first-passage time distributions include `mean()`, `var()`, `rand()`,
`pdf()`, `logpdf()`, `cdf()`, and `quantile()`. From here, the transition rates of the
distribution can be fit to data. See the notebooks directory for additional tutorials and
parameter recovery exercises.

Because the first-passage time distributions here are sub-types of the continuous univariate
distribution from[`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl),
plotting functions from [`StatsPlots.jl`](https://github.com/JuliaPlots/StatsPlots.jl)
work out of the box. For example, to plot the probability density function and the
cumulative distribution function from the above sample, you can run:

```julia
julia> using Plots, StatsPlots
julia> plot(fp, label="PDF")
julia> plot!(fp, func=cdf, label="CDF")
```

Additional examples are provided in the notebooks directory. Verification of correct
sampling, at least for certain statistical models, is provided in the
`SimulationBasedCalibration.jl` script in the `notebooks` directory.


