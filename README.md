# FirstPassageTools

[![Build Status](https://github.com/garrett-m-smith/FirstPassageTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/garrett-m-smith/FirstPassageTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/garrett-m-smith/FirstPassageTools.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/garrett-m-smith/FirstPassageTools.jl)

`FirstPassageTools` is a Julia package for setting up first-passage time distributions for
continuous-time, discrete-state Markov processes. The first-passage time distributions can
then be fit to empirical first-passage time data using [Turing.jl](https://turing.ml).

# Installation

Currently, FirstPassageTools.jl can be installed directly from Github:

```julia
] add https://github.com/garrett-m-smith/FirstPassageTools.jl
```

Soon, it will be available via Julia's central package installation registry.


# Usage

To set up a first-passage time distribution, one needs to provide two transition rate
matrices and a vector with the initial probability distribution over the transient states.
The first matrix determines the transition rates between transient states of the system. The
second determines the transition rates from the transient states to one or more absorbing
states. The rows of the transient matrix T should sum to the values given in the absorbing
matrix.

(Currently, only the exit-time distribution is implemented, the probability of
reaching any absorbing state by time t.)

For a single transient state and a single absorbing state, the first-passage time
distribution is equivalent to the exponential distribution:

```julia
> using FirstPassageTools
> T = [-1;;]  # transition rate matrices need to be 2-dimensional
> A = [1;;]
> p0 = [1]
> fp = fpdistribution(T, A, p0)
```

From here, the distribution can be fit to data. Tutorials are coming soon!

Available methods for first-passage time distributions include `mean()`, `var()`, `rand()`,
`pdf()`, `logpdf()`, `cdf()`, and `quantile()`.
