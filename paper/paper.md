---
title: "FirstPassageTools.jl: A Julia package for first-passage time problems for
continous-time, discrete-state Markov processes."
tags:
  - Julia
  - Stochastic processes
  - Markov processes
  - First-passage time problems
authors:
  - name: Garrett Smith
    orcid: 0000-0001-7238-3942
    affiliation: University of Potsdam, Potsdam, Germany
date: 28 April 2022
bibliography: paper.bib

---

# Summary

In stochastic systems, one is often interested in how long it takes a system to reach a
particular state, the *first-passage time* from an initial transient state to an absorbing
state. Once the system reaches an absorbing state, it cannot transition out of it again.
There are examples of such processes in many scientific fields, for example, radioactive
decay [@vankampen2007stochastic], biological processes like the passage of molecules through
a cell membrane [@iyer-biswas2016first], or the transition between products and reactants in
chemical reactions [@swinburne2020rare]. First-passage times have also found use in
cognitive science: [@diederich2003simple] use discretized continuous-time, continuous-space
random walks to model perceptual decision making, and [@smith2021software] propose a
framework for modeling word-by-word reading times using continuous-time, discrete-state
Markov processes. `FirstPassageTools.jl` is a Julia package for sampling from and fitting
such first-passage time distributions for continuous-time, discrete-state Markov processes. 

In general, the dynamics of continuous-time, discrete-state Markov processes are governed by
a transition rate matrix that specifies the rate per unit time that the process can jump
from one state to another. `FirstPassageTools.jl` focuses on cases where the structure of
the transition rate matrices (number of states, allowed transitions, etc.) is fixed in
advance. Transition rate matrices can be read from a .csv file via the `setup()` function,
or they can be constructed by hand. This makes it possible to fit the transition rates of a
first-passage time distribution from a hypothesized process to empirical data, and compare
the fit to the data between alternative generating processes. `FirstPassageTools.jl` builds
on the well-established `Distributions.jl` and `Turing.jl` packages, making Julia tools for
data analysis and visualization accessible for a new type of probability distribution.

# Statement of need

To the best of my knowledge, there are currently no other packages that provide an interface
to first-passage time statistics for general continuous-time, discrete-state, absorbing
Markov processes. There are a number of other software packages that focus on specific
subtypes of this process, however, they do not have the same generality or focus on
parameter estimation as `FirstPassageTools.jl`. They are:

- [`MarkovChains.jl`](https://github.com/mfornino/MarkovChains.jl) and
  [`DiscreteMarkovChains.jl`](https://github.com/Maelstrom6/DiscreteMarkovChains.jl):
  Simulating state-to-state transitions of Markov chains
- [`EMpht.jl`](https://github.com/Pat-Laub/EMpht.jl): Fitting phase-type models, a specific
  type of continuous-time, discrete state Markov process, using the expectation maximization
  algorithm
- [`DiffModels.jl`](https://github.com/DrugowitschLab/DiffModels.jl): Simulating
  one-dimensiona, continuous-time, continuous-state drift-diffusion models with two
  absorbing boundaries
- [`PhaseTypeR`](https://rivasiker.github.io/PhaseTypeR/index.html): Simulate and fit
  phase-type distributions in R
- `DifferentialEquations.jl`: General simulation framework that includes continuous-time
  jump processes [@rackauckas2017differentialequations].

`FirstPassageTools.jl` aims to fill a gap in the existing software landscape for researchers
seeking to implement, fit, and compare first-passage time models with very few assumptions
about the allowed transitions between discrete transient states.

Currently, the focus is on models with a single absorbing state; however, basic functions
are implemented for systems with more than one absorbing state, e.g., a function for
calculating the splitting probabilities (the probability of being absorbed in a particular
absorbing state given an initial state) and conditional first-passsage time densities (the
probability density for how long it takes to be absorbed in a particular absorbing state
conditional on not being absorbed in another state).

Installation instructions and a basic example are provided in the repository's [`README`
file](https://github.com/garrett-m-smith/FirstPassageTools.jl). Further examples of how to
set up, sample, and fit first-passage time distributions are provided in the `notebooks`
directory of the repository.

# Acknowledgements

I would like to thank Prof. Dr. Shravan Vasishth, members of the Vasishth lab at the
University of Potsdam, and the audience and reviewers of 54th Annual Meeting of the Society
for Mathematical Psychology for helpful feedback. This work was partially funded by the
Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – Project ID 317633480
– SFB 1287

# References
