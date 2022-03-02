# utilities.jl
# Tested with Julia 1.6 on macOS 11.6.4
# To do: define struct for storing info about first passage problem, functions for setting
# up transition matrices, getting dimension names
# Notes:
# 1. make sure that the p0 stored in the structure has the same number of dimensions as the
#    transient matrix
# 2. If implementing all-conditions-at-once, block-style matrix combination of mulitple
# conditions, will need to make sure to rescale the block for each condition by the number
# of states invovled in that condition.

using LinearAlgebra

"""
    rescale!(W::Matrix, τ)

Rescale a matrix a W-matrix (van Kampen, 2007, Ch. V.2) by a constant τ.
"""
function rescale!(W::Matrix{Float64}, τ)
    W[diagind(W)] .= 0.0
    W *= τ
    W[diagind(W)] = -sum(W, dims=1)
    return W
end

# New idea for how to set up a system: .csv file w/ colums: cond,from,to(,rate)
# The rate column could be left off for now, since all we want to do is estimate τ. That
# means that the rate can just be set to #(states per condition) 