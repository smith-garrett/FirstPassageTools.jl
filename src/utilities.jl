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
using CSV
using DataFrames

# New idea for how to set up a system: .csv file w/ colums: cond,from,to(,rate)
# The rate column could be left off for now, since all we want to do is estimate τ. That
# means that the rate can just be set to #(states per condition) 

"""
    setup(path::String)

Setup a first-passage time problem by specifying a path to a .csv file containing the
transition.

The .csv file should have the column names "condition", "from", and "to". 
"""
function setup(path::String; return_names=false)
    df = CSV.read(path, DataFrame, stringtype=String)

    nconditions = length(unique(df.condition))
    statespercondition = Dict()
    for cond = unique(df.condition)
        statespercondition[cond] = nrow(df[isequal.(df.condition, cond), :])
    end

    statenames = []
    ctr = 1
    nametoidx = Dict()
    rows = []
    columns = []
    for row = 1:nrow(df)
        if df[row, :from] ∉ statenames
            push!(statenames, df[row, :from])
            nametoidx[df[row, :from]] = ctr
            ctr += 1
        end

        if df[row, :to] ∉ statenames
            push!(statenames, df[row, :to])
            nametoidx[df[row, :to]] = ctr
            ctr += 1
        end
        push!(rows, [df[row, :condition], nametoidx[df[row, :to]]])
        push!(columns, [df[row, :condition], nametoidx[df[row, :from]]])
    end

    nstates = length(statenames)
    W = zeros(nstates, nstates)
    for col in columns
        for row in rows
            W[row[2], col[2]] = statespercondition[row[1]]
        end
    end
    W = setdiagonal!(W)

    @assert 0 in W[diagind(W)] "No absorbing states provided. Check transitions."

    if return_names
        return W, statenames
    else
        return W
    end
end

"""
    setdiagonal!(W::Matrix)

Set the diagonal correctly for W-matrices (van Kampen, 2007, Ch. V.2)
"""
function setdiagonal!(W::Matrix)
    W[diagind(W)] .= 0.0
    W[diagind(W)] = -sum(W, dims=1)
    return W
end

"""
    rescale!(W::Matrix, τ)

Rescale a matrix a W-matrix (van Kampen, 2007, Ch. V.2) by a constant τ.
"""
function rescale!(W::Matrix{Float64}, τ)
    W *= τ
    return setdiagonal!(W)
end
