# utilities.jl

using LinearAlgebra
using CSV
using DataFrames

"""
    setup(path::String)

Setup a first-passage time problem by specifying a path to a .csv file containing the
transition.

The .csv file should have the column names "condition", "from", and "to". Currently, the
transition rates are set to be equal to the number of states in each condition. This keeps
the processing times comparable between conditions, and when fitting the scaling parameter
τ, that τ represents the transition rate per jump.
"""
function setup(path::String)
    df = CSV.read(path, DataFrame, stringtype=String)
    @assert names(df) == ["from", "to"] "Incorrect column names in .csv file."

    statenames = []
    ctr = 1
    nametoidx = Dict()
    rows = []
    columns = []
    # Getting a mapping from state name to its index in the W matrix
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
        push!(rows, nametoidx[df[row, :to]])
        push!(columns, nametoidx[df[row, :from]])
    end

    # Setting transition rates
    nstates = length(statenames)
    W = zeros(nstates, nstates)
    for (col, row) in zip(columns, rows)
        W[row, col] = nstates
    end
    W = setdiagonal!(W)

    @assert 0 in W[diagind(W)] "No absorbing states provided."

    # Separating transient and absorbing submatrices
    (transient, absorbing) = getdims(W)
    T = W[transient, transient]
    A = W[absorbing, transient]

    return T, A
end

"""
    setdiagonal!(W::Matrix)

Set the diagonal correctly for W-matrices (van Kampen, 2007, Ch. V.2).
"""
function setdiagonal!(W::Matrix)
    W[diagind(W)] .= 0.0
    W[diagind(W)] = -sum(W, dims=1)
    return W
end

"""
    rescale!(W::Matrix, τ)

Rescale a matrix a matrix by a constant τ.
"""
function rescale!(W::Matrix{Float64}, τ)
    W *= τ
end

"""
    getdims(W::Matrix)

Takes a W-matrix and returns the indices of the transient and absorbing states.
"""
function getdims(W::Matrix)
    transient = []
    absorbing =[]
    for i = 1:size(W, 1)
        if W[i, i] == 0
            push!(absorbing, i)
        else
            push!(transient, i)
        end
    end
    return transient, absorbing
end
