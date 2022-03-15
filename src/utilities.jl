# utilities.jl
# Tested with Julia 1.6 on macOS 11.6.4

# Notes:
# 1. make sure that the p0 stored in the structure has the same number of dimensions as the
#    transient matrix

# To do: fix rescale!() fn. to work for both square (W and T) and non-square (A) matrices

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

TO DO: return state names, separated by transient/absorbing.
"""
function setup(path::String)#; return_names=false)
    df = CSV.read(path, DataFrame, stringtype=String)
    @assert names(df) == ["condition", "from", "to"] "Incorrect column names in .csv file."

    matrixlist = []
    allnames = []
    for grp in groupby(df, :condition)
        curr_cond = grp[1, :condition]
	    statenames = []
	    ctr = 1
	    nametoidx = Dict()
	    rows = []
	    columns = []
        # Getting a mapping from state name to its index in the W matrix
	    for row = 1:nrow(grp)
	        if grp[row, :from] ∉ statenames
	            push!(statenames, grp[row, :from])
	            nametoidx[grp[row, :from]] = ctr
	            ctr += 1
	        end
	
	        if grp[row, :to] ∉ statenames
	            push!(statenames, grp[row, :to])
	            nametoidx[grp[row, :to]] = ctr
	            ctr += 1
	        end
            push!(rows, nametoidx[grp[row, :to]])
            push!(columns, nametoidx[grp[row, :from]])
	    end
	
        # Setting transition rates
	    nstates = length(statenames)
	    W = zeros(nstates, nstates)
        for (col, row) in zip(columns, rows)
            W[row, col] = nstates
        end
	    W = setdiagonal!(W)

        @assert 0 in W[diagind(W)] "No absorbing states provided for condition $(curr_cond)"
        push!(matrixlist, W)
        push!(allnames, statenames...)
    end

    # Separating transient and absorbing submatrices
    full = cat(matrixlist..., dims=(1, 2))
    (transient, absorbing) = getdims(full)
    T = full[transient, transient]
    A = full[absorbing, transient]

#    if return_names
#        return T, A, allnames
#    else
    return T, A
#    end
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
