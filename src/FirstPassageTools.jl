module FirstPassageTools

include("utilities.jl")
export rescale!
export setup

include("firstpassagedistributions.jl")
export splittingprobabilities
export fpdistribution
export @distr_support
export minimum
export maximum
export mean
export var
export pdf
export logpdf
export cdf
export quantile
export rand
export quasistationary

end
