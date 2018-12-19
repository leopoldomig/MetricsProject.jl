module MetricsProject

using LinearAlgebra, Statistics, NLsolve

using Random, Distributions
using Optim
using Optim: converged, maximum, maximizer, minimizer, iterations

# Including the modules
include("Berry1994.jl")
include("sdc.jl")

# Exporting the functions of the modules
export Berry1994, sdc

end # module
