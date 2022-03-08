module SumSpace

using LinearAlgebra

export solvesvd, collocation_points, riemann, evaluate, expansion,
            supporter_functions, interpolate_supporter_functions

include("frame.jl")
include("cft.jl")


end # module
