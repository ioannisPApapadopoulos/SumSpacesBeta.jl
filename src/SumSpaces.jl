module SumSpaces

using LinearAlgebra

include("frame.jl")
include("cft.jl")
include("assembly.jl")

export solvesvd, collocation_points, riemann, evaluate, expansion, framematrix, dualframematrix,
            supporter_functions, interpolate_supporter_functions, columns_supporter_functions,
            idmap_sum2dual, idmap_append2sum, idmap_append2dual, hilbertmap, diffmap, fractionalhelmholtzmap,
            sum_space, appended_sum_space, dual_sum_space


function sum_space(eT, ewU, u, xx, N)
    yy = (
        eT[xx,1:N+2] * u[1:N+2]
        + ewU[xx,1:N+1] * u[N+3:2*N+3]
    )
    return yy
end

function appended_sum_space(eT, ewU, yU_2, yU_1, ywT0, ywT1, u, xx, N)
    yy = (
        eT[xx,1:N+2] * u[1:N+2]
        + ewU[xx,1:N+1] * u[N+3:2*N+3]
        + yU_2(xx) .* u[2*N+4]
        + yU_1(xx) .* u[2*N+5]
        + ywT0(xx) .* u[2*N+6]
        + ywT1(xx) .* u[2*N+7]
    )
    return yy
end

function dual_sum_space(eU, ewT, u, xx, N)
    yy = (
        eU[xx,1:N+4] * u[1:N+4]
        + ewT[xx,1:N+3] * u[N+5:2*N+7]
    )
    return yy
end

end # module
