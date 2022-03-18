module SumSpaces

using LinearAlgebra

include("frame.jl")
include("cft.jl")
include("assembly.jl")

export solvesvd, collocation_points, riemann, evaluate, expansion_sum_space, framematrix, dualframematrix,          
            supporter_functions, interpolate_supporter_functions, columns_supporter_functions, inverse_fourier_transform, fractional_heat_fourier_solve,
            idmap_sum2dual, idmap_append2sum, idmap_append2dual, hilbertmap, diffmap, fractionalhelmholtzmap,
            affinetransform, sum_space, appended_sum_space, dual_sum_space



# Affine transform to scale and shift polys. 
function affinetransform(a,b,x)
    y = 2/(b-a) .* (x.-(a+b)/2)
end

function sum_space(eT, ewU, u, xx, N; a=[-1.,1.])
    elements = length(a)-1
    yy = (
        eT[affinetransform(a[1],a[2],xx),1:N+2] * u[1:N+2]
        + ewU[affinetransform(a[1],a[2],xx),1:N+1] * u[N+3:2*N+3]
    )
    for e in 2:elements
        yy += (
            eT[affinetransform(a[e],a[e+1],xx),2:N+2] * u[2*(e-1)*N+2*e:(2*e-1)*N+2*e]
            + ewU[affinetransform(a[e],a[e+1],xx),1:N+1] * u[(2*e-1)*N+2*e+1:2*e*N+2*e+1]
        )
    end
    return yy
end

function appended_sum_space(eT, ewU, yU_2, yU_1, ywT0, ywT1, u, xx, N; a=[-1.,1.])
    yy = (
        eT[affinetransform(a[1],a[2],xx),1:N+2] * u[1:N+2]
        + ewU[affinetransform(a[1],a[2],xx),1:N+1] * u[N+3:2*N+3]
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
