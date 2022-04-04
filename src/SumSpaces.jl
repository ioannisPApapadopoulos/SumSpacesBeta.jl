module SumSpaces

using SpecialFunctions, LinearAlgebra, BlockBandedMatrices, BlockArrays, 
    ClassicalOrthogonalPolynomials, DomainSets, StaticArrays, ContinuumArrays, QuasiArrays,
    FillArrays, ArrayLayouts, LazyBandedMatrices, LazyArrays, FFTW, Interpolations, InfiniteArrays


import QuasiArrays: DefaultQuasiArrayStyle, cardinality
import Base: in, axes, getindex, broadcasted, tail, +, -, *, /, \, convert, OneTo, show, summary, ==, oneto
import ContinuumArrays: Weight, grid, ℵ₁, ℵ₀, @simplify, Basis, basis, @simplify, Identity, AbstractAffineQuasiVector, AbstractQuasiArray, AbstractQuasiMatrix
import ClassicalOrthogonalPolynomials: checkpoints, ldiv, paddeddata, jacobimatrix, orthogonalityweight, sqrtx2, Hilbert, Derivative
import BlockArrays: block, blockindex, Block, _BlockedUnitRange, BlockSlice
import BlockBandedMatrices: BlockTridiagonal, AbstractBlockBandedMatrix, blockbandwidths, subblockbandwidths
import InfiniteArrays: OneToInf

include("frame.jl")
include("cft.jl")
include("assembly.jl")
include("extendedchebyshev.jl")
include("sumspace.jl")

export solvesvd, collocation_points, riemann, evaluate, expansion_sum_space, framematrix, dualframematrix, split_block_helmholtz_matrix, split_block_helmholtz_vector,      
            supporter_functions, interpolate_supporter_functions, columns_supporter_functions, inverse_fourier_transform, fractional_heat_fourier_solve,
            idmap_sum2dual, idmap_append2sum, idmap_append2dual, hilbertmap, diffmap, fractionalhelmholtzmap, fractionallaplacemap,
            affinetransform, sum_space, appended_sum_space, dual_sum_space, 
            dual_sum_space2, sqrt_laplace_wT0, sqrt_laplace_U_1,
            ExtendedChebyshev, ExtendedChebyshevT, ExtendedChebyshevU, extendedchebyshevt, ExtendedWeightedChebyshevT, ExtendedWeightedChebyshevU,
            SumSpace, SumSpaceP, SumSpaceD, 
            Block, Derivative



# Affine transform to scale and shift polys. 
function affinetransform(a,b,x)
    y = 2/(b-a) .* (x.-(a+b)/2)
end

function sum_space(eT, ewU, u, xx, N; a=[-1.,1.])
    el_no = length(a)-1
    
    # Constant term
    yy = eT[affinetransform(a[1],a[2],xx),1] * u[1]

    for e in 1:el_no
        yy += (
            eT[affinetransform(a[e],a[e+1],xx),2:N+2] * u[2*(e-1)*N+2*e:(2*e-1)*N+2*e]
            + ewU[affinetransform(a[e],a[e+1],xx),1:N+1] * u[(2*e-1)*N+2*e+1:2*e*N+2*e+1]
        )
    end
    return yy
end

function appended_sum_space(eT, ewU, yU0, yU_1, ywT0, ywT1, u, xx, N; a=[-1.,1.])
    el_no = length(a)-1
    yy = eT[affinetransform(a[1],a[2],xx),1] * u[1]

    #FIXME: Explain why order is changed
    for e in 1:el_no
        yy += (
            eT[affinetransform(a[e],a[e+1],xx),2:N+2] * u[2*(e-1)*N+2*e:(2*e-1)*N+2*e]
            + ewU[affinetransform(a[e],a[e+1],xx),1:N+1] * u[(2*e-1)*N+2*e+1:2*e*N+2*e+1]
            + yU_1[e](xx) .* u[end-4*(el_no-e+1)+1]
            + yU0[e](xx) .* u[end-4*(el_no-e+1)+2]
            + ywT0[e](xx) .* u[end-4*(el_no-e+1)+3]
            + ywT1[e](xx) .* u[end-4*(el_no-e+1)+4]

        )
    end
    return yy
end

function dual_sum_space(eU, ewT, u, xx, N)
    yy = (
        eU[xx,1:N+4] * u[1:N+4]
        + ewT[xx,1:N+3] * u[N+5:2*N+7]
    )
    return yy
end

function dual_sum_space2(eU, ewT, u, xx, N)
    yy = (
        eU[xx,2:N+3] * u[1:N+2]
        + ewT[xx,1:N+2] * u[N+3:end]
    )
    return yy
end

sqrt_laplace_wT0 = xx -> [abs(x) <= 1 ? 0.5*log(4)-Base.MathConstants.eulergamma : 0.5*log(4)-Base.MathConstants.eulergamma-asinh(sqrt(x^2-1)) for x in xx]
sqrt_laplace_U_1 = xx -> [abs(x) <= 1 ? -asin(x) : -sign(x)*pi/2 for x in xx]

end # module
