module SumSpaces

using SpecialFunctions, LinearAlgebra, BlockBandedMatrices, BlockArrays, 
    ClassicalOrthogonalPolynomials, StaticArrays, ContinuumArrays, DomainSets,
    FillArrays, LazyBandedMatrices, LazyArrays, FFTW, Interpolations, InfiniteArrays,
    QuasiArrays

# import QuasiArrays: DefaultQuasiArrayStyle, cardinality
import Base: in, axes, getindex, broadcasted, tail, +, -, *, /, \, convert, OneTo, show, summary, ==, oneto
import ContinuumArrays: Weight, grid, ℵ₁, ℵ₀, @simplify, Basis, basis, @simplify, Identity, AbstractAffineQuasiVector, AbstractQuasiArray, AbstractQuasiMatrix
import ClassicalOrthogonalPolynomials: sqrtx2, Hilbert, Derivative #checkpoints, ldiv, paddeddata, jacobimatrix, orthogonalityweight, 
import BlockArrays: block, blockindex, Block, _BlockedUnitRange#, BlockSlice
import BlockBandedMatrices: _BandedBlockBandedMatrix #BlockTridiagonal, AbstractBlockBandedMatrix, blockbandwidths, subblockbandwidths
import InfiniteArrays: OneToInf

include("frame.jl")
include("cft.jl")
include("extendedchebyshev.jl")
include("sumspace.jl")
include("element-sumspace.jl")

export solvesvd, collocation_points, riemann, evaluate, framematrix,
            supporter_functions, interpolate_supporter_functions, coefficient_supporter_functions, inverse_fourier_transform, fractional_heat_fourier_solve, fft_supporter_functions,
            affinetransform, sqrt_laplace_wT0, sqrt_laplace_U_1,
            ExtendedChebyshev, ExtendedChebyshevT, ExtendedChebyshevU, extendedchebyshevt, ExtendedWeightedChebyshevT, ExtendedWeightedChebyshevU,
            SumSpace, SumSpaceP, SumSpaceD, AppendedSumSpace, ElementSumSpace, ElementAppendedSumSpace,
            Block, Derivative, Fill, BlockArray,
            coefficient_interlace, coefficient_stack, Id_Sp_Sd
            



# Affine transform to scale and shift polys. 
function affinetransform(a,b,x)
    y = 2/(b-a) .* (x.-(a+b)/2)
end

sqrt_laplace_wT0 = xx -> [abs(x) <= 1 ? 0.5*log(4)-Base.MathConstants.eulergamma : 0.5*log(4)-Base.MathConstants.eulergamma-asinh(sqrt(x^2-1)) for x in xx]
sqrt_laplace_U_1 = xx -> [abs(x) <= 1 ? -asin(x) : -sign(x)*pi/2 for x in xx]

end # module
