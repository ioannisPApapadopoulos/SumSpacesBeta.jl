"""
SumSpace{kind,T}()

is a quasi-matrix representing the sum space of the specified kind (1 or 2)
on -∞..∞.
"""

struct SumSpace{kind,T} <: Basis{T} end
SumSpace{kind}() where kind = SumSpace{kind,Float64}()

AbstractQuasiArray{T}(::SumSpace{kind}) where {T,kind} = SumSpace{kind,T}()
AbstractQuasiMatrix{T}(::SumSpace{kind}) where {T,kind} = SumSpace{kind,T}()
 
const SumSpaceP = SumSpace{1}
const SumSpaceD = SumSpace{2}

axes(S::SumSpace) = (Inclusion(ℝ), _BlockedUnitRange(1:2:∞))


SumSpace() = SumSpace{1}()

==(a::SumSpace{kind}, b::SumSpace{kind}) where kind = true
==(a::SumSpace, b::SumSpace) = false

function getindex(S::SumSpaceP{T}, x::Real, j::Int)::T where T
    isodd(j) && return ExtendedChebyshevT{T}()[x, (j ÷ 2)+1]
    ExtendedWeightedChebyshevU{T}()[x, j ÷ 2]
end

function getindex(S::SumSpaceD{T}, x::Real, j::Int)::T where T
    isodd(j) && return ExtendedChebyshevU{T}()[x, (j ÷ 2)+1]
    ExtendedWeightedChebyshevT{T}()[x, j ÷ 2]
end

###
# Operators
###

# SumSpaceD() \(D*SumSpaceP())
@simplify function *(D::Derivative{<:Any,<:Real}, S::SumSpaceP)
    T = promote_type(eltype(D),eltype(S))
    A = _BandedMatrix((zero(T):∞)', ℵ₀, -1,1)
    ApplyQuasiMatrix(*, SumSpaceD{T}(), A)
end

# Credit to Timon Gutleb for the below implementation of the identity mapping
function \(Sd::SumSpaceD, Sp::SumSpaceP)
    T = promote_type(eltype(Sp), eltype(Sd))
    halfvec = mortar(Fill.(1/2,Fill(2,∞)))
    d = Diagonal((-1).^(2:∞))*halfvec
    zs = mortar(Zeros.(Fill(2,∞)))
    ld = Diagonal((-1).^(1:∞))*halfvec
    dat = BlockBroadcastArray(hcat,d,zs,zs,zs,ld,zs)
    dat = BlockVcat([-1.,0.,0.,0.,0.,1.]', dat)
    A = BlockBandedMatrices._BandedBlockBandedMatrix(dat', (axes(dat,1),axes(dat,1)), (2,0), (1,0))
    return A
end