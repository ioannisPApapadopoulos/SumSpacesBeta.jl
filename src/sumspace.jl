"""
SumSpace{kind,T}()

is a quasi-matrix representing the sum space of the specified kind (1 or 2)
on (-∞,∞).
"""

struct SumSpace{kind,V,T} <: Basis{T} 
    I::V
end
SumSpace{kind, T}(I::Vector{T}) where {kind, T} = SumSpace{kind,Vector{T},T}(I::Vector{T})
SumSpace{kind}(I::Vector{Float64}) where kind = SumSpace{kind,Float64}(I::Vector{Float64})
SumSpace{kind}() where kind = SumSpace{kind}([-1.,1.])

# SumSpaceP{T}(I) where T = SumSpace{1}(I)
# SumSpaceP(I) = SumSpaceP{Float64}(I)
# SumSpaceP() = SumSpaceP([-1.,1.])
# AbstractQuasiVector{T}(::SumSpace{kind}) where {kind,I,T} = SumSpace{kind,I,T}(I)
# AbstractQuasiMatrix{T}(::SumSpace{kind}) where {kind,I,T} = SumSpace{kind,I,T}(I)
 
# function SumSpaceP(I)
#     SumSpace{1}(I)
# end
# SumSpaceP() = SumSpaceP([-1.,1.])

# const SumSpaceP = SumSpace{1}
# const SumSpaceD = SumSpace{2}


# SumSpaceP{V,T}() where {V,T} = SumSpace{1,V,T}()
# SumSpaceP{T}() where T = SumSpaceP{SVector([-1.,1.]), T}()
# SumSpace() = SumSpace{Float64}()

# UltrasphericalArc{V}(h, a, T::TT, U::UU) where {V,TT,UU} = UltrasphericalArc{V,TT,UU}(h,a,T,U)


function SumSpaceP()
    SumSpace{1}()
end

function SumSpaceD()
    SumSpace{2}()
end


axes(S::SumSpace) = (Inclusion(ℝ), _BlockedUnitRange(1:2:∞))


# SumSpace() = SumSpace{1}()

function ==(a::SumSpace{kind}, b::SumSpace{kind}) where kind
    if a.I == b.I 
        return true
    else
        return false
    end
end

==(a::SumSpace, b::SumSpace) = false

function getindex(S::SumSpace{1, V, T}, x::Real, j::Int)::T where {V, T}
    isodd(j) && return ExtendedChebyshevT{T}()[x, (j ÷ 2)+1]
    ExtendedWeightedChebyshevU{T}()[x, j ÷ 2]
end

function getindex(S::SumSpace{2, V, T}, x::Real, j::Int)::T where {V, T}
    isodd(j) && return ExtendedChebyshevU{T}()[x, (j ÷ 2)+1]
    ExtendedWeightedChebyshevT{T}()[x, j ÷ 2]
end

###
# Operators
###

# Credit to Timon Gutleb for the below implementation of the identity mapping
# Identity Sp -> Sd
function \(Sd::SumSpace{2}, Sp::SumSpace{1})
    Sd.I != Sp.I && error("Sum spaces bases not centred on same element")
    T = promote_type(eltype(Sp), eltype(Sd))
    halfvec = mortar(Fill.(1/2,Fill(2,∞)))
    d = Diagonal((-1).^(2:∞))*halfvec
    zs = mortar(Zeros.(Fill(2,∞)))
    ld = Diagonal((-1).^(1:∞))*halfvec
    dat = BlockBroadcastArray(hcat,d,zs,zs,zs,ld,zs)
    dat = BlockVcat([-1.,0.,0.,0.,0.,1.]', dat)
    A = BlockBandedMatrices._BandedBlockBandedMatrix(dat', (axes(dat,1),axes(dat,1)), (2,0), (1,0))
    
    # Timon's version, also gives what I want but with the incorrect block structure
    # halfvec = mortar(Ones.(Fill(2,∞)))
    # d = Diagonal((-1).^(1:∞)./(2 .^[0;ones(∞)]))*halfvec
    # zs = mortar(Zeros.(Fill(2,∞)))
    # ld = Diagonal((-1).^(2:∞)./(2 .^[0;ones(∞)]))*halfvec
    # dat = BlockBroadcastArray(hcat,d,zs,ld)
    # A = BlockBandedMatrices._BandedBlockBandedMatrix(dat', (axes(dat,1),axes(dat,1)), (2,0), (0,0))
    return A
end

# Derivative Sp -> Sd
function *(D::Derivative{<:Real}, Sp::SumSpace{1})
    T = eltype(Sp)
    (a,b) = Sp.I

    fracvec = mortar(Fill.(1. /(b-a),Fill(2,∞)))
    zs = mortar(Zeros.(Fill(2,∞)))
    ld = Diagonal(((-1).^(1:∞)) .* (2 .* ((1:∞) .÷ 2))[2:∞] )*fracvec

    dat = BlockBroadcastArray(hcat,zs,ld)
    dat = BlockVcat(Fill(0,2)', dat)
    A = BlockBandedMatrices._BandedBlockBandedMatrix(dat', (axes(dat,1),axes(dat,1)), (1,0), (0,0))
    return ApplyQuasiMatrix(*, SumSpace{2,T}(Sp.I), A)
end

# Hilbert: Sp -> Sd 
# function *(H::Hilbert{<:Any,<:Any,<:Any}, Sp::SumSpace{1})
#     T = eltype(Sp)
#     halfvec = mortar(Fill.(1/2,Fill(2,∞)))
#     d = Diagonal((-1).^(2:∞))*halfvec
#     zs = mortar(Zeros.(Fill(2,∞)))
#     ld = Diagonal((-1).^(1:∞))*halfvec
#     # dat = BlockBroadcastArray(hcat,onevec,2*onevec,3*onevec,4*onevec,5*onevec,6*onevec,7*onevec,8*onevec,9*onevec)
#     dat = BlockBroadcastArray(hcat,d,zs,ld,zs,zs,zs,ld,zs,d)
#     dat = BlockVcat(Fill(0,9)', dat)
#     A = BlockBandedMatrices._BandedBlockBandedMatrix(dat', (axes(dat,1),axes(dat,1)), (2,0), (1,1))
#     return ApplyQuasiMatrix(*, SumSpace{2,T}(Sp.I), A)
# end

# Hilbert: Sp -> Sp 
function *(H::Hilbert{<:Any,<:Any,<:Any}, Sp::SumSpace{1})
    T = eltype(Sp)
    onevec = mortar(Fill.(1. ,Fill(2,∞)))
    dat = BlockBroadcastArray(hcat,-onevec,zs,onevec)
    dat = BlockVcat(Fill(0,3)', dat)
    A = BlockBandedMatrices._BandedBlockBandedMatrix(dat', (axes(dat,1),axes(dat,1)), (0,0), (1,1))
    return ApplyQuasiMatrix(*, SumSpace{1,T}(Sp.I), A)
end