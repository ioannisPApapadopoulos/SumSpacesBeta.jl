"""
SumSpace{kind,T}()

is a quasi-matrix representing the sum space of the specified kind (1 or 2)
on (-∞,∞).
"""
###
# Primal and dual sum space
###

struct SumSpace{kind,E,T} <: Basis{T} 
    I::E
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


# SumSpaceP{E,T}() where {E,T} = SumSpace{1,E,T}()
# SumSpaceP{T}() where T = SumSpaceP{SVector([-1.,1.]), T}()
# SumSpace() = SumSpace{Float64}()

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

function getindex(S::SumSpace{1, E, T}, x::Real, j::Int)::T where {E, T}
    y = affinetransform(S.I[1],S.I[2], x)
    isodd(j) && return ExtendedChebyshevT{T}()[y, (j ÷ 2)+1]
    ExtendedWeightedChebyshevU{T}()[y, j ÷ 2]
end

function getindex(S::SumSpace{2, E, T}, x::Real, j::Int)::T where {E, T}
    y = affinetransform(S.I[1],S.I[2], x)
    isodd(j) && return ExtendedChebyshevU{T}()[y, (j ÷ 2)+1]
    ExtendedWeightedChebyshevT{T}()[y, j ÷ 2]
end


###
# Appended sum space
###

struct AppendedSumSpace{AA, CC,E, T} <: Basis{T} 
    A::AA
    C::CC
    I::E
end
AppendedSumSpace{E,T}(A, C, I::Vector{Float64}) where {E,T} = AppendedSumSpace{Any,Any,E,T}(A, C, I::Vector{Float64})
AppendedSumSpace{E}(A, C, I::Vector{Float64}) where E = AppendedSumSpace{E,Float64}(A, C, I::Vector{Float64})
AppendedSumSpace(A, C) = AppendedSumSpace{Vector{Float64}}(A, C, [-1.,1.])


axes(ASp::AppendedSumSpace) = (Inclusion(ℝ), OneToInf())


function getindex(ASp::AppendedSumSpace{AA, CC, E, T}, x::Real, j::Int)::T where {AA, CC, E, T}

    if 1<=j<=4
        return ASp.A[j][1](x)
    else
        return SumSpace{1, E, T}(ASp.I)[x,j-4]
    end
end

### 
# Element primal sum space
###

struct ElementSumSpace{kind,E,T} <: Basis{T} 
    I::E
end
ElementSumSpace{kind, T}(I::Vector{T}) where {kind, T} = ElementSumSpace{kind,Vector{T},T}(I::Vector{T})
ElementSumSpace{kind}(I::Vector{Float64}) where kind = ElementSumSpace{kind,Float64}(I::Vector{Float64})
ElementSumSpace{kind}() where kind = ElementSumSpace{kind}([-1.,1.])

axes(ES::ElementSumSpace) = (Inclusion(ℝ), _BlockedUnitRange(1:(length(ES.I)-1):∞))

function ==(a::ElementSumSpace{kind}, b::ElementSumSpace{kind}) where kind
    if a.I == b.I 
        return true
    else
        return false
    end
end
==(a::ElementSumSpace, b::ElementSumSpace) = false

function getindex(ES::ElementSumSpace{1, E, T}, x::Real, j::Int)::T where {E, T}

    el_no = length(ES.I)-1
    if j == 1
        return SumSpace{1}()[x, 1]
    else
        ind = (j-2) ÷ el_no + 1
        i = isodd(ind) ? (ind ÷ 2)+1 : ind ÷ 2
        el = (j-1) - ((j-2) ÷ el_no)*el_no

        y = affinetransform(ES.I[el],ES.I[el+1], x)
        if isodd(ind)
            return ExtendedWeightedChebyshevU{T}()[y, i]
        else
            return ExtendedChebyshevT{T}()[y, i+1]
        end

    
    end
end

# function getindex(S::SumSpace{2, E, T}, x::Real, j::Int)::T where {E, T}
#     y = affinetransform(S.I[1],S.I[2], x)
#     isodd(j) && return ExtendedChebyshevU{T}()[y, (j ÷ 2)+1]
#     ExtendedWeightedChebyshevT{T}()[y, j ÷ 2]
# end


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

# Identity ASp -> Sd
function \(Sd::SumSpace{2}, ASp::AppendedSumSpace)
    Sd.I != ASp.I && error("Sum spaces bases not centred on same element")
    T = promote_type(eltype(ASp), eltype(Sd))
    
    # This should work but it hangs when attempting BlockHcat
    # Bm = Sd \ SumSpace{1,Vector{T},T}(ASp.I)   
    # Id = Eye((axes(Bm,1),))
    # B = BlockVcat(ASp.C[4][1], mortar(Zeros.(Fill(2,∞))))
    # Id = BlockHcat(B, Id)
    # B = BlockVcat(ASp.C[3][1], mortar(Zeros.(Fill(2,∞))))
    # Id = BlockHcat(B, Id)
    # B = BlockVcat(ASp.C[2][1], mortar(Zeros.(Fill(2,∞))))
    # Id = BlockHcat(B, Id)
    # B = BlockVcat(ASp.C[1][1], mortar(Zeros.(Fill(2,∞))))
    # Id = BlockHcat(B, Id)
    # return Bm * BId


    # FIXME: Temporary hack in finite-dimensional indexing
    N = Int64(1e2)
    Bm = (Sd \ SumSpace{1,Vector{T},T}(ASp.I))[1:2N+7,1:2N+3]    
    B = BlockBroadcastArray(hcat, ASp.C[1][1],ASp.C[2][1],ASp.C[3][1],ASp.C[4][1])[1:end,1:end]
    zs = Zeros(∞,4)
    B = vcat(B, zs)
    Id = hcat(B[1:2N+3,:], I[1:2N+3,1:2N+3])
    A  = Bm * Id

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
    onevec = mortar(Fill.(convert(T, π), Fill(2,∞)))
    zs = mortar(Zeros.(Fill(2,∞)))
    dat = BlockBroadcastArray(hcat,-onevec,zs,onevec)
    dat = BlockVcat(Fill(0,3)', dat)
    A = BlockBandedMatrices._BandedBlockBandedMatrix(dat', (axes(dat,1),axes(dat,1)), (0,0), (1,1))
    return ApplyQuasiMatrix(*, SumSpace{1,T}(Sp.I), A)
end