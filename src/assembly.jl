using BandedMatrices: _BandedMatrix
using ClassicalOrthogonalPolynomials: Hcat, Vcat, ∞, ℵ₀, Ones, Zeros

# Identity mapping from appended sum space to dual sum space
function idmap_append2sum(N, yu_2, yu_1, ywt0, ywt1)
    Id = I[1:2*N+3,1:2*N+3]
    Id = hcat(Id, yu_2)
    Id = hcat(Id, yu_1)
    Id = hcat(Id, ywt0)
    Id = hcat(Id, ywt1)
    return Id
end

# Hilbert matrix for appended sum space
function hilbertmap(N, Tp)
    H = Matrix{Tp}(undef,2*N+3,2*N+3)
    H[:,:].= zero(Tp)
    dc1 = N+2; dc2 = N+1
    H[dc1+1:dc1+dc2, 2:dc1] = Matrix{Tp}(-I, dc2, dc2)
    H[2:dc1,dc1+1:end] = Matrix{Tp}(I, dc1-1, dc1-1)
    return H
end

# Differentiation matrix for appended sum space
function diffmap(N, Tp)
    dc1 = N+2; dc2 = N+1
    C = Matrix{Tp}(undef, 2*N+7,2*N+3)
    C .= zero(Tp)
    C[3:dc1+1,2:dc1] = Diagonal(1.:dc2)
    C[N+6:end-1,N+3:end] = Diagonal(-1.:-1:-dc2)
    C[end, N+3:end] = zeros(1,dc2)
    return C
end

# Identity map from appended sum space to dual sum space
function idmap_sum2dual(N, Tp)
    dc1 = N+2; dc2 = N+1
    TU = _BandedMatrix(Vcat(-Ones{Tp}(1,∞)/2,Zeros{Tp}(1,∞),Hcat(Ones{Tp}(1,1),Ones{Tp}(1,∞)/2)), ℵ₀, 0,2)
    L = _BandedMatrix(Vcat(Ones{Tp}(1,∞)/2,Zeros{Tp}(1,∞),-Ones{Tp}(1,∞)/2), ℵ₀, 2,0)
    B = Matrix{Tp}(undef, 2*N+7,2*N+3)
    B .= zero(Tp)
    B[3:dc1+2,1:dc1] = TU[1:dc1,1:dc1]
    B[2,2] = -one(Tp)/2
    B[1,1] = -one(Tp)
    B[N+5:end,N+3:end] = L[1:N+3,1:N+1]   
    return B 
end

function idmap_append2dual(N, yu_2, yu_1, ywt0, ywt1, Tp)
    Id = idmap_append2sum(N, yu_2, yu_1, ywt0, ywt1)
    B = idmap_sum2dual(N, Tp)
    return B * Id
end


function fractionalhelmholtzmap(Δt, N, Tp)
    H = hilbertmap(N, Tp)
    C = diffmap(N, Tp)
    B = idmap_sum2dual(N, Tp)
    
    scale = 1.
    D =  Δt * C * H + B
    D = hcat(D, zeros(size(D,1), 4))
    D[1,2*N+4] = scale
    D[2,2*N+5] = scale
    D[N+5,2*N+6] = scale
    D[N+6,2*N+7] = scale

    return D
end