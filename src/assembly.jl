using BandedMatrices: _BandedMatrix
using ClassicalOrthogonalPolynomials: Hcat, Vcat, ∞, ℵ₀, Ones, Zeros

# Identity mapping from appended sum space to dual sum space
function idmap_append2sum(N, yu_2, yu_1, ywt0, ywt1; el_no=1)
    Id = I[1:2*N+3 + (el_no-1)*(2*N+2),1:2*N+3+(el_no-1)*(2*N+2)]

    Id = hcat(Id, yu_2...)
    Id = hcat(Id, yu_1...)
    Id = hcat(Id, ywt0...)
    Id = hcat(Id, ywt1...)
    return Id
end

# Hilbert matrix from sum space to sum space
function hilbertmap(N; el_no=1, Tp=Float64)
    H = Matrix{Tp}(undef,2*N+3+(el_no-1)*(2*N+2),2*N+3+(el_no-1)*(2*N+2))
    H[:,:].= zero(Tp)
    # dc1 = N+2; dc2 = N+1
    H[N+3:2*N+3, 2:N+2] = Matrix{Tp}(-I, N+1, N+1)
    H[2:N+2,N+3:2*N+3] = Matrix{Tp}(I, N+1, N+1)

    for e in 1:el_no
        H[(2*e-1)*N+2*e+1:2*e*N+2*e+1,2*(e-1)*N+2*e:(2*e-1)*N+2*e] = Matrix{Tp}(-I, N+1, N+1)
        H[2*(e-1)*N+2*e:(2*e-1)*N+2*e,(2*e-1)*N+2*e+1:2*e*N+2*e+1] = Matrix{Tp}(I, N+1, N+1)
    end
    return H
end

# Differentiation matrix from sum space to dual sum space
function diffmap(N; a=[-1.,1.], el_no =1, Tp=Float64)
    # dc1 = N+2; dc2 = N+1
    C = Matrix{Tp}(undef, 2*N+7+(el_no-1)*(2*N+6),2*N+3+(el_no-1)*(2*N+2))
    C .= zero(Tp)
    
    el_size = a[2:end] .- a[1:end-1]

    C[3:N+3,2:N+2] = (2. / el_size[1]) .* Diagonal(1.:N+1)
    C[N+6:2*N+6,N+3:2*N+3] = (2. / el_size[1]) .* Diagonal(-1.:-1:-(N+1))
    # C[2*N+7, N+3:2*N+3] = zeros(1,dc2)
    for e in 2:el_no
        C[2*(e-1)*N+6*e-3:(2*e-1)*N+6*e-3,2*(e-1)*N+2*e:(2*e-1)*N+2*e] = (2. / el_size[e]) .* Diagonal(1.:N+1)
        C[(2*e-1)*N+6*e:2*e*N+6*e,(2*e-1)*N+2*e+1:2*e*N+2*e+1] = (2. / el_size[e]) .* Diagonal(-1.:-1:-(N+1))
    end
    return C
end

# Identity map from sum space to dual sum space
function idmap_sum2dual(N; el_no = 1, Tp = Float64)
    # dc1 = N+2; dc2 = N+1
    TU = _BandedMatrix(Vcat(-Ones{Tp}(1,∞)/2,Zeros{Tp}(1,∞),Hcat(Ones{Tp}(1,1),Ones{Tp}(1,∞)/2)), ℵ₀, 0,2)
    L = _BandedMatrix(Vcat(Ones{Tp}(1,∞)/2,Zeros{Tp}(1,∞),-Ones{Tp}(1,∞)/2), ℵ₀, 2,0)
    B = Matrix{Tp}(undef, 2*N+7+(el_no-1)*(2*N+6),2*N+3+(el_no-1)*(2*N+2))
    B .= zero(Tp)
    B[3:N+4,1:N+2] = TU[1:N+2,1:N+2]
    B[2,2] = -one(Tp)/2
    B[1,1] = -one(Tp)
    B[N+5:2*N+7,N+3:2*N+3] = L[1:N+3,1:N+1]   
    for e in 2:el_no
        B[2*(e-1)*N+6*e-4:(2*e-1)*N+6*e-2,2*(e-1)*N+2*e:(2*e-1)*N+2*e] = -L[1:N+3,1:N+1]  
        B[(2*e-1)*N+6*e-1:2*e*N+6*e+1,(2*e-1)*N+2*e+1:2*e*N+2*e+1] = L[1:N+3,1:N+1]  
    end
    return B 
end

# Identity map from appended sum space to dual sum space
function idmap_append2dual(N, yu_2, yu_1, ywt0, ywt1; el_no=1, Tp=Float64)
    Id = idmap_append2sum(N, yu_2, yu_1, ywt0, ywt1, el_no=el_no)
    B = idmap_sum2dual(N, el_no=el_no, Tp=Tp)
    return B * Id
end


function fractionalhelmholtzmap(λ, μ, N; a=[-1.,1.], el_no=1, Tp=Tp)
    H = hilbertmap(N, el_no=el_no, Tp=Tp)
    C = diffmap(N, a=a, el_no=el_no, Tp=Float64)
    B = idmap_sum2dual(N, el_no=el_no, Tp=Tp)
    
    scale = 1.

    D =  λ.*B + μ.*B*H + C*H 
    D = hcat(D, zeros(size(D,1), 4))
    D[1,2*N+4] = scale
    D[2,2*N+5] = scale
    D[N+5,2*N+6] = scale
    D[N+6,2*N+7] = scale

    return D
end