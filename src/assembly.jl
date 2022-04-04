# Identity mapping from appended sum space to dual sum space
function idmap_append2sum(N, yu0, yu_1, ywt0, ywt1; el_no=1)
    Id = I[1:2*N+3 + (el_no-1)*(2*N+2),1:2*N+3+(el_no-1)*(2*N+2)]

    for e in 1:el_no
        # Id = hcat(Id, yu_1[e])
        # Id = hcat(Id, yu0[e])
        # Id = hcat(Id, ywt0[e])
        # Id = hcat(Id, ywt1[e])
        
        Id = hcat(yu0[e],Id)
        Id = hcat(ywt1[e],Id)
        Id = hcat(yu_1[e],Id)
        Id = hcat(ywt0[e],Id)

    end
    return Id
end

# Hilbert matrix from sum space to sum space
function hilbertmap(N; el_no=1, Tp=Float64)
    H = Matrix{Tp}(undef,1 + el_no*(2*N+2),1 + el_no*(2*N+2))
    H[:,:].= zero(Tp)

    for e in 1:el_no
        H[(2*e-1)*N+2*e+1:2*e*N+2*e+1,2*(e-1)*N+2*e:(2*e-1)*N+2*e] = Matrix{Tp}(-I, N+1, N+1)
        H[2*(e-1)*N+2*e:(2*e-1)*N+2*e,(2*e-1)*N+2*e+1:2*e*N+2*e+1] = Matrix{Tp}(I, N+1, N+1)
    end
    return H
end

# Differentiation matrix from sum space to dual sum space
function diffmap(N; a=[-1.,1.], el_no =1, Tp=Float64)
    # dc1 = N+2; dc2 = N+1
    C = Matrix{Tp}(undef, 1 + el_no*(2*N+6), 1 + el_no*(2*N+2))
    C .= zero(Tp)
    
    el_size = a[2:end] .- a[1:end-1]

    for e in 1:el_no
        C[2*(e-1)*N+6*e-3:(2*e-1)*N+6*e-3,2*(e-1)*N+2*e:(2*e-1)*N+2*e] = (2. / el_size[e]) .* Diagonal(1.:N+1)
        C[(2*e-1)*N+6*e:2*e*N+6*e,(2*e-1)*N+2*e+1:2*e*N+2*e+1] = (2. / el_size[e]) .* Diagonal(-1.:-1:-(N+1))
    end
    return C
end

# Identity map from sum space to dual sum space
function idmap_sum2dual(N; el_no = 1, Tp = Float64)
    # dc1 = N+2; dc2 = N+1
    # TU = _BandedMatrix(Vcat(-Ones{Tp}(1,∞)/2,Zeros{Tp}(1,∞),Hcat(Ones{Tp}(1,1),Ones{Tp}(1,∞)/2)), ℵ₀, 0,2)
    L = _BandedMatrix(Vcat(Ones{Tp}(1,∞)/2,Zeros{Tp}(1,∞),-Ones{Tp}(1,∞)/2), ℵ₀, 2,0)
    B = Matrix{Tp}(undef, 1+el_no*(2*N+6),1+el_no*(2*N+2))
    B .= zero(Tp)
    # B[3:N+4,1:N+2] = TU[1:N+2,1:N+2]
    # B[2,2] = -one(Tp)/2
    # B[1,1] = -one(Tp)
    # B[N+5:2*N+7,N+3:2*N+3] = L[1:N+3,1:N+1]   
    for e in 1:el_no
        B[2*(e-1)*N+6*e-4:(2*e-1)*N+6*e-2,2*(e-1)*N+2*e:(2*e-1)*N+2*e] = -L[1:N+3,1:N+1]  
        B[(2*e-1)*N+6*e-1:2*e*N+6*e+1,(2*e-1)*N+2*e+1:2*e*N+2*e+1] = L[1:N+3,1:N+1]  
    end
    B[1,1] = -one(Tp)
    B[3,1] = one(Tp)
    return B 
end

# Identity map from appended sum space to dual sum space
function idmap_append2dual(N, yu0, yu_1, ywt0, ywt1; el_no=1, Tp=Float64)
    Id = idmap_append2sum(N, yu0, yu_1, ywt0, ywt1, el_no=el_no)
    B = idmap_sum2dual(N, el_no=el_no, Tp=Tp)
    return B * Id
end


function fractionalhelmholtzmap(λ, μ, N; a=[-1.,1.],Tp=Float64)
    el_no = length(a)-1

    H = hilbertmap(N, el_no=el_no, Tp=Tp)
    C = diffmap(N, a=a, el_no=el_no, Tp=Float64)
    B = idmap_sum2dual(N, el_no=el_no, Tp=Tp)
    
    scale = 1.

    D =  λ.*B + μ.*B*H + C*H 
    D = hcat(D, zeros(size(D,1), 4*el_no))
    D[2,end-4*el_no+1] = scale
    D[3,end-4*el_no+2] = scale
    D[N+5,end-4*el_no+3] = scale
    D[N+6,end-4*el_no+4] = scale

    for e in 2:el_no
        D[2*(e-1)*N+6*(e-1)+2,end-4*(el_no-e+1)+1] = scale
        D[2*(e-1)*N+6*(e-1)+3,end-4*(el_no-e+1)+2] = scale
        D[(2*e-1)*N+6*e-1,end-4*(el_no-e+1)+3] = scale
        D[(2*e-1)*N+6*e,end-4*(el_no-e+1)+4] = scale
    end

    return D
end

function fractionallaplacemap(N; a=[-1.,1.], Tp=Float64)

    length(a) > 2 && error("Multiple elements not implemented yet")

    D = fractionalhelmholtzmap(zero(Tp), zero(Tp), N, a=a, Tp=Tp)
    D = [D[1:N+3, :]' D[N+5:end,:]']'
    D = D[1:end-1, :]
    D = D[2:end, :]
    D = [D[:, 1:2N+4] D[:, 2N+6]]
    D = D[:,2:end]
    
    return D
end



function split_block_helmholtz_matrix(D, el_no)
    block_length = size(D)[1] ÷ el_no
    Dt = D[1:block_length+1,1:block_length-3]
    Dt = hcat(Dt, D[1:block_length+1, end-el_no*4+1])
    Dt = hcat(Dt, D[1:block_length+1, end-el_no*4+2])
    Dt = hcat(Dt, D[1:block_length+1, end-el_no*4+3])
    Dt = hcat(Dt, D[1:block_length+1, end-el_no*4+4])

    Ds = [Dt]


    for e in 2 : el_no
        Dt = D[(e-1)*block_length+2:e*block_length+1,(e-1)*block_length-4*e+6:e*block_length-4*e+1]
        Dt = hcat(Dt, D[(e-1)*block_length+2:e*block_length+1, end-(el_no-e+1)*4+1])
        Dt = hcat(Dt, D[(e-1)*block_length+2:e*block_length+1, end-(el_no-e+1)*4+2])
        Dt = hcat(Dt, D[(e-1)*block_length+2:e*block_length+1, end-(el_no-e+1)*4+3])
        Dt = hcat(Dt, D[(e-1)*block_length+2:e*block_length+1, end-(el_no-e+1)*4+4])

        append!(Ds,  [Dt])
    end
    return Ds
end

function split_block_helmholtz_vector(v, el_no)
    block_length = length(v) ÷ el_no
    
    vs = [v[1:block_length+1]]
    for e in 2 : el_no
        append!(vs,  [v[(e-1)*block_length+2:e*block_length+1]])
    end
    return vs
end