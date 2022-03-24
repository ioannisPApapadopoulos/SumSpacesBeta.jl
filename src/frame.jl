# Custom SVD solver for least squares problems
function solvesvd(A, b; tol=1e-7)
    U,S,V = svd(A)
    S⁺ = inv.(S)[:]
    S⁺[S⁺ .> 1/tol] .= 0
    c = V * Diagonal(S⁺) * U' * b
    return c
end

# Construct collocation points
function collocation_points(M, Me; endpoints=5, innergap = 0, outergap=1e-5)
    Tp = Float64
    x = Array{Tp}(undef,M+2*Me)
    x[1:M] = LinRange{Tp}(-1+innergap,1-innergap,M)
    x[M+1:M+Me] = LinRange{Tp}(-endpoints,-1-outergap,Me)
    x[M+1+Me:M+2*Me] = LinRange{Tp}(1+outergap,endpoints,Me)
    return x
end

# Convert function evaluation to Riemann sum
function riemann(x, f)
    y = sort(x)
    h = 0.5 .* (
            append!(y[2:end], y[end]) .- y
         .+ y .- append!([y[1]], y[1:end-1])
    )
    b = sqrt.(h).*f(y)
    return b
end

# Just function evaluation
function evaluate(x, f)
    y = sort(x)
    b = f(y)
    return b
end

# Fit low order expansion to higher order expansion
function expansion_sum_space(c, Nn, N, el_no, constant)
    # N1 = N+2; N2 = N+1; Nn1 = Nn+2; Nn2 = Nn+1; 
    v = zeros(1 + el_no*(2*N+2))
    
    if constant == true
        v[1] = c[1]
        for e in 1:el_no
            v[2*(e-1)*N+2*e:2*(e-1)*N+2*e+Nn] = c[2*(e-1)*Nn+2*e:(2*e-1)*Nn+2*e]
            v[(2*e-1)*N+2*e+1:(2*e-1)*N+2*e+1+Nn] = c[(2*e-1)*Nn+2*e+1:2*e*Nn+2*e+1]
        end
    else
        for e in 1:el_no
            v[2*(e-1)*N+2*e:2*(e-1)*N+2*e+Nn] = c[2*(e-1)*Nn+2*e-1:(2*e-1)*Nn+2*e-1]
            v[(2*e-1)*N+2*e+1:(2*e-1)*N+2*e+1+Nn] = c[(2*e-1)*Nn+2*e:2*e*Nn+2*e]
        end
    end
    return v
end

# Construct Least Squares matrix for sum space
function framematrix(x, eT, ewU, Nn, M, Me)
    Tp = eltype(eT)
    eltype(x) == Tp ? x=[x] : x=x

    el = length(x)

    A = Matrix{Tp}(undef, M+2*Me, 2*Nn+3 + (el-1)*(2*Nn+2))
    A .= zero(Tp)
    # Form columns of Least Squares matrix. 

    A[:,1] = riemann(x[1], x -> eT[x,1])


    for els in 1:el
        for iter in 1:Nn+1
            A[:,2*(els-1)*Nn+2*els-1+iter] = riemann(x[els], x -> eT[x,iter+1])
            A[:,(2*els-1)*Nn+2*els+iter] = riemann(x[els], x -> ewU[x,iter])
        end
    end
    return A
end

# Construct Least Squares matrix for dual sum space
function dualframematrix(x, eU, ewT, Nn, M, Me)
    Tp = eltype(eU)
    
    A = Matrix{Tp}(undef, M+2*Me, 2*Nn+7)
    A .= zero(Tp)
    # Form columns of Least Squares matrix. 
    for iter in 1:Nn+4
        A[:,iter] = riemann(x, x -> eU[x,iter])
    end
    for iter in 1:Nn+3
        A[:,Nn+4+iter] = riemann(x, x -> ewT[x,iter])
    end
    return A
end