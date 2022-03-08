# Custom SVD solver for least squares problems
function solvesvd(A, b; tol=1e-7)
    U,S,V = svd(A)
    S⁺ = inv.(S)[:]
    S⁺[S⁺ .> 1/tol] .= 0
    c = V * Diagonal(S⁺) * U' * b
    return c
end

# Construct collocation points
function collocation_points(M, Me; endpoints=5, gap=1e-5)
    Tp = Float64
    x = Array{Tp}(undef,M+2*Me)
    x[1:M] = LinRange{Tp}(-1,1,M)
    x[M+1:M+Me] = LinRange{Tp}(-endpoints,-1-gap,Me)
    x[M+1+Me:M+2*Me] = LinRange{Tp}(1+gap,endpoints,Me)
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
function expansion(N, Nn, c)
    v = zeros(2*N+3)
    v[1:Nn+2] = c[1:Nn+2]
    v[N+3:N+3+Nn] = c[Nn+3:end]
    return v
end

# Construct Least Squares matrix
function framematrix(x, eT, ewU, Nn, M, Me)
    Tp = eltype(eT)
    
    A = Matrix{Tp}(undef, M+2*Me, 2*Nn+3)
    A .= zero(Tp)
    # Form columns of Least Squares matrix. 
    for iter in 1:Nn+2
        A[:,iter] = riemann(x, x -> eT[x,iter])
    end
    for iter in 1:Nn+1
        A[:,Nn+2+iter] = riemann(x, x -> ewU[x,iter])
    end
    return A
end