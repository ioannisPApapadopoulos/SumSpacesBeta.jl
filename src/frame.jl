# Custom SVD solver for least squares problems
function solvesvd(A, b; tol=1e-7, block=true)
    
    Am = A[1:end,1:end] # BlockBandedMatrices cannot do SVD

    U,S,V = svd(Am)
    S⁺ = inv.(S)[:]
    S⁺[S⁺ .> 1/tol] .= 0
    c = V * Diagonal(S⁺) * U' * b

    # Rearrange into block structure
    if block == true
        # FIXME: Get correct block structure in output vector
        el_no = length(A[1, Block.(2)])
        if el_no == 1
            c = BlockArray(c, vcat(1,Fill(2,(length(c)-1)÷2)))
        else
            c = BlockArray(c, vcat(1,Fill(el_no,(length(c)-1)÷el_no)))
        end
    end
    return c
end

# Construct collocation points
function collocation_points(M, Me; endpoints=5, innergap = 0, outergap=1e-5)
    Tp = Float64
    x = Array{Tp}(undef,M+2*Me)
    x[1:M] = LinRange{Tp}(-1+innergap,1-innergap,M)
    # x[1:M] = LinRange{Tp}(0+innergap,1-innergap,M)
    # x[1:M] = cos.(π.*x[1:M])
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
function expansion_sum_space(c, N, el_no)
    v = zeros(1 + el_no*(N-1))
    if el_no == 1
        v = BlockArray(v, vcat(1,Fill(2,(length(v)-1)÷2)))
    else
        v = BlockArray(v, vcat(1,Fill(el_no,(length(v)-1)÷el_no)))
    end
    
    v[1:length(c)] = c
    return v
end

# Construct Least Squares matrix for sum space
function framematrix(x, Sp, Nn, M, Me)
    Tp = eltype(Sp)
    el = length(Sp.I) - 1
    if typeof(Sp) == SumSpace{1, Vector{Tp}, Tp}
        rows = [M+2*Me]; cols = vcat([1], Fill(2, el*(Nn+1)))
    elseif typeof(Sp) == ElementSumSpace{1, Vector{Tp}, Tp}
        rows = [M+2*Me]; cols = vcat([1], Fill(el, (2*Nn+2)))
    else
        error("Use either SumSpace{1} or ElementSumSpace{1}")
    end

    # Create correct block structure
    A = BlockBandedMatrix(Zeros(sum(rows),sum(cols)), rows, cols, (sum(rows),sum(cols)))
    
    # Form columns of Least Squares matrix.
    A[:,Block.(1:length(cols))] = riemann(x, x->Sp[x, Block.(1:length(cols))])
    return A
end

# Construct Least Squares matrix for dual sum space
function dualframematrix(x, Sd, Nn, M, Me)
    Tp = eltype(Sd)
    el = length(Sd.I) - 1
    if typeof(Sd) == SumSpace{2, Vector{Tp}, Tp}
        rows = [M+2*Me]; cols = vcat([1], Fill(2, el*(Nn+3)))
    elseif typeof(Sd) == ElementSumSpace{2, Vector{Tp}, Tp}
        rows = [M+2*Me]; cols = vcat([1], Fill(el, (2*Nn+6)))
    else
        error("Use either SumSpace{2} or ElementSumSpace{2}")
    end

    # Create correct block structure
    A = BlockBandedMatrix(Zeros(sum(rows),sum(cols)), rows, cols, (sum(rows),sum(cols)))
    
    # Form columns of Least Squares matrix.
    A[:,Block.(1:length(cols))] = riemann(x, x->Sd[x, Block.(1:length(cols))])
    return A
end