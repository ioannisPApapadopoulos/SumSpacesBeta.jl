function solvesvd(A, b; tol=1e-7)
    U,S,V = svd(A)
    S⁺ = inv.(S)[:]
    S⁺[S⁺ .> 1/tol] .= 0
    c = V * Diagonal(S⁺) * U' * b
    return c
end