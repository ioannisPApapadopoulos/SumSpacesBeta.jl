using Revise
using SumSpaces, LinearAlgebra
using Plots

N = 7;
μ = 0; η = 0;
fxλ = (x, λ) -> (sqrt(π)*exp(λ/4)).*ExtendedWeightedChebyshevU()[x,5]

# for λ in [-1.0]
λ = -1.
    a = [-5,-3,-1.,1,3,5]; K = length(a)-1

    eSp = ElementSumSpace{1}(a)
    eSd = ElementSumSpace{2}(a)

    fa = x -> fxλ(x, λ)

    M = max(N^2,6001)  # Number of collocation points in [-1,1]
    Me = M #÷ 10  # Number of collocation points in [-2,-1) and (1,2].
    x = collocation_points(M, Me, a=a, endpoints=[-25,25]) # Collocation points
    A = framematrix(x, eSp, N, normtype="evaluate") 
    f = A[1:end,1:end] \ evaluate(x, fa)

    uS = fft_supporter_functions(λ, μ, η, a=a, N=N, W=1e4, δ=1e-2, stabilise=true, correction=false);
    cuS = coefficient_supporter_functions(A, x, uS, 2N+3, normtype="evaluate") 

    # Plot sanity check
    xx = -10:0.001:10
    plot(xx, uS[1][3](xx))
    y = eSp[xx,1:length(cuS[1][1])]*cuS[4][3]
    plot!(xx, y)

    # Create appended sum space
    ASp = ElementAppendedSumSpace(uS, cuS, a)

    # Create matrix for element 1
    Id = (eSd \ ASp)[1:1+K*(2N+6),1:1+K*(2N+6)]

    x = axes(eSp, 1); H = inv.( x .- x')
    Hm = (1/π).*((eSp\(H*eSp))[1:2N+3,1:2N+3])    # Hilbert: Sp -> Sp
    Cm = [(eSd\(Derivative(x)*eSp)[j])[1:2N+7,1:2N+3] for j in 1:K]# Derivative: Sp -> Sd
    Bm = (eSd\eSp)[1:2N+7,1:2N+3]                 # Identity: Sp -> Sd


    Dm =  [λ.*Bm + μ.*Bm*Hm + η.*Cm[j] + Cm[j]*Hm for j in 1:K]     # Helmholtz-like operator: Sp -> Sd   
    Dm = [hcat(zeros(size(Dm[j],1), 4),Dm[j]) for j in 1:K] # Adding 4 columns to construct: ASp -> Sd
    for j in 1:K
        Dm[j][2:3,1:2] = I[1:2,1:2]; Dm[j][end-1:end,3:4] = I[1:2,1:2]
        if j == 1
            # In first element permute T0 column to start
            Dm[j] = [Dm[j][:,5] Dm[j][:,1:4] Dm[j][:,6:end]] 
        else
            # In the rest delete the T0 column and row
            Dm[j] = [Dm[j][:,1:4] Dm[j][:,6:end]]
            Dm[j] = Dm[j][2:end,:] 
        end
    end
    
    u = []
    fd = [f[1]' zeros(K*4)' f[2:end]']'
    fd = Id*fd
    # Multiply RHS with λ
    fd = BlockArray(fd, vcat(1,Fill(K,(length(fd)-1)÷K)))
    # Rearrange coefficients element-wise
    fd = coefficient_stack(fd, N, K, appended=true)
    for j = 1:K
        # Solve for each element seperately and append to form global
        # vector of coefficients
        append!(u, Dm[j]\fd[Block.(j)])
    end
    # Rearrange coefficients back to interlaced
    u = coefficient_interlace(u, N, K, appended=true)
    fd = coefficient_interlace(fd[1:end],N, K, appended=true)
    
    xx = Array(-10.:0.01:10)
    yy = ASp[xx,1:length(u)]*u
    p = plot(xx,yy, title="Wave Propogation, λ=$λ", 
            label="Sum space - 5 elements", 
            legend=:topleft)

# end