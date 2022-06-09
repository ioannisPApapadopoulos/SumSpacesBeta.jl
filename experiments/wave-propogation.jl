using Revise
using SumSpaces, LinearAlgebra, Interpolations
using Plots

N = 7;
μ = 0; η = 0;
fxλ = (x, λ) -> (sqrt(π)*exp(λ/4)).*ExtendedWeightedChebyshevU()[x,5]
solns = []
# for λ in [-1.0]
# for λ in [-1,-10,-35]
λ = -2e1
    a = [-5,-3,-1.,1,3,5]; K = length(a)-1

    eSp = ElementSumSpace{1}(a)
    eSd = ElementSumSpace{2}(a)

    fa = x -> fxλ(x, λ)

    M = max(N^2,6001)  # Number of collocation points in [-1,1]
    Me = M #÷ 10  # Number of collocation points in [-2,-1) and (1,2].
    x = collocation_points(M, Me, a=a, endpoints=[-25,25]) # Collocation points
    A = framematrix(x, eSp, N, normtype="evaluate") 
    f = A[1:end,1:end] \ evaluate(x, fa)

    λ1 = λ + 1e2*eps()*im
    if abs(λ) < 1 
        (xx, s, ywT0, yU_1, ywT1, yU0) = supporter_functions(λ1, μ, η, W=1e2*abs(λ), δ=1e-3*abs(λ), a=a, N=N, stabilise=true);
    else
        (xx, s, ywT0, yU_1, ywT1, yU0) = supporter_functions(λ1, μ, η, W=1e2, δ=1e-3, a=a, N=N, stabilise=true);
    end
    f1 = interpolate((xx,), real.(ywT0[1])[:], Gridded(Linear()))
    f2 = interpolate((xx,), real.(yU_1[1])[:], Gridded(Linear()))
    f3 = interpolate((xx,), real.(ywT1[1])[:], Gridded(Linear()))
    f4 = interpolate((xx,), real.(yU0[1])[:], Gridded(Linear()))
    # c = f1(-97π/(2λ)); ywT0[1] = ywT0[1].-c*sign.(λ.*xx).*sin.(λ.*xx)
    # c = f2(-100π/(λ)); yU_1[1] = yU_1[1].+c*sign.(λ.*xx).*cos.(λ.*xx)

    c = f1(-13π/(2λ)); ywT0[1] = ywT0[1].-c*sign.(λ.*xx).*sin.(λ.*xx)
    c = f2(-100π/(λ)); yU_1[1] = yU_1[1].+c*sign.(λ.*xx).*cos.(λ.*xx)
    c = f3(-100π/(-λ)); ywT1[1] = ywT1[1].+c*sign.(-λ.*xx).*cos.(-λ.*xx)
    # c = f4(-97π/(-2λ)); yU0[1] = yU0[1].-c*sign.(-λ.*xx).*sin.(λ.*xx)
    c = f4(-97π/(-2λ)); yU0[1] = yU0[1].-c*sign.(-λ.*xx).*sin.(-λ.*xx)
    uS = interpolate_supporter_functions(xx, xx, [1.0], ywT0, yU_1, ywT1, yU0, a=a)

    # uS = fft_supporter_functions(λ1, μ, η, a=a, N=N, W=1e2, δ=1e-3, stabilise=true, correction=false);
    # cuS = coefficient_supporter_functions(A, x, uS, 2N+3, normtype="evaluate") 

    # Plot sanity check
    xx = -20:0.01:20
    plot(xx, uS[1][3](xx)) 
    plot(xx, uS[2][3](xx))  
    plot(xx, uS[3][3](xx)) 
    plot(xx, uS[4][3](xx))    

    # d =  (a[1:end-1] + a[2:end]) ./ 2
    # # for k = 1:K
    # j = 1
    # # c = findmin(uS[j][3](xx))[1]
    # c = uS[j][3](-133π/20)
    # plot(xx, uS[j][3](xx))#.-c*sign.(xx).*sin.(xx.*10))
    
    # j = 2
    # c = findmin(uS[j][3](xx))[1]
    # c = uS[j][3](-6π)
    # plot(xx, uS[j][3](xx).+c*sign.(xx).*cos.(xx))    
    # y = eSp[xx,1:length(cuS[1][1])]*cuS[4][3]
    # plot!(xx, y)

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
    
    append!(solns, [u])

    xx = Array(-10.:0.01:10)
    yy = ASp[xx,1:length(u)]*u
    p = plot(xx,yy, title="Wave Propogation, λ=$λ", 
            label="Sum space - 5 elements", 
            legend=:topleft)

# end  