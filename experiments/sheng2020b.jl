using Revise
using SumSpacesBeta, Plots
using SpecialFunctions, HypergeometricFunctions
using LinearAlgebra
using MathLink
using DelimitedFiles
using LaTeXStrings


errors = []
# append!(errors, [[1e-5,31]])
# writedlm("errors-inf.txt",errors)

for N in [3,5,7,11,13,15,21,25,31,41] #, ,31
# for N in [41]
# N = 41  # Truncation degree
    λ = 1; μ = 1; η = 1# Constants

    a = [-5,-3,-1.,1,3,5]
    # a = [6,-2.,2,6]
    K = length(a)-1

    eSp = ElementSumSpace{1}(a)
    eSd = ElementSumSpace{2}(a)

    ua = x ->  exp.(-x.^2)

    # fa = x -> ((λ .- 2η.*x) .* exp.(-x.^2) 
    #             .+ 2/gamma(1/2) * _₁F₁.(1,1/2,-x.^2)
    #             .+ μ * exp.(-x.^2) .* abs.(x) .* erfi.(abs.(x)) ./ x
    # )

    fa = z -> [x==0 ? 2.128379167095513 : ((λ .- 2η.*x) .* exp.(-x.^2) 
    .+ 2/gamma(1/2) * _₁F₁.(1,1/2,-x.^2)
    .+ μ * exp.(-x.^2) .* abs.(x) .* erfi.(abs.(x)) ./ x
    ) for x in z]


    M = max(N^2,6001)  # Number of collocation points in [-1,1]
    Me = M #÷ 10  # Number of collocation points in [-2,-1) and (1,2].
    x = collocation_points(M, Me, a=a, endpoints=[-25,25]) # Collocation points
    A = framematrix(x, eSp, N, normtype="evaluate") 
    f = A[1:end,1:end] \ evaluate(x, fa)
    uc = A[1:end,1:end] \ evaluate(x, ua)

    # xx = -25:0.01:25
    # plot(xx,ua(xx))
    # plot!(xx, eSp[xx,1:length(uc)]*uc)
    # norm(ua(xx).-eSp[xx,1:length(uc)]*uc, Inf)



    # Compute support functions
    supp = readdlm("uS-lmbda-$λ-mu-$μ-eta-$η/uS-N-$N.txt")
    x1 = []; x2 = [];
    ywT0 = []; yU_1 = []; ywT1 = []; yU0 = []
    append!(ywT0, [supp[3,:]]); append!(yU_1, [supp[4,:]]); append!(ywT1, [supp[5,1:2000210]]); append!(yU0, [supp[6,1:2000210]]); 
    x1 = supp[1,:]; x2 = supp[2,1:2000210];
    uS = fft_supporter_functions(λ, μ, η, a=a, N=N, W=1e4, δ=1e-2, stabilise=true, correction=true,x1=x1,x2=x2,ywT0=ywT0,yU_1=yU_1,ywT1=ywT1,yU0=yU0); # Actual functions
    # Element primal sum space coefficients
    cuS = coefficient_supporter_functions(A, x, uS, 2N+3, normtype="evaluate") 


    # Plot sanity check
    xx = -10:0.00001:10
    plot(xx, uS[4][3](xx))
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
    fc = []
    fd = [f[1]' zeros(K*4)' f[2:end]']'
    fd = Id*fd
    uc1 = [uc[1]' zeros(K*4)' uc[2:end]']'
    # Multiply RHS with λ
    fd = BlockArray(fd, vcat(1,Fill(K,(length(fd)-1)÷K)))
    uc1 = BlockArray(uc1, vcat(1,Fill(K,(length(fd)-1)÷K)))
    plot!(xx, eSd[xx,1:length(fd)]*fd)
    # Rearrange coefficients element-wise
    fd = coefficient_stack(fd, N, K, appended=true)
    uc1 = coefficient_stack(uc1, N, K, appended=true)
    for j = 1:K
        # Solve for each element seperately and append to form global
        # vector of coefficients
        append!(u, Dm[j]\fd[Block.(j)])
        append!(fc, Dm[j]*uc1[Block.(j)])
    end
    # Rearrange coefficients back to interlaced
    u = coefficient_interlace(u, N, K, appended=true)
    fc = coefficient_interlace(fc, N, K, appended=true)
    fd = coefficient_interlace(fd[1:end],N, K, appended=true)
    # u[1] = u[1] - ASp[1e2,1:length(u)]'*u[1:end]
    # u[1] = - ASp[0,2:length(u)]'*u[2:end] +ua(0)
    # u[1] = 0

    # p = plot() 
    xx = Array(-5.:0.01:5); #xx = [xx[1:length(xx)÷2]' xx[length(xx)÷2+2:end]']'
    # ylim = [0.,1.]
    append!(errors,[[norm(ua(xx) .- ASp[xx,1:length(u)]*u, Inf),N]])
    writedlm("errors-full-inf.txt", errors)
end

xx = Array(-5.:0.01:5)
yy = ASp[xx,1:length(u)]*u

# append!(errors, norm(dx(xx), Inf))
p = plot(xx,yy, title="Sheng2020", label="Sum space - 5 elements", legend=:topleft)
p = plot!(xx, ua(xx), linewidth=4, ylabel=L"$u(x)$", xlabel=L"$x$", title="Solution", ytickfontsize=10,xlabelfontsize=15,ylabelfontsize=15,legend=false)
p = plot(xx,eSd[xx,1:length(fd)]*fd, title="Right-hand side", label="RHS", legend=false)
p = plot(xx, fa(xx), linewidth=4, ylabel=L"$y$", xlabel=L"$x$", title="Example", ytickfontsize=10,xlabelfontsize=15,ylabelfontsize=15,legend=:topleft, label="Right-hand side")
p = plot!(xx, ua(xx), linewidth=4, ylabel=L"$y$", xlabel=L"$x$", title="Example", ytickfontsize=10,xlabelfontsize=15,ylabelfontsize=15,legend=:topleft, label="Solution")
savefig(p, "example.pdf")
display(p)

p = plot(xx, 
    abs.(ua(xx) .- ASp[xx,1:length(u)]*u), 
    # yaxis=:log,
    ylabel="Error",
    xlabel="x",
    title="Error plot of solution",
    legend=false)

findmax(abs.(ua(xx).-ASp[xx,1:length(u)]*u))[2]
xx[findmax(abs.(ua(xx).-ASp[xx,1:length(u)]*u))[2]]


p = plot(xx, 
    abs.(fa(xx) .- eSp[xx,1:length(f)]*f),
    ylabel="Error",
    xlabel="x",
    title="Error plot of right-hand side",
    legend=false)


p = plot(spy(Dm[1], markersize=2,color=:darktest), 
        xtickfontsize=12, ytickfontsize=12,xlabelfontsize=15,ylabelfontsize=15,
        title= L"$\mathrm{Spy \ plot \ of} \ L^{I_1} \; (n=41)$")
savefig(p, "spy-sheng-1.pdf")
p = plot(spy(Dm[2], markersize=2,color=:darktest), 
        xtickfontsize=12, ytickfontsize=12,xlabelfontsize=15,ylabelfontsize=15,
        title= L"$\mathrm{Spy \ plot \ of} \ L^{I_k}, k \geq 2 \; (n=41)$")
savefig(p, "spy-sheng-2.pdf")