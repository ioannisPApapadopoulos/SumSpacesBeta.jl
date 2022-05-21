using Revise
using SumSpaces
using LaTeXStrings, Plots
using DelimitedFiles
using LinearAlgebra

N = 61;

# fa = x -> [abs(y) < 1 ? 1 : (abs(y)==1 ? 1/2 : 0) for y in x]
fa = x -> [abs(y) < 1 ? 1 : 0 for y in x]

a = [-5,-3,-1.,1,3,5]; eSp = ElementSumSpace{1}(a); eSd = ElementSumSpace{2}(a);

M = max(N^2,6001)  # Number of collocation points in [-1,1]
Me = M #÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me, a=a, endpoints=[-10,10],innergap=1e-2) # Collocation points
A = framematrix(x, eSp, N, normtype="evaluate") 
f = A[1:end,1:end] \ evaluate(x, fa)

xx = -5:0.011:5;
plot(xx, fa(xx))
ff = eSp[xx,1:length(f)]*f;
plot!(xx, ff)
norm(abs.(fa(xx).-ff), Inf)
wf = findmax(abs.(fa(xx).-eSp[xx,1:length(f)]*f))[2]
# abs.(fa(xx).-eSp[xx,1:length(f)]*f)[401]
xx[wf]


xdual = collocation_points(M, Me, a=a, endpoints=[-10,10], innergap=1e-2)
Adual = dualframematrix(xdual, eSd, N, normtype="evaluate")
fdual = Adual[1:end,1:end] \ evaluate(xdual, fa)
plot(xx, fa(xx))
ff = eSd[xx,1:length(fdual)]*fdual; ff[isnan.(ff)] .= 0;
plot!(xx,ff)
# ff[401]=0.5; ff[601]=0.5;
norm(abs.(fa(xx).-ff), Inf)
findmax(abs.(fa(xx).-ff))[2]
# abs.(fa(xx).-eSp[xx,1:length(f)]*f)[401]
xx[findmax(abs.(fa(xx).-ff))[2]]


p = plot(xx, 
    fa(xx),
    label="Exact",
    xlabel=L"$x$",
    ylabel=L"$y$",
    title="Right-hand side",
    ytickfontsize=10,xlabelfontsize=15,ylabelfontsize=15)
p = plot!(xx, 
    eSp[xx,1:length(f)]*f,
    xlabel=L"$x$",
    ylabel=L"$y$",
    title="Right-hand side",
    label="Approximation",
    ytickfontsize=10,xlabelfontsize=15,ylabelfontsize=15)
savefig(p,"discontinous-rhs1.pdf")


p = plot(xx, 
    abs.(fa(xx) .- eSp[xx,1:length(f)]*f),
    ylabel="Error",
    xlabel="x",
    title="Error plot of right-hand side",
    legend=false)

### Collect errors

errors = []; xx = -5:0.01:5;
Nn = [3,5,7,11,15,21,31,41,51,61,71,81,91,101]
# for N in Nn

#     M = max(N^2,6001)  # Number of collocation points in [-1,1]
#     Me = M #÷ 10  # Number of collocation points in [-2,-1) and (1,2].
#     x = collocation_points(M, Me, a=a, endpoints=[-10,10]) # Collocation points
#     A = framematrix(x, eSp, N, normtype="evaluate") 
#     f = A[1:end,1:end] \ evaluate(x, fa)
#     append!(errors,norm(fa(xx).-eSp[xx,1:length(f)]*f, Inf)) 
# end
for N in Nn

    M = max(N^2,6001)  # Number of collocation points in [-1,1]
    Me = M #÷ 10  # Number of collocation points in [-2,-1) and (1,2].
    x = collocation_points(M, Me, a=a, endpoints=[-10,10],innergap=1e-10) # Collocation points
    A = framematrix(x, eSp, N, normtype="evaluate") 
    f = A[1:end,1:end] \ evaluate(x, fa)
    append!(errors,norm(fa(xx).-eSp[xx,1:length(f)]*f, Inf)) 
end
plot(Nn, log10.(errors), legend=:none, 
    title="Error",
    markers=:circle,
    xlabel=L"$n$",
    ylabel=L"$\log_{10}(\infty\mathrm{-norm \;\; error})$",
    # ylim=[log10(1e-15),log10(1e-3)],
    xtickfontsize=10, ytickfontsize=10,xlabelfontsize=15,ylabelfontsize=15,
    linewidth=2,
    marker=:dot,
    markersize=5)
savefig("errors-rhs-infty.pdf")



 ### Solve a fractional PDE

λ = 1; μ = 0; η = 0# Constants
K = length(a)-1

# Compute support functions
supp = readdlm("uS-lmbda-$λ-mu-$μ-eta-$η/uS-N-$N.txt")
x1 = []; x2 = [];
ywT0 = []; yU_1 = []; ywT1 = []; yU0 = []
append!(ywT0, [supp[3,:]]); append!(yU_1, [supp[4,:]]); append!(ywT1, [supp[5,1:2000210]]); append!(yU0, [supp[6,1:2000210]]); 
x1 = supp[1,:]; x2 = supp[2,1:2000210];
uS = fft_supporter_functions(λ, μ, η, a=a, N=N, W=1e4, δ=1e-2, stabilise=true, correction=true,x1=x1,x2=x2,ywT0=ywT0,yU_1=yU_1,ywT1=ywT1,yU0=yU0); # Actual functions
# Element primal sum space coefficients
cuS = coefficient_supporter_functions(A, x, uS, 2N+3) 


# Create appended sum space
ASp = ElementAppendedSumSpace(uS, cuS, a)

# Create matrix for element 1
Id = (eSd \ ASp)[1:1+K*(2N+6),1:1+K*(2N+6)]
# Id[:,2:4*K] .= 0

x = axes(eSp, 1); H = inv.( x .- x')
# Hm = (1/π).*((eSp\(H*eSp))[1:2N+3,1:2N+3])    # Hilbert: Sp -> Sp
# Cm = [(eSd\(Derivative(x)*eSp)[j])[1:2N+7,1:2N+3] for j in 1:K]# Derivative: Sp -> Sd
Hm = (1/π).*(((H*eSp))[1:2N+3,1:2N+3])    # Hilbert: Sp -> Sp
Cm = [((Derivative(x)*eSp)[j])[1:2N+7,1:2N+3] for j in 1:K]# Derivative: Sp -> Sd
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

plot(xx, eSp[xx,1:length(f)]*f[1:end])
plot!(xx, eSd[xx,1:length(fd)]*fd[1:end])
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

xx = Array(-5.:0.01:5)
yy = ASp[xx,1:length(u)]*u

# append!(errors, norm(dx(xx), Inf))
p = plot(xx, yy, ylabel=L"$y$", xlabel=L"$x$", title="Discontinuous data", ytickfontsize=10,xlabelfontsize=15,ylabelfontsize=15,legend=:topleft, label="Solution")
savefig(p, "example-rhs-jump.pdf")
display(p)