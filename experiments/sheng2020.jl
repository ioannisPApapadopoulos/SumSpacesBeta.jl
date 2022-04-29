using Revise
using SumSpaces, Plots
using SpecialFunctions, HypergeometricFunctions
using LinearAlgebra

N = 31  # Truncation degree
λ = 1; μ = 0; η = 0# Constants

a = [-5,-3,-1.,1,3,5]
K = length(a)-1

eSp = ElementSumSpace{1}(a)
eSd = ElementSumSpace{2}(a)

ua = x ->  exp.(-x.^2)

fa = x -> 2/gamma(1/2) * _₁F₁.(1,1/2,-x.^2) .+ λ * exp.(-x.^2)

M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me, a=a, endpoints=[-5*K,5*K]) # Collocation points
A = framematrix(x, eSp, N) 
f = A[1:end,1:end] \ riemann(x, fa)


xx = -25:0.01:25
plot(xx,fa(xx))
plot!(xx, eSp[xx,1:length(f)]*f)
norm(fa(xx).-eSp[xx,1:length(f)]*f, Inf)



# Compute support functions
uS = fft_supporter_functions(λ, μ, η, a=a, N=N, stabilise=true) # Actual functions
# Element primal sum space coefficients
cuS = coefficient_supporter_functions(A, x, uS, 2N+3) 


# Plot sanity check
xx = -10:0.01:10
plot(xx, uS[4][3](xx))
y = eSp[xx,1:length(cuS[1][1])]*cuS[1][1]
plot!(xx, y)

# Create appended sum space
ASp = ElementAppendedSumSpace(uS, cuS, a)

# Create matrix for element 1
Id = (eSd \ ASp)[1:1+K*(2N+6),1:1+K*(2N+6)]

x = axes(eSp, 1); H = inv.( x .- x')
Hm = (1/π).*((eSp\(H*eSp))[1:2N+3,1:2N+3])    # Hilbert: Sp -> Sp
Cm = [(eSd\(Derivative(x)*eSp)[j])[1:2N+7,1:2N+3] for j in 1:K]# Derivative: Sp -> Sd
Bm = (eSd\eSp)[1:2N+7,1:2N+3]                 # Identity: Sp -> Sd


Dm =  [λ.*Bm + μ.*Bm*Hm + Cm[j]*Hm for j in 1:K]     # Helmholtz-like operator: Sp -> Sd   
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
plot!(xx, eSd[xx,1:length(fd)]*fd)
# Rearrange coefficients element-wise
fd = coefficient_stack(fd, N, K, appended=true)
for j = 1:K
    # Solve for each element seperately and append to form global
    # vector of coefficients
    append!(u, Dm[j]\fd[Block.(j)])
end
# Rearrange coefficients back to interlaced
u = coefficient_interlace(u, N, K, appended=true)
u[1] = u[1] - ASp[1e2,1:length(u)]'*u

uc = x -> ASp[x,1:length(u)]'*u
p = plot() 
xx = -5.:0.01:5.
ylim = [0.,1.]
errors = []
norm(ua(xx).-ASp[xx,1:length(u)]*u, Inf)
findmax(ua(xx).-ASp[xx,1:length(u)]*u)[2]
xx[2]

xx = -20:0.01:20
yy = ASp[xx,1:length(u)]*u
# append!(errors, norm(dx(xx), Inf))
p = plot(xx,yy, title="Sheng2020", label="Sum space - 5 elements", legend=:topleft, ylim=ylim)
p = plot!(xx, ua(xx), label="True solution", legend=:topleft)
display(p)

plot(xx, abs.(ua(xx) .- ASp[xx,1:length(u)]*u), yaxis=:log)
p = plot()
p = plot(xx, abs.(fa(xx) .- eSp[xx,1:length(f)]*f))