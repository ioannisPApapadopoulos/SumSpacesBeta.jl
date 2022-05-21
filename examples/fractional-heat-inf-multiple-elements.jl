using Revise
using SumSpaces, Plots
using QuadGK
import LinearAlgebra: I, norm
using LaTeXStrings, DelimitedFiles
using PyPlots

"""
Solve the fractional heat equation with 3 elements at [-3,-1] ∪ [-1,1] ∪ [1,3].
"""


N = 5 # Truncation degree
λ = 1e2; μ = 0; η = 0; Δt = 1/λ # Constants

a = [-5,-3,-1,1.,3,5] # 3 elements at [-3,-1] ∪ [-1,1] ∪ [1,3]
el_no = length(a)-1

eSp = ElementSumSpace{1}(a) # Primal element sum space
eSd = ElementSumSpace{2}(a) # Dual element sum space

M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10 + 1  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me, a=a, endpoints=[-20.,20]) # Collocation points

Nn = N#min(N,7) # Truncation degree to approximate the supporter functions

A = framematrix(x, eSp, Nn) # Blocked frame matrix

# Compute support functions
# @time uS = fft_supporter_functions(λ, μ, η, a=a); # Actual functions
supp = readdlm("uS-lmbda-$λ-mu-$μ-eta-$η/uS-N-1.txt")
x1 = []; x2 = [];
ywT0 = []; yU_1 = []; ywT1 = []; yU0 = []
append!(ywT0, [supp[3,:]]); append!(yU_1, [supp[4,:]]); append!(ywT1, [supp[5,1:2000210]]); append!(yU0, [supp[6,1:2000210]]); 
x1 = supp[1,:]; x2 = supp[2,1:2000210];
uS = fft_supporter_functions(λ, μ, η, a=a, N=N, W=1e4, δ=1e-2, correction=true,x1=x1,x2=x2,ywT0=ywT0,yU_1=yU_1,ywT1=ywT1,yU0=yU0); # Actual functions

# Element primal sum space coefficients
@time cuS = coefficient_supporter_functions(A, x, uS, 2N+3); 

# Plot sanity check
xx = -10:0.01:10
plot(xx, uS[1][1](xx))
y = eSp[xx,1:length(cuS[1][1])]*cuS[1][1]
plot!(xx, y)

# Create appended sum space
ASp = ElementAppendedSumSpace(uS, cuS, a)

# Create matrix for element 1
Id = (eSd \ ASp)[1:2N+7+(el_no-1)*(2N+6),1:2N+7+(el_no-1)*(2N+6)]

x = axes(eSp, 1); H = inv.( x .- x')
Hm = (1/π).*((eSp\(H*eSp))[1:2N+3,1:2N+3])    # Hilbert: Sp -> Sp
Cm = [(eSd\(Derivative(x)*eSp)[j])[1:2N+7,1:2N+3] for j in 1:el_no]# Derivative: Sp -> Sd
Bm = (eSd\eSp)[1:2N+7,1:2N+3]                 # Identity: Sp -> Sd


Dm =  [λ.*Bm + μ.*Bm*Hm + Cm[j]*Hm for j in 1:el_no]     # Helmholtz-like operator: Sp -> Sd   
Dm = [hcat(zeros(size(Dm[j],1), 4),Dm[j]) for j in 1:el_no] # Adding 4 columns to construct: ASp -> Sd
for j in 1:el_no
    Dm[j][2:5,1:4] = I[1:4,1:4]
    if j == 1
        # In first element permute T0 column to start
        Dm[j] = [Dm[j][:,5] Dm[j][:,1:4] Dm[j][:,6:end]] 
    else
        # In the rest delete the T0 column and row
        Dm[j] = [Dm[j][:,1:4] Dm[j][:,6:end]]
        Dm[j] = Dm[j][2:end,:] 
    end
end

# Initial condition, u₀ = √(1-x²)
u₀ = zeros(1+el_no*(2N+6))
u₀ = BlockArray(u₀, vcat(1,Fill(el_no,(length(u₀)-1)÷el_no)))
u₀[Block.(6)][3] = 1.
u = [u₀]

# u0 = x -> 1. ./ ((x.^2 .+ 1) )
# u₀ = zeros(2*N+7)
# u₀[6] = 1.
# x = collocation_points(M, Me, a=a, endpoints=[-20.,20])
# A = framematrix(x, eSp, N, norm="riemann")
# u₀ = solvesvd(A, riemann(x, u0))#, tol=1e-3)
# u₀₀= zeros(1+el_no*(2N+6)); u₀₀[1]=u₀[1]; u₀₀[2+4*el_no:4*el_no+length(u₀)]=u₀[2:end]
# u = [u₀₀]

# Run solve loop for time-stepping
timesteps=100
@time for k = 1:timesteps
    u1 = []
    
    # Map from ASp to Sd
    v = Id * u[k]
    # Multiply RHS with λ
    v = λ.*v
    v = BlockArray(v, vcat(1,Fill(el_no,(length(v)-1)÷el_no)))
    
    # Rearrange coefficients element-wise
    v = coefficient_stack(v, N, el_no, appended=true)

    for j = 1:el_no
        # Solve for each element seperately and append to form global
        # vector of coefficients
        append!(u1, Dm[j]\v[Block.(j)])
     end

    
     # Rearrange coefficients back to interlaced
    u1 = coefficient_interlace(u1, N, el_no, appended=true)
    # u1[1] = u1[1] - ASp[2e2,1:length(u1)]'*u1

    # if mod(k,5) == 0
    #     y = x->ASp[x,1:length(u1)]*u1
    #     u1 = solvesvd(A, riemann(x, y), tol=1e-5)
    #     u1 = [u1[1]' zeros(4*el_no)' u1[2:end]']'
    # end


    # Append solution to list for different time-steps
    append!(u,  [u1])

    # Normalise constant so that it decays to zero
    u[k+1][1] = u[k+1][1] - ASp[2e2,1:length(u[k+1])]'*u[k+1]

end


# Plot solution
p = plot() 
xx = fx[-20 .< fx .< 20]
xlim = [xx[1],xx[end]]; ylim = [-0.02,1]
y = (x,t) -> (1 + t) ./ ((x.^2 .+ (1+t).^2))
d = (x,t,u) -> abs.(y(x,t) .- ASp[x,1:length(u)]*u)
errors = []

# anim = @animate  for k = 2: timesteps+1
for k = 2:timesteps+1
    t = Δt*(k-1)
    
    tdisplay = round(t, digits=2)
    yy = ASp[xx,1:length(u[k])]*u[k]
    
    dx = x->d(x,t,u[k])
    # append!(errors, sqrt(quadgk(dx, -5, 5)[1]))
    append!(errors, norm(real.(fv[k-1][-20 .< fx .< 20])-ASp[xx,1:length(u[k])]*u[k], Inf))

    # p = plot(xx,yy, title="time=$tdisplay (s)", label="Sum space - 5 elements", legend=:topleft, xlim=xlim, ylim=ylim)
    # p = plot!(xx, real.(fv[k-1][-20 .< fx .< 20]), label="Fourier solution", legend=:topleft, xlim=xlim, ylim=ylim)
    # p = plot!(xx, y(xx, t), label="True solution", legend=:topleft, xlim=xlim, ylim=ylim)
    # sleep(0.001)
    # display(p)
end
# gif(anim, "anim_fps10.gif", fps = 10)

# plot(2:length(errors)+1, errors5)
p = plot(2:length(errors)+1, errors, legend=:none, 
    title="Error",
    # yaxis=:log,
    markers=:circle,
    xlabel=L"$k$",
    xtickfontsize=12, ytickfontsize=12,xlabelfontsize=15,ylabelfontsize=15,
    ylabel=L"$\mathrm{FFT}-\mathbf{S}^{\mathbf{I},\!\!\!\!+}_5\!\!\!\!\!(x) \mathbf{u}_k)\Vert_\infty$")
savefig(p, "errors-infty-W0.pdf")
 
# xx = -10:0.01:10
# xlim = [xx[1],xx[end]]; ylim = [-0.02,1]
# p = plot()
# for k = [1,51,101]
#     t = Δt*(k-1)
    
#     tdisplay = round(t, digits=2)
#     yy = ASp[xx,1:length(u[k])]*u[k]
    
#     if t ≈ 0 || t ≈ 0.5 || t ≈ 1
#         p = plot!(xx,yy, title=L"$\mathrm{5\ elements}$", 
#                 label=L"$\mathrm{time}=$"*"$tdisplay"*L"$\ \mathrm{(s)}$", 
#                 legendfontsize = 10, legend=:topleft, xlim=xlim, ylim=ylim,
#                 xlabel=L"$x$",
#                 xtickfontsize=12, ytickfontsize=12,xlabelfontsize=15,ylabelfontsize=15,
#                 ylabel=L"$\mathbf{S}^{\mathbf{I},\!\!\!\!+}_5\!\!\!\!\!(x) \mathbf{u}$")
#     end
#     # p = plot!(xx, y(xx, t), label="True solution", legend=:topleft, xlim=xlim, ylim=ylim)
#     # sleep(0.001)
#     display(p)
# end  
# savefig(p, "ic1.pdf")


# p = plot(spy(Dm[1], markersize=4,color=:darktest), 
#         xtickfontsize=12, ytickfontsize=12,xlabelfontsize=15,ylabelfontsize=15,
#         title= L"$\mathrm{Spy \ plot \ of} \ \lambda E + A^{I_1} \; (N=11)$")
# savefig(p, "spy-1.pdf")
# p = plot(spy(Dm[2], markersize=4,color=:darktest), 
#         xtickfontsize=12, ytickfontsize=12,xlabelfontsize=15,ylabelfontsize=15,
#         title= L"$\mathrm{Spy \ plot \ of} \ \lambda E + A^{I_k}, k \geq 2 \; (N=11)$")
# savefig(p, "spy-2.pdf")

uaa(x,t) = ASp[x,1:length(u[1])]'*u[Int64(round(t/Δt))][1:end] 
x = -3:0.01:3; t = 0.01:0.01:1.01;
X = repeat(reshape(x, 1, :), length(t), 1)
T = repeat(t, 1, length(x))
Z = map(uaa, X, T)
p = contour(x,t,Z,fill=true, rev=true,
        xlabel=L"$x$",ylabel=L"$t$",
        xtickfontsize=8, ytickfontsize=8,xlabelfontsize=15,ylabelfontsize=15)
savefig(p, "W0-frac-heat-contour.pdf")