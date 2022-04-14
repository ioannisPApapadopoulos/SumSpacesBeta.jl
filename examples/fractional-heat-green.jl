using Revise
using SumSpaces, Plots
import LinearAlgebra: I

"""
Solve the fractional heat equation with one element at [-1,1]. 
"""

N = 101 # Truncation degree
λ = 1e2; μ = 0; η = 0; Δt = 1/λ # Constants
a = [-50.,50]

Sp = SumSpaceP(a) # Primal sum space
Sd = SumSpaceD(a) # Dual sum space

M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me, endpoints=20) # Collocation points
Nn = min(N,7)

A = framematrix(x, Sp, Nn, M, Me) # Blocked frame matrix

# Compute support functions
uS = fft_supporter_functions(λ, μ, η, a=a) # Actual functions
# Primal sum space coefficients
cuS = coefficient_supporter_functions(A, x, uS, N) 

# Create appended sum space
ASp = AppendedSumSpace(uS, cuS, a)

# Sanity check plot
# xx = -5:0.01:5
# plot(xx, uS[1][1](xx))
# plot!(xx, Sp[xx,Block.(1:N+2)]*cuS[1][1])


x = axes(Sp, 1); H = inv.( x .- x')
Hm = (1/π).*((Sp\(H*Sp))[1:2N+3,1:2N+3])    # Hilbert: Sp -> Sp
Cm = (Sd\(Derivative(x)*Sp))[1:2N+7,1:2N+3] # Derivative: Sp -> Sd
Bm = (Sd\Sp)[1:2N+7,1:2N+3]                 # Identity: Sp -> Sd

Id = (Sd \ ASp)[1:2N+7,1:2N+7]  # Identity: ASp -> Sd


Dm =  λ.*Bm + μ.*Bm*Hm + Cm*Hm     # Helmholtz-like operator: Sp -> Sd   
Dm = hcat(zeros(size(Dm,1), 4),Dm) # Adding 4 columns to construct: ASp -> Sd
Dm[2:5,1:4] = I[1:4,1:4]
Dm = [Dm[:,5] Dm[:,1:4] Dm[:,6:end]] # permute T0 column to start

# Initial condition, u₀ = √(1-x²)
u0 = x -> 1. ./ (x.^2 .+ 1) 
# u₀ = zeros(2*N+7)
# u₀[6] = 1.
x = collocation_points(M, Me, endpoints=20)
A = framematrix(x, Sp, N, M, Me)

# u₀ = solvesvd(A, riemann(x, u0))
# u₀₀= zeros(2N+7); u₀₀[1]=u₀[1]; u₀₀[6:end]=u₀[2:end]
# u = [u₀₀]

T = chebyshevt(-50..50)
u1 = (T \ u0.(axes(T,1)))[1:N+2]
u11 = zeros(2N+7); u11[1] = u1[1]; u11[7:2:end] = u1[2:end]
u = [u11]
# Run solve loop for time-stepping
timesteps=50
for k = 1:timesteps

    # Map from Sp to Sd
    v = Id * u[k]
    # Multiply RHS with λ
    v = λ.*v

    # Solve for coefficients in ASp and append to list
    u1 = Dm \ v
    append!(u,  [u1])

    # Normalise constant so that it decays to zero
    u[k+1][1] = u[k+1][1] - ASp[1e2,1:length(u[k+1])]'*u[k+1]
end


# Plot solution
p = plot() 
xx = -15:0.01:15
xlim = [xx[1],xx[end]]; ylim = [-0.1,1]
# anim = @animate  for k = 2: timesteps+1
for k = 1:  timesteps+1
    t = Δt*(k-1)
    y = x -> (1 + t) ./ (x.^2 .+ (1+t).^2)
    tdisplay = round(t, digits=2)
    yy = ASp[xx,1:length(u[k])]*u[k]
    p = plot(xx,yy, title="time=$tdisplay (s)", label="Sum space - 1 element", legend=:topleft, xlim=xlim, ylim=ylim)
    p = plot!(xx, y(xx), label="True solution", legend=:topleft, xlim=xlim, ylim=ylim)
    sleep(0.1)
    display(p)
end  
# gif(anim, "anim_fps10.gif", fps = 10)





