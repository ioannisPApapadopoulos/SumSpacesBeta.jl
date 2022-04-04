using Revise
using SumSpaces, Plots
import LinearAlgebra: I

N = 5 # Truncation degree
λ = 1e2; μ = 0; Δt = 1/λ # Constants

Sp = SumSpaceP() # Primal sum space
Sd = SumSpaceD() # Dual sum space

M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me) # Collocation points
Nn = min(N,7)

A = framematrix([x], Sp, Nn, M, Me) # Blocked frame matrix

# Helper functions for computing support functions
(w, yU_1, yU0, ywT0, ywT1) = supporter_functions(λ, μ)
(yU_1, yU0, ywT0, ywT1) = interpolate_supporter_functions(w, yU_1, yU0, ywT0, ywT1)
(yu_1, yu0, ywt0, ywt1) = columns_supporter_functions(A, x, yU_1, yU0, ywT0, ywT1, Nn, N)

# Create appended sum space
ASp = AppendedSumSpace((ywT0, yU_1, ywT1, yU0), (ywt0, yu_1, ywt1, yu0))

# Quick test
xx = -5:0.01:5
plot(xx, ywT0[1](xx))
plot!(xx, Sp[xx,Block.(1:N+2)]*ywt0[1])


x = axes(Sp, 1); H = inv.( x .- x')
Hm = (1/π).*((Sp\(H*Sp))[1:2N+3,1:2N+3])    # Hilbert: Sp -> Sp
Cm = (Sd\(Derivative(x)*Sp))[1:2N+7,1:2N+3] # Derivative: Sp -> Sd
Bm = (Sd\Sp)[1:2N+7,1:2N+3]                 # Identity: Sp -> Sd

Id = (Sd \ ASp)[1:2N+7,1:2N+7]  # Identity: ASp -> Sd


Dm =  λ.*Bm + μ.*Bm*Hm + Cm*Hm     # Helmholtz-like operator: Sp -> Sd   
Dm = hcat(zeros(size(Dm,1), 4),Dm) # Adding 4 columns to construct: ASp -> Sd
Dm[2:5,1:4] = I[1:4,1:4]

# Initial condition, u₀ = √(1-x²)
u₀ = zeros(2*N+7)
u₀[6] = 1.
u = [u₀]

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
    u[k+1][5] = u[k+1][5] - ASp[1e2,1:length(u[k+1])]'*u[k+1]
end


# Plot solution
p = plot()
for k = 1: timesteps+1
    t = Δt*(k-1)
    xx = -5:0.001:5
    yy = ASp[xx,1:length(u[k])]*u[k]

    xlim = [xx[1],xx[end]]; ylim = [-0.1,1]
    p = plot(xx,yy, label="time=$t (s)", legend=:topleft, xlim=xlim, ylim=ylim)
    sleep(0.1)
    display(p)
end  