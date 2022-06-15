using Revise
using SumSpacesBeta
using Plots

"""
In this script we solve the Helmholtz equation 

-u''(x) - k²u(x) = √(1-x^2), u(x) → 0 as |x| → ∞

via two solves of the fractional Helmholtz equation. This makes sense as

(√|Δ|-k)(√|Δ|+k) = (d/dx) H (d/dx) H - k² = -(d/dx)(d/dx) - k²= -Δ-k². 
"""

N = 5 # Truncation degree
k = 5 # Constants

Sp = SumSpaceP() # Primal sum space
Sd = SumSpaceD() # Dual sum space

M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me) # Collocation points
Nn = min(N,7)
A = framematrix(x, Sp, Nn, M, Me) # Blocked frame matrix

# Compute support functions and expansions
uSp = fft_supporter_functions(k, 0) # Actual functions
cuSp = coefficient_supporter_functions(A, x, uSp, Nn, N) 

uSm = fft_supporter_functions(-k, 0) # Actual functions
cuSm = coefficient_supporter_functions(A, x, uSm, Nn, N) 

# Create appended sum space
ASpp = AppendedSumSpace(uSp, cuSp)
ASpm = AppendedSumSpace(uSm, cuSm)

x = axes(Sp, 1); H = inv.( x .- x')
Hm = (1/π).*((Sp\(H*Sp))[1:2N+3,1:2N+3])    # Hilbert: Sp -> Sp
Cm = (Sd\(Derivative(x)*Sp))[1:2N+7,1:2N+3] # Derivative: Sp -> Sd
Bm = (Sd\Sp)[1:2N+7,1:2N+3]                 # Identity: Sp -> Sd

Id = (Sd \ ASpp)[1:2N+7,1:2N+7]  # Identity: ASpp -> Sd

Dm = []
append!(Dm, [k.*Bm + Cm*Hm])    
append!(Dm, [-k.*Bm + Cm*Hm])  
Dm = [hcat(zeros(size(Dm[j],1), 4),Dm[j]) for j = 1:2] # Adding 4 columns to construct: ASp -> Sd
Dm[1][2:5,1:4] = I[1:4,1:4]; Dm[2][2:5,1:4] = I[1:4,1:4]
Dm = [[Dm[j][:,5] Dm[j][:,1:4] Dm[j][:,6:end]] for j = 1:2] # permute T0 column to start

# Initial condition, u₀ = √(1-x²)
u₀ = zeros(2*N+7)
u₀[6] = 1.


u1 = Dm[1] \ u₀       # First solve
u2 = Dm[2] \ (Id*u1)  # Map back to dual sum space and second solve

# Plot solution
p = plot()
xx = -5:0.001:5
yy = ASpp[xx,1:length(u1)]*u1

xlim = [xx[1],xx[end]]; ylim = [-0.2,0.4]
p = plot(xx,yy, label="solution", legend=:topleft, xlim=xlim, ylim=ylim)
display(p)
