using Revise
using ClassicalOrthogonalPolynomials, SumSpacesBeta, Plots
import LinearAlgebra: I

"""
Solve the fractional heat equation with one element at [-1,1]. 
"""

N = 61 # Truncation degree
λ = 0; μ = 0; η = 0 # Constants
a = [-1.,1]

Sp = SumSpaceP(a) # Primal sum space
Sd = SumSpaceD(a) # Dual sum space

M = max(N^2,8001)  # Number of collocation points in [-1,1]
Me = M # Number of collocation points in [-2,-1) and (1,2].
Nn = N

x = axes(Sp, 1); H = inv.( x .- x')
Hm = (1/π).*((Sp\(H*Sp))[1:2N+3,1:2N+3])    # Hilbert: Sp -> Sp
Cm = (Sd\(Derivative(x)*Sp))[1:2N+7,1:2N+3] # Derivative: Sp -> Sd

Id = (Sd \ ASp)[1:2N+7,1:2N+7]  # Identity: ASp -> Sd
# Id[:,2:5] .= 0

Dm = Cm*Hm     # Helmholtz-like operator: Sp -> Sd   
Dm = Dm[4:end-2,2:end]
# Initial condition, u₀ = √(1-x²)
u0 = x -> Sd[x,2]
x = collocation_points(M, Me,endpoints=[-10.,10.],innergap=1e-5)
A = dualframematrix(x, Sd, Nn-1) # Blocked frame matrix
A = A[:,4:end]
# A = A[:,1:2:end]
u₀ = solvesvd(A, riemann(x, u0),tol=1e-10)
# v = zeros(size(Dm,1)); v[1:2:end] = u₀
v = u₀
# x = collocation_points(M, Me,endpoints=[-5.,5.],innergap=1e-5)
# A = framematrix(x, Sp, Nn-3) # Blocked frame matrix
# u₀ = solvesvd(A, riemann(x, u0))


u  = Dm \ v
u = vcat(zeros(1), u)
# u[k+1][1] = u[k+1][1] - ASp[50.,1:length(u[k+1])]'*u[k+1]



# Plot solution
p = plot() 
xx = -10:0.01:10
xlim = [xx[1],xx[end]]; ylim = [-1,1]
# anim = @animate  for k = 2: timesteps+1

yy = Sp[xx,1:length(u)]*u
yyd = Sd[xx,4:length(v)+3]*v
# yyp = Sp[xx,1:length(u₀)]*u₀
p = plot(xx,yy, title="u_V0", label="Sum space - 1 element", legend=:topleft)#, xlim=xlim, ylim=ylim)
p = plot(xx,yyd,  label="approx-Sd", legend=:topleft)
# p = plot(xx,yyp,  label="approx-Sp", legend=:topleft)
p = plot!(xx, Sd[xx,2])





