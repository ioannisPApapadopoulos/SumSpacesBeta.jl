using Revise
using SumSpaces
using ClassicalOrthogonalPolynomials, Plots

"""
In this script we solve the Poisson equation -u''(x) = (4x^2-1)√(1-x^2), u(x) → 0 as |x| → ∞
via two solves of the half-Laplacian. This makes sense as

√|Δ| √|Δ| = (d/dx) H (d/dx) H = -(d/dx)(d/dx) = -Δ. 
"""

N = 5
λ = 0
μ = 0

ewU = ExtendedWeightedChebyshevU()
eT = ExtendedChebyshevT()
ewT = ExtendedWeightedChebyshevT()
eU = ExtendedChebyshevU()

M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me, innergap=1e-10) # Collocation points
Nn = N
A = framematrix(x, eT, ewU, Nn, M, Me)

# (w, yU0, yU_1, ywT0, ywT1) = supporter_functions(λ, μ)
# (yU0, yU_1, ywT0, ywT1) = interpolate_supporter_functions(w, yU0, yU_1, ywT0, ywT1)

yU_1 = [sqrt_laplace_U_1]
yU0 = [x -> zeros(length(x))]
ywT0 = [sqrt_laplace_wT0]
ywT1 = [x -> zeros(length(x))]


(yu0, yu_1, ywt0, ywt1) = columns_supporter_functions(A, x, yU0, yU_1, ywT0, ywT1, Nn, N)

xx = -5:0.01:5
plot(xx,sum_space(eT, ewU, ywt0[1], xx, N))
plot!(xx, ywT0[1](xx))

Id = idmap_append2dual(N, yu0, yu_1, ywt0, ywt1)
D = fractionallaplacemap(N)

f = x -> ewU[x,3]

Ad = dualframematrix(x, eU, ewT, Nn+3, Nn+2, M, Me)[:,2:end]
fc = solvesvd(Ad, riemann(x, f); tol=1e-6)
xx = -5:0.001:5
plot(xx,f(xx))
plot!(xx, dual_sum_space2(eU, ewT, fc, xx, Nn))

u = D \ fc
u = vcat([0.], u)
u = vcat(vcat(vcat(u[1:end-1], [0.]), u[end]), [0.])

plot(xx, appended_sum_space(eT, ewU, yU0, yU_1, ywT0, ywT1, u, xx, N; a=[-1.,1.]))

Id = idmap_append2dual(N, yu0, yu_1, ywt0, ywt1; el_no=1)
du = Id * u
plot!(xx, dual_sum_space(eU, ewT, du, xx, Nn))
du = [du[1:N+3]' du[N+5:end]']'
du = du[2:end-1]


u2 = D \ du
u2 = vcat([0.], u2)
u2 = vcat(vcat(vcat(u2[1:end-1], [0.]), u2[end]), [0.])

y = xx -> [abs(x) < 1 ? -(1 - x^2)^(5/2)/5 : 0. for x in xx]
plot(xx,y(xx))
plot!(xx, appended_sum_space(eT, ewU, yU0, yU_1, ywT0, ywT1, u2, xx, N; a=[-1.,1.]))