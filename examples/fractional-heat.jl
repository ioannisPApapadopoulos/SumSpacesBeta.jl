using SumSpaces
using ClassicalOrthogonalPolynomials, Plots

N = 11
Tp = Float64
λ = 1e2
μ = 0
Δt = 1/λ

ewU = ExtendedWeightedChebyshevU{Tp}()
eT = ExtendedChebyshevT{Tp}()

M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me) # Collocation points
Nn = min(N,5)
A = framematrix(x, eT, ewU, Nn, M, Me)

(w, yU_2, yU_1, ywT0, ywT1) = supporter_functions(λ, μ)
(yU_2, yU_1, ywT0, ywT1) = interpolate_supporter_functions(w, yU_2, yU_1, ywT0, ywT1)
(yu_2, yu_1, ywt0, ywt1) = columns_supporter_functions(A, x, yU_2, yU_1, ywT0, ywT1, N+2, N+1, Nn+2, Nn+1)

Id = idmap_append2dual(N, yu_2, yu_1, ywt0, ywt1, Tp)
D = fractionalhelmholtzmap(λ, μ, N, Tp)

u₀ = zeros(2*N+7)
u₀[N+3] = 1.
u = [u₀]

timesteps=50
# Ap = framematrix(x, eT, ewU, N, M, Me)
for k = 1:timesteps
    v = Id * u[k]
    v = λ.*v

    u1 = D \ v
    # f = x -> appended_sum_space(eT, ewU, yU_2, yU_1, ywT0, ywT1, u1, x, N)
    # u1 = solvesvd(Ap, riemann(x, f); tol=1e-2)
    # u1 = vcat(u1, zeros(4))
    
    append!(u,  [u1])
    u[k+1][1] = u[k+1][1] - appended_sum_space(eT, ewU, yU_2, yU_1, ywT0, ywT1, u[k+1], [1e2], N)[1]
    # f = x -> appended_sum_space(eT, ewU, yU_2, yU_1, ywT0, ywT1, u[k+1], x, N)
end



# xx1 = w[-5 .< w .< 5]

p = plot()
for k = 2: timesteps+1
    t = Δt*(k-1)
    xx = -5:0.001:5
    yy = appended_sum_space(eT, ewU, yU_2, yU_1, ywT0, ywT1, u[k], xx, N)
    xlim = [xx[1],xx[end]]; ylim = [-0.1,1]
    p = plot(xx,yy, label="time=$t (s)", legend=:topleft, xlim=xlim, ylim=ylim)
    # p = plot!(xx1,real.(fv[k-1][-5 .< w .< 5]), label="time=$t (s)", legend=:topleft, xlim=xlim, ylim=ylim)

    sleep(0.1)
    display(p)
end  