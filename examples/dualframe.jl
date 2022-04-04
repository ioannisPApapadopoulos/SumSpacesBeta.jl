using SumSpaces
using ClassicalOrthogonalPolynomials, Plots

N = 21
Tp = Float64
λ = 1e2
μ = 0
Δt = 1/λ

ewU = ExtendedWeightedChebyshevU{Tp}()
eT = ExtendedChebyshevT{Tp}()
ewT = ExtendedWeightedChebyshevT{Tp}()
eU = ExtendedChebyshevU{Tp}()


M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me, innergap=1e-10) # Collocation points
Nn = min(N,5)
A = dualframematrix(x, eU, ewT, Nn, M, Me)

(w, yU_1, yU0, ywT0, ywT1) = supporter_functions(λ, μ)
(yU_1, yU0, ywT0, ywT1) = interpolate_supporter_functions(w, yU_1, yU0, ywT0, ywT1)
(yu_1, yu0, ywt0, ywt1) = columns_supporter_functions(A, x, yU_1, yU0, ywT0, ywT1, N+4, N+3, Nn+4, Nn+3)

xx = -5:0.01:5
plot(xx, yU0(xx))
plot!(xx, dual_sum_space(eU, ewT, yu0, xx, N))


Id = idmap_sum2dual(N, Tp)
D = fractionalhelmholtzmap(λ, μ, N, Tp)
# D = D[1:end-1,:]
# D = vcat(zeros(1,size(D,2)),D)
# endpoint = 10
# D[1,1:N+2] = eT[endpoint,1:N+2]
# D[1,N+3:2*N+3] = ewU[endpoint,1:N+1]
# D[1,2*N+4] = yU0(endpoint)
# D[1,2*N+5] = yU_1(endpoint)
# D[1,2*N+6] = ywT0(endpoint)
# D[1,2*N+7] = ywT1(endpoint)

u₀ = zeros(2*N+7)
u₀[N+3] = 1.
u = [u₀]

timesteps=50
for k = 1:timesteps
    v = Id * u[k][1:end-4]
    v = (v 
        + u[k][end-3].* yu0 
        + u[k][end-2].* yu_1 
        + u[k][end-1].* ywt0 
        + u[k][end].* ywt1)
    # v = v[1:end-1]
    # v = vcat([0.],v)
    v = λ .* v
    append!(u,  [D \ v])
    u[k+1][1] = u[k+1][1] - appended_sum_space(eT, ewU, yU0, yU_1, ywT0, ywT1, u[k+1], [1e2], N)[1]
end

p = plot()
# xx1 = w[-5 .< w .< 5]
for k = 1: timesteps+1
    t = Δt*(k-1)
    xx = -5:0.001:5
    yy = appended_sum_space(eT, ewU, yU0, yU_1, ywT0, ywT1, u[k], xx, N)
    xlim = [xx[1],xx[end]]; ylim = [-0.1,1]
    p = plot(xx,yy, label="time=$t (s)", legend=:topleft, xlim=xlim, ylim=ylim)
    # p = plot!(xx1,real.(fv[k-1][-5 .< w .< 5]), label="time=$t (s)", legend=:topleft, xlim=xlim, ylim=ylim)
    sleep(0.1)
    display(p)
end  