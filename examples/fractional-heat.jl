using SumSpaces
using ClassicalOrthogonalPolynomials, Plots

N = 21
Tp = Float64
Δt = 1e-2

ewU = ExtendedWeightedChebyshevU{Tp}()
eT = ExtendedChebyshevT{Tp}()

M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me) # Collocation points
Nn = min(N,9)
A = framematrix(x, eT, ewU, Nn, M, Me)

(w, yU_2, yU_1, ywT0, ywT1) = supporter_functions(Δt)
(yU_2, yU_1, ywT0, ywT1) = interpolate_supporter_functions(w, yU_2, yU_1, ywT0, ywT1)
(yu_2, yu_1, ywt0, ywt1) = columns_supporter_functions(A, x, yU_2, yU_1, ywT0, ywT1, N, Nn)

Id = idmap_append2dual(N, yu_2, yu_1, ywt0, ywt1, Tp)
D = fractionalhelmholtzmap(Δt, N, Tp)

u₀ = zeros(2*N+7)
u₀[N+3] = 1.
u = [u₀]

timesteps=50
for k = 1:timesteps
    v = Id * u[k]
    append!(u,  [D \ v])
end

p = plot()
for k = 1 : timesteps+1
    t = Δt*(k-1)
    xx = -5:0.001:5
    yy = sum_space(eT, ewU, yU_2, yU_1, ywT0, ywT1, u[k], xx, N)
    xlim = [xx[1],xx[end]]; ylim = [-0.1,1]
    p = plot(xx,yy, label="time=$t (s)", legend=:topleft, xlim=xlim, ylim=ylim)
    sleep(0.1)
    display(p)
end  