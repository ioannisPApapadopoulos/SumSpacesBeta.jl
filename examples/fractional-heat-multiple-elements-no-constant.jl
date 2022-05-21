using Revise
using SumSpaces
using ClassicalOrthogonalPolynomials, Plots

N = 5
Tp = Float64
λ = 1e2
μ = 0
Δt = 1/λ

ewU = ExtendedWeightedChebyshevU{Tp}()
eT = ExtendedChebyshevT{Tp}()
a = [-3,-1.,1.,3]
el_no = length(a)-1


M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me, endpoints=10) # Collocation points
x1 = affinetransform(a[1],a[2],x)
x2 = affinetransform(a[2],a[3],x)
x3 = affinetransform(a[3],a[4],x)

Nn = min(N,11)
A = framematrix([x1,x2,x3], eT, ewU, Nn, M, Me)
A = A[:,2:end]

(w, yU_1, yU0, ywT0, ywT1) = supporter_functions(λ, μ, a=a)
(yU_1, yU0, ywT0, ywT1) = interpolate_supporter_functions(w, yU_1, yU0, ywT0, ywT1)
(yu_1, yu0, ywt0, ywt1) = columns_supporter_functions(A, x, yU_1, yU0, ywT0, ywT1, Nn, N, constant=false)

xx = -5:0.01:5
plot(xx, sum_space(eT, ewU, yu0[1], xx, N, a=a))
plot!(xx, yU0[1](xx))

Id = idmap_append2dual(N, yu0, yu_1, ywt0, ywt1, el_no=el_no, Tp=Tp)
D = fractionalhelmholtzmap(λ, μ, N, a=a,Tp=Tp)
D = D[2:end, 2:end]
D[2,end-11] = 1.

# Ds = split_block_helmholtz_matrix(D, el_no)

u₀ = zeros(1+el_no*(2*N+6))
u₀[3N+5] = 1. # Coefficient corresponding to sqrt(1-x^2) on |x| ≤ 1. 
u = [u₀]

timesteps=100
for k = 1:timesteps
    v = Id * u[k]
    v = λ.*v

    u1 = D \ v[2:end]
    
    # append!(u,  [u1])
    append!(u,  [vcat([0.],u1)])
    # u[k+1][1] = u[k+1][1] - appended_sum_space(eT, ewU, yU0, yU_1, ywT0, ywT1, u[k+1], [1e3], N)[1]
end

# for k = 1:timesteps
#     v = Id * u[k]
#     v = λ.*v
#     vs = split_block_helmholtz_vector(v, el_no)

#     us = []
#     uss = []
#     for e in 1:el_no
#         ut = Ds[e] \ vs[e]
#         append!(us, ut[1:end-4])
#         append!(uss, ut[end-3:end])
#     end
    
#     append!(u,  [append!(us, uss)])
#     # append!(u,  [uss])
#     u[k+1][1] = u[k+1][1] - appended_sum_space(eT, ewU, yU0, yU_1, ywT0, ywT1, u[k+1], [1e3], N)[1]
# end

b = 20.
# xx1 = w[-b .< w .< b]

p = plot()
# anim = @animate  for k = 2: timesteps+1
for k = 1: timesteps+1
    t = round(Δt*(k-1), digits=2)
    xx = -10:0.001:10
    yy = appended_sum_space(eT, ewU, yU0, yU_1, ywT0, ywT1, u[k], xx, N,a=a)
    xlim = [xx[1],xx[end]]; ylim = [-0.1,1]
    p = plot(xx,yy, title="time=$t (s)", label="Sum space - 3 elements", legend=:topleft, xlim=xlim, ylim=ylim)
    # p = plot!(xx1,real.(fv[k-1][-b .< w .< b]), label="Fourier transform method", legend=:topleft, xlim=xlim, ylim=ylim)

    sleep(0.02)
    display(p)
end  
# gif(anim, "anim_fps10.gif", fps = 10)