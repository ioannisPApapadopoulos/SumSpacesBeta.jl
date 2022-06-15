using SumSpacesBeta, SpecialFunctions
using Plots

Δt = 1e-2
W=1000; δ=0.001
ω=range(-W, W,step=δ); ω = ω[1:end-1]
timesteps = 100

Fwu0 = (k, n) -> [m ≈ 0. ? ComplexF64(pi/2) : ComplexF64(pi * besselj.(1, abs.(m)) ./ (abs.(m) .* (1 .+ Δt .* abs.(m)).^n)) for m in k]
# sFwu0 = (k, n) -> [m ≈ 0. ? ComplexF64(5*pi/2) : ComplexF64(pi * besselj.(1, 5*abs.(m)) ./ (abs.(m) .* (1 .+ Δt .* abs.(m)).^n)) for m in k]

@time (fx, fv) = fractional_heat_fourier_solve(Fwu0, ω, timesteps)


# xx = x[-20 .< x .< 20]
# p = plot()
# for n = 1: timesteps
#     t = Δt*n
#     xlim = [-20,20]; ylim = [-0.1,1]
#     p = plot(xx,real.(fv[n][-20 .< x .< 20]), label="time=$t (s)", legend=:topleft, xlim=xlim, ylim=ylim)
#     sleep(0.1)
#     display(p)
# end    