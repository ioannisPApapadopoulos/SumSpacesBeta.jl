using SumSpaces, SpecialFunctions
using Plots

Δt = 1e-2
t0=-1000; dt=0.001
t=range(t0,-t0,step=dt)

# Fwu0 = k -> pi * besselj.(1, abs.(k)) ./ (abs.(k) .* (1 .+ Δt .* abs.(k)))

Fwu0 = (k, n) -> [m ≈ 0. ? ComplexF64(pi/2) : ComplexF64(pi * besselj.(1, abs.(m)) ./ (abs.(m) .* (1 .+ Δt .* abs.(m)).^n)) for m in k]
fv = []
for n = 1:10
    Fwu0n = k -> Fwu0(k, n)
    (x, f) = inverse_fourier_transform(Fwu0n, t)
    append!(fv, [f])
end

xx = x[-5 .< x .< 5]
p = plot()
for n = 1: 10
    t = Δt*n
    xlim = [-5,5]; ylim = [-0.1,1]
    p = plot(xx,real.(fv[n][-5 .< x .< 5]), label="time=$t (s)", legend=:topleft, xlim=xlim, ylim=ylim)
    sleep(0.1)
    display(p)
end  