using FFTW, SpecialFunctions, Interpolations


function supporter_functions(Δt; t0=-1000., dt=0.001)
    t=range(t0,-t0,step=dt)
    
    ## Define function F[wT0] / (1 .+ Δt .* abs.(k))
    FwT0 = k -> pi * besselj.(0, abs.(k)) ./ (1 .+ Δt .* abs.(k))
    ## Define function F[wT1] / (1 .+ Δt .* abs.(k))
    FwT1 = k -> -im * pi * besselj.(1, k) ./ (1 .+ Δt .* abs.(k))
    ## Define function F[wT1] / (sign.(k) .+ Δt .* k )
    FU_2 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64(- im * pi * besselj.(1, m) ./ (-im * sign.(m) .- im * Δt .* m )) for m in k]
    ## Define function F[wT0] / (sign.(k) .+ Δt .* k )
    FU_1 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64(pi * besselj.(0, abs.(m)) ./ (-im *sign.(m) .- im * Δt .* m )) for m in k]

    w = ifftshift(fftfreq(length(t), 1/dt) * 2 * pi)
    
    ywT0 = ifftshift(ifft(FwT0(t)))
    ywT0 = ywT0 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    ywT1 = ifftshift(ifft(FwT1(t)))
    ywT1 = ywT1 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    yU_2 = ifftshift(ifft(FU_2(t)))
    yU_2 = yU_2 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    yU_1 = ifftshift(ifft(FU_1(t)))
    yU_1 = yU_1 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    return (w[2:end], yU_2[2:end], yU_1[2:end], ywT0[2:end], ywT1[2:end])
end

function interpolate_supporter_functions(w, yU_2, yU_1, ywT0, ywT1)
    yU_2 = interpolate((w,), real.(yU_2), Gridded(Linear()))
    yU_1 = interpolate((w,), real.(yU_1), Gridded(Linear()))
    ywT0 = interpolate((w,), real.(ywT0), Gridded(Linear()))
    ywT1 = interpolate((w,), real.(ywT1), Gridded(Linear()))
    return (yU_2, yU_1, ywT0, ywT1)
end