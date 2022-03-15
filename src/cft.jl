using FFTW, SpecialFunctions, Interpolations

function inverse_fourier_transform(F, t; t0=-1000, dt=0.001)
    w = ifftshift(fftfreq(length(t), 1/dt) * 2 * pi)
    
    f = ifftshift(ifft(F(t)))
    f = f .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    return (w[2:end], f[2:end])
end


function supporter_functions(Δt; t0=-1000., dt=0.001)
    t=range(t0,-t0,step=dt)
    
    ## Define function F[wT0] / (1 .+ Δt .* abs.(k))
    FwT0 = k -> pi * besselj.(0, abs.(k)) ./ (1 .+ Δt .* abs.(k))
    ## Define function F[wT1] / (1 .+ Δt .* abs.(k))
    FwT1 = k -> -im * pi * besselj.(1, k) ./ (1 .+ Δt .* abs.(k))
    ## Define function F[wT1] / (-i*sign.(k) - i*Δt .* k)
    FU_2 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64(- im * pi * besselj.(1, m) ./ (-im * sign.(m) .- im * Δt .* m )) for m in k]
    ## Define function F[wT0] / (-i*sign.(k) - i*Δt .* k)
    FU_1 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64(pi * besselj.(0, abs.(m)) ./ (-im *sign.(m) .- im * Δt .* m )) for m in k]

    # Experimental for decaying support
    # FU_2 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64((- im * pi * besselj.(1, m) + 2 *im * sin.(m) ./ abs.(m)) ./ (-im * sign.(m) .- im * Δt .* m )) for m in k]


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

function columns_supporter_functions(A, x, yU_2, yU_1, ywT0, ywT1, N1, N2, Nn1, Nn2; tol=1e-6)
    yu_1 = solvesvd(A, riemann(x, yU_1); tol=tol)
    yu_2 = solvesvd(A, riemann(x, yU_2); tol=tol)
    ywt0 = solvesvd(A, riemann(x, ywT0); tol=tol)
    ywt1 = solvesvd(A, riemann(x, ywT1); tol=tol)
    yu_1 = expansion(N1, N2, Nn1, Nn2, yu_1)
    yu_2 = expansion(N1, N2, Nn1, Nn2, yu_2)
    ywt0 = expansion(N1, N2, Nn1, Nn2, ywt0)
    ywt1 = expansion(N1, N2, Nn1, Nn2, ywt1)
    return (yu_2, yu_1, ywt0, ywt1)
end