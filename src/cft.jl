using FFTW, SpecialFunctions, Interpolations

function inverse_fourier_transform(F, t; t0=-1000, dt=0.001)
    w = ifftshift(fftfreq(length(t), 1/dt) * 2 * pi)
    
    f = ifftshift(ifft(F(t)))
    f = f .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    return (w[2:end], f[2:end])
end


function supporter_functions(λ, μ; t0=-1000., dt=0.001, a=[-1.,1.])
    t=range(t0,-t0,step=dt)
    
    ## Define function F[wT0] / (λ - i*μ*sgn(k)+  abs.(k))
    FwT0 = k -> pi * besselj.(0, abs.(k)) ./ (λ .- im.*μ.*sign.(k) .+ abs.(k))
    ## Define function F[wT1] / (λ - i*μ*sgn(k)+  abs.(k))
    FwT1 = k -> -im * pi * besselj.(1, k) ./ (λ .- im.*μ.*sign.(k) .+ abs.(k))
    ## Define function F[wT1] / (-i*λ*sign.(k) - μ - i*k)
    FU_2 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64(- im * pi * besselj.(1, m) ./ (-im * λ * sign.(m) .- μ .- im * m )) for m in k]
    ## Define function F[wT0] / (-i*λ*sign.(k) - μ - i*k)
    FU_1 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64(pi * besselj.(0, abs.(m)) ./ (-im * λ * sign.(m) .- μ .- im * m )) for m in k]

    # Experimental for decaying support
    # FU_2 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64((- im * pi * besselj.(1, m) + 2 *im * sin.(m) ./ abs.(m)) ./ (-im * sign.(m) .- im * Δt .* m )) for m in k]

    # scaling if element is not [-1,1]
    sFwT0 = k -> 0.5*(a[2]-a[1]).*exp.(-im.*k.*(a[2]^2-a[1]^2)/4).*FwT0((a[2]-a[1])/2 .*k)
    sFwT1 = k -> 0.5*(a[2]-a[1]).*exp.(-im.*k.*(a[2]^2-a[1]^2)/4).*FwT1((a[2]-a[1])/2 .*k)
    sFU_2 = k -> 0.5*(a[2]-a[1]).*exp.(-im.*k.*(a[2]^2-a[1]^2)/4).*FU_2((a[2]-a[1])/2 .*k)
    sFU_1 = k -> 0.5*(a[2]-a[1]).*exp.(-im.*k.*(a[2]^2-a[1]^2)/4).*FU_1((a[2]-a[1])/2 .*k)



    w = ifftshift(fftfreq(length(t), 1/dt) * 2 * pi)
    
    ywT0 = ifftshift(ifft(sFwT0(t)))
    ywT0 = ywT0 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    ywT1 = ifftshift(ifft(sFwT1(t)))
    ywT1 = ywT1 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    yU_2 = ifftshift(ifft(sFU_2(t)))
    yU_2 = yU_2 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    yU_1 = ifftshift(ifft(sFU_1(t)))
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

function fractional_heat_fourier_solve(F, t, timesteps)
    
    Fn = k -> F(k, 1)
    (x, f) = inverse_fourier_transform(Fn, t)
    
    fv = [f]
    for n = 2:timesteps
        Fn = k -> F(k, n)
        (_, f) = inverse_fourier_transform(Fn, t)
        append!(fv, [f])
    end
    return (x, fv)
end