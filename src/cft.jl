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
    ## Define function F[U0] / (λ - i*μ*sgn(k)+  abs.(k))
    FU0 = k -> pi * besselj.(1, abs.(k)) ./ (λ .- im.*μ.*sign.(k) .+ abs.(k))
    ## Define function F[U-1] / (λ - i*μ*sgn(k)+  abs.(k))
    FU_1 = k -> [m ≈ 0 ? ComplexF64(0.) : ComplexF64(im*pi*m.*besselj.(0,abs.(m)) ./ abs.(m)) ./ (λ .- im.*μ.*sign.(m) .+ abs.(m)) for m in k]
    
    # Old code to find support solutions via a Hilbert transform first. 
    ## Define function F[wT1] / (-i*λ*sign.(k) - μ - i*k)
    # FU_2 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64(- im * pi * besselj.(1, m) ./ (-im * λ * sign.(m) .- μ .- im * m )) for m in k]
    ## Define function F[wT0] / (-i*λ*sign.(k) - μ - i*k)
    # FU_1 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64(pi * besselj.(0, abs.(m)) ./ (-im * λ * sign.(m) .- μ .- im * m )) for m in k]
    # Experimental for decaying support
    # FU_2 = k -> [m ≈ 0. ? ComplexF64(0.) : ComplexF64((- im * pi * besselj.(1, m) + 2 *im * sin.(m) ./ abs.(m)) ./ (-im * sign.(m) .- im * Δt .* m )) for m in k]

    # Scaling if element is not [-1,1]

    ## FIXME: Should be able to just do one FFT for each differently sized element. 
    ## Experiments shows differences between translated elements that are the same size...
    
    c = 2. ./ (a[2:end] - a[1:end-1]); d = -c .* (a[1:end-1] + a[2:end]) ./ 2
    el_no = length(c)

    # sFwT0 = k -> (1/c).*FwT0((1/c) .*k)
    # sFwT1 = k -> (1/c).*FwT1((1/c) .*k)
    # sFU_2 = k -> (1/c).*FU_2((1/c) .*k)
    # sFU_1 = k -> (1/c).*FU_1((1/c) .*k)

    w = ifftshift((fftfreq(length(t), 1/dt)) * 2 * pi) 
    ywT0 = []; ywT1 = []; yU0 = []; yU_1 = []
    
    for els = 1:el_no
        ci = c[els]; di = d[els]
        sFwT0 = k -> (1/ci).*exp.(im.*k.*(di/ci)).*FwT0((1/ci) .*k)
        sFwT1 = k -> (1/ci).*exp.(im.*k.*(di/ci)).*FwT1((1/ci) .*k)
        sFU0 = k -> (1/ci).*exp.(im.*k.*(di/ci)).*FU0((1/ci) .*k)
        sFU_1 = k -> (1/ci).*exp.(im.*k.*(di/ci)).*FU_1((1/ci) .*k)
        append!(ywT0, [cifft(sFwT0, t, dt, t0, w)[2:end]])
        append!(ywT1, [cifft(sFwT1, t, dt, t0, w)[2:end]])
        append!(yU0, [cifft(sFU0, t, dt, t0, w)[2:end]])
        append!(yU_1, [cifft(sFU_1, t, dt, t0, w)[2:end]])
    end
    
    # w = ifftshift((fftfreq(length(t), 1/dt)) * 2 * pi) 
    
    # ywT0 = ifftshift(ifft(sFwT0(t)))
    # ywT0 = ywT0 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    # ywT1 = ifftshift(ifft(sFwT1(t)))
    # ywT1 = ywT1 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    # yU_2 = ifftshift(ifft(sFU_2(t)))
    # yU_2 = yU_2 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    # yU_1 = ifftshift(ifft(sFU_1(t)))
    # yU_1 = yU_1 .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    return (w[2:end], yU0, yU_1, ywT0, ywT1)
end

function cifft(f, t, dt, t0, w)
    yf= ifftshift(ifft(f(t)))
    yf = yf .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))
    return yf
end



function interpolate_supporter_functions(w, yU0, yU_1, ywT0, ywT1)
    el_no = length(yU0)
    yyU0 = [interpolate((w,), real.(yU0[j]), Gridded(Linear())) for j in 1:el_no]
    yyU_1 = [interpolate((w,), real.(yU_1[j]), Gridded(Linear())) for j in 1:el_no]
    yywT0 = [interpolate((w,), real.(ywT0[j]), Gridded(Linear())) for j in 1:el_no]
    yywT1 = [interpolate((w,), real.(ywT1[j]), Gridded(Linear())) for j in 1:el_no]
    return (yyU0, yyU_1, yywT0, yywT1)
end

function columns_supporter_functions(A, x, yU0, yU_1, ywT0, ywT1, Nn, N; tol=1e-6, constant=true)
    el_no = length(yU0)
    yu_1 = [solvesvd(A, riemann(x, yU_1[j]); tol=tol) for j in 1:el_no]
    yu0 = [solvesvd(A, riemann(x, yU0[j]); tol=tol) for j in 1:el_no]
    ywt0 = [solvesvd(A, riemann(x, ywT0[j]); tol=tol) for j in 1:el_no]
    ywt1 = [solvesvd(A, riemann(x, ywT1[j]); tol=tol) for j in 1:el_no]
    yyu_1 = [expansion_sum_space(yu_1[j], Nn, N, el_no, constant) for j in 1:el_no]
    yyu0 = [expansion_sum_space(yu0[j], Nn, N, el_no, constant) for j in 1:el_no]
    yywt0 = [expansion_sum_space(ywt0[j], Nn, N, el_no, constant) for j in 1:el_no]
    yywt1 = [expansion_sum_space(ywt1[j], Nn, N, el_no, constant) for j in 1:el_no]
    return (yyu0, yyu_1, yywt0, yywt1)
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