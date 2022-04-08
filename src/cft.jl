function inverse_fourier_transform(F, t; t0=-1000, dt=0.001)
    w = ifftshift(fftfreq(length(t), 1/dt) * 2 * pi)
    
    f = ifftshift(ifft(F(t)))
    f = f .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))

    return (w[2:end], f[2:end])
end


function supporter_functions(λ, μ, η, N; t0=-1000., dt=0.001, a=[-1.,1.])
    t=range(t0,-t0,step=dt)

    # FIXME: Should probably find a nicer way around this
    if μ ≈ 0 && real(λ) < 0 && λ ∈ t
        λ += 1e2*eps() * π
        @warn "μ ≈ 0, λ < 0, and λ ∈ Sample, we slightly perturb λ to avoid NaNs in the FFT computation in supporter_functions in cft.jl"
    end

    fm = k -> (λ .- im.*μ.*sign.(k) .+im.*η.*k .+ abs.(k))
    hfm = k -> (- λ.*im.*sign.(k) .- μ .+ η.*abs.(k) .- im.*k)
    
    # For certain values of λ, μ and η, fm can be equal to 0. This causes
    # NaNs in the ifft routine. To avoid this we watch if fm=0 and if it does
    # then we average the Fourier transform around 0 to get an approximate of the limit.
    # FIXME: These list comprehensions are slow, can we speed it up? 

    ## Define function F[wT0] / (λ - i*μ*sgn(k)+ iηk + abs.(k))
    # λ = μ = 0, η = -1 FwT0(0) → ∞; λ = μ = 0, η = 1 FwT0(0) → ∞ ± i∞
    # NaNs if λ = 0, but rest can be set to 0.
    tFwT0 = k -> pi * besselj.(0, abs.(k)) ./ fm(k)
    if λ ≈ 0
        FwT0 = k -> [fm(m) ≈ 0 ? (tFwT0(m-eps())+tFwT0(m+eps()))/2 : tFwT0(m)  for m in k]
    else
        FwT0 = k -> tFwT0(k)
    end
    
    ## Define function F[wT1] / (λ - i*μ*sgn(k)+ iηk + abs.(k))
    # λ = μ = 0, η = -1 FwT1(0) → c; λ = μ = 0, η = 1 FwT1(0) → c ± ic, where c finite
    # λ = μ = η = 0, FwT1(0) = ± c
    # NaNs if λ = 0, but rest can be set to 0.
    tFwT1 = k -> (-im).^(N+2) * pi * besselj.(N+2, k) ./ fm(k)
    if λ ≈ 0
        FwT1 = k -> [fm(m) ≈ 0 ? (tFwT1(m-eps())+tFwT1(m+eps()))/2 : tFwT1(m)  for m in k]
    else
        FwT1 = k -> tFwT1(k)
    end
    ## Define function F[U0] / (λ - i*μ*sgn(k)+ iηk + abs.(k))
    # λ = μ = 0, η = -1 FU0(0) → c ∓ ic; λ = μ = 0, η = 1 FU0(0) → c ± ic, where c finite
    # λ = μ = η = 0, FU0(0) = c
    # NaNs if λ = 0, but rest can be set to 0.
    
    tFU0 =  k -> ( (-im).^(N+2)*pi.*besselj.(N+2, k) ) ./ hfm(k)
    FU0 = k -> [hfm(m) ≈ 0 ? (tFU0(m-eps())+tFU0(m+eps()))/2 : tFU0(m)  for m in k]

    # FU0 = k -> [m ≈ 0 ? (tFU0(-eps())+tFU0(eps()))/2 : tFU0(m)  for m in k]
    
    ## Define function F[U-1] / (λ - i*μ*sgn(k)+ iηk + abs.(k))
    tFU_1 = k -> ( im*pi*k.*besselj.(0,abs.(k)) ./ abs.(k)) ./ fm(k)
    # FU_1 = k -> [m ≈ 0 ? ComplexF64(0.) : ComplexF64(im*pi*m.*besselj.(0,abs.(m)) ./ abs.(m)) ./ fm(m) for m in k]
    FU_1 = k -> [m ≈ 0 ? (tFU_1(-eps())+tFU_1(eps()))/2 : tFU_1(m)  for m in k]
    
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

    return (w[2:end], ywT0, yU_1, ywT1, yU0) 
end

function cifft(f, t, dt, t0, w)
    yf= ifftshift(ifft(f(t)))
    yf = yf .* dt .* exp.(-im .*w .*t0) .* length(t) ./ ((2*pi))
    return yf
end

function interpolate_supporter_functions(w, ywT0, yU_1, ywT1, yU0)
    el_no = length(yU_1)
    yU_1 = [interpolate((w,), real.(yU_1[j]), Gridded(Linear())) for j in 1:el_no]
    yU0 = [interpolate((w,), real.(yU0[j]), Gridded(Linear())) for j in 1:el_no]
    ywT0 = [interpolate((w,), real.(ywT0[j]), Gridded(Linear())) for j in 1:el_no]
    ywT1 = [interpolate((w,), real.(ywT1[j]), Gridded(Linear())) for j in 1:el_no]
    return (ywT0, yU_1, ywT1, yU0)
end

function fft_supporter_functions(λ, μ, η, N; t0=-1000., dt=0.001, a=[-1.,1.])
    # Special case analytical expressions
    if λ == μ == η ≈ 0
        ywT0 = []; ywT1 = []; yU0 = []; yU_1 = []
        for els = 1 : length(a)-1
            append!(ywT0, [x->half_laplace_wT0(affinetransform(a[els], a[els+1], x))])
            append!(ywT1, [x->half_laplace_wT1(affinetransform(a[els], a[els+1], x))])
            append!(yU0, [x->half_laplace_U0(affinetransform(a[els], a[els+1], x))])
            append!(yU_1, [x->half_laplace_U_1(affinetransform(a[els], a[els+1], x))])
            return (ywT0, yU_1, ywT1, yU0)
        end
    end 
    
    (w, ywT0, yU_1, ywT1, yU0) = supporter_functions(λ, μ, η, N, t0=t0, dt=dt, a=a)
    uS = interpolate_supporter_functions(w, ywT0, yU_1, ywT1, yU0)
    return uS
end

function coefficient_supporter_functions(A, x, uS, N; tol=1e-6)
    (ywT0, yU_1, ywT1, yU0) = uS
    el_no = length(yU0)
    yu_1 = [solvesvd(A, riemann(x, yU_1[j]); tol=tol) for j in 1:el_no]
    yu0 = [solvesvd(A, riemann(x, yU0[j]); tol=tol) for j in 1:el_no]
    ywt0 = [solvesvd(A, riemann(x, ywT0[j]); tol=tol) for j in 1:el_no]
    ywt1 = [solvesvd(A, riemann(x, ywT1[j]); tol=tol) for j in 1:el_no]
    yu_1 = [expansion_sum_space(yu_1[j], N, el_no) for j in 1:el_no]
    yu0 = [expansion_sum_space(yu0[j],  N, el_no) for j in 1:el_no]
    ywt0 = [expansion_sum_space(ywt0[j], N, el_no) for j in 1:el_no]
    ywt1 = [expansion_sum_space(ywt1[j],  N, el_no) for j in 1:el_no]
    return (ywt0, yu_1, ywt1, yu0)
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