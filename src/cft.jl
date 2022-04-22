function supporter_functions(λ, μ, η; W=1000., δ=0.001, a=[-1.,1.])
    ω=range(-W, W, step=δ)
    ω = ω[1:end-1]

    # FIXME: Should probably find a nicer way around this
    if μ ≈ 0 && real(λ) < 0 && λ ∈ ω
        λ += 1e2*eps() * π
        @warn "μ ≈ 0, λ < 0, and λ ∈ Sample, we slightly perturb λ to avoid NaNs in the FFT computation in supporter_functions in cft.jl"
    end

    fmultiplier = k -> (λ .- im.*μ.*sign.(k) .+im.*η.*k .+ abs.(k))
    
    # For certain values of λ, μ and η, fmultiplier can be equal to 0. This causes
    # NaNs in the ifft routine. To avoid this we watch if fmultiplier=0 and if it does
    # then we average the Fourier transform around 0 to get an approximate of the limit.
    # FIXME: These list comprehensions are slow, can we speed it up? 

    ## Define function F[wT0] / (λ - i*μ*sgn(k)+ iηk + abs.(k))
    tFwT0 = (k,κ) -> pi * besselj.(0, abs.(k)) ./ fmultiplier(κ)
    if λ ≈ 0
        FwT0 = (k,κ) -> [fmultiplier(m) ≈ 0 ? (tFwT0(m-eps(),mm-eps())+tFwT0(m+eps(),mm+eps()))/2 : tFwT0(m,mm)  for (m,mm) in zip(k,κ)]
    else
        FwT0 = (k,κ) -> tFwT0(k,κ)
    end
    
    ## Define function F[wT1] / (λ - i*μ*sgn(k)+ iηk + abs.(k))
    tFwT1 = (k, κ) -> -im * pi * besselj.(1, k) ./ fmultiplier(κ)
    if λ ≈ 0
        FwT1 = (k, κ) -> [fmultiplier(m) ≈ 0 ? (tFwT1(m-eps(),mm-eps())+tFwT1(m+eps(),mm+eps()))/2 : tFwT1(m,mm)  for (m,mm) in zip(k,κ)]
    else
        FwT1 = (k,κ) -> tFwT1(k,κ)
    end
    
    ## Define function F[U0] / (λ - i*μ*sgn(k)+ iηk + abs.(k))
    tFU0 = (k,κ) -> (pi * besselj.(1, abs.(k)) + 2 .*sin.(k)./k - 2 .* sin.(abs.(k)) ./ abs.(k)) ./ fmultiplier(κ)
    FU0 = (k,κ) -> [m ≈ 0 ? (tFU0(m-eps(),mm-eps())+tFU0(m+eps(),mm+eps()))/2 : tFU0(m,mm)  for (m,mm) in zip(k,κ)]
    
    ## Define function F[U-1] / (λ - i*μ*sgn(k)+ iηk + abs.(k))
    tFU_1 = (k,κ) -> ( im*pi*k.*besselj.(0,abs.(k)) ./ abs.(k)) ./ fmultiplier(κ)
    FU_1 = (k,κ) -> [m ≈ 0 ? (tFU_1(m-eps(),mm-eps())+tFU_1(m+eps(),mm+eps()))/2 : tFU_1(m,mm)  for (m,mm) in zip(k,κ)]
    
    ## Different element supporter functions on the same sized elements are simply translations
    ## of the supporter functions on a reference element s*[-1,1] where s is the scaling.
    ## Hence, we only want to do 4 FFTs per differently sized element to compute
    ## the supporter functions on the reference elements s*[-1,1]. These are then translated in
    ## interpolate_supporter_functions. 

    s = unique(2. ./ (a[2:end] - a[1:end-1]))
    
    x = ifftshift((fftfreq(length(ω), 1/δ)) * 2 * pi) 
    
    # Compute reference supporter functions
    ywT0 = []; ywT1 = []; yU0 = []; yU_1 = []
    for ss in s
        sFwT0 = k -> (1/s).*FwT0((1/s) .*k, k)
        sFwT1 = k -> (1/s).*FwT1((1/s) .*k, k)
        sFU0 = k -> (1/s).*FU0((1/s) .*k, k)
        sFU_1 = k -> (1/s).*FU_1((1/s) .*k, k)
        append!(ywT0, [cifft(sFwT0, ω, δ, W, x)])
        append!(ywT1, [cifft(sFwT1, ω, δ, W, x)])
        append!(yU0, [cifft(sFU0, ω, δ, W, x)])
        append!(yU_1, [cifft(sFU_1, ω, δ, W, x)])
    end

    return (x, s, ywT0, yU_1, ywT1, yU0) 
end

function cifft(f, ω, δ, W, x)
    yf= ifftshift(ifft(f(ω)))
    N = length(ω)
    yf = (δ .* N .* exp.(-im .*x .*W)  ./ (2*pi)) .* yf
    return yf
end

function interpolate_supporter_functions(x, s, ywT0, yU_1, ywT1, yU0; a=[-1.,1.])

    ## Scale and translate the reference supporter functions during the interpolation
    ## for each element. 

    el_no = length(a)-1
    c = (a[2:end] - a[1:end-1]) ./ 2; d =  (a[1:end-1] + a[2:end]) ./ 2


    yU_1 = [interpolate((x .+ d[j],), real.(yU_1[findall(x->x==c[j],s)[1]])[:], Gridded(Linear())) for j in 1:el_no]
    yU0 =  [interpolate((x .+ d[j],), real.(yU0[findall(x->x==c[j],s)[1]])[:], Gridded(Linear())) for j in 1:el_no]
    ywT0 = [interpolate((x .+ d[j],), real.(ywT0[findall(x->x==c[j],s)[1]])[:], Gridded(Linear())) for j in 1:el_no]
    ywT1 = [interpolate((x .+ d[j],), real.(ywT1[findall(x->x==c[j],s)[1]])[:], Gridded(Linear())) for j in 1:el_no]
    return (ywT0, yU_1, ywT1, yU0)
end

function fft_supporter_functions(λ, μ, η; W=1000., δ=0.001, a=[-1.,1.])
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
    
    (x, s, ywT0, yU_1, ywT1, yU0) = supporter_functions(λ, μ, η, W=W, δ=δ, a=a)
    uS = interpolate_supporter_functions(x, s, ywT0, yU_1, ywT1, yU0, a=a)
    return uS
end

function coefficient_supporter_functions(A, x, uS, N; tol=1e-6)
    (ywT0, yU_1, ywT1, yU0) = uS
    el_no = length(yU0)
    yu_1 = [solvesvd(A, riemann(x, yU_1[j]); tol=tol) for j in 1:el_no]
    yu0 = [solvesvd(A, riemann(x, yU0[j]); tol=tol) for j in 1:el_no]
    ywt0 = [solvesvd(A, riemann(x, ywT0[j]); tol=tol) for j in 1:el_no]
    ywt1 = [solvesvd(A, riemann(x, ywT1[j]); tol=tol) for j in 1:el_no]
    yu_1 = [expansion_sum_space(yu_1[j],  N, el_no) for j in 1:el_no]
    yu0 = [expansion_sum_space(yu0[j], N, el_no) for j in 1:el_no]
    ywt0 = [expansion_sum_space(ywt0[j], N, el_no) for j in 1:el_no]
    ywt1 = [expansion_sum_space(ywt1[j],  N, el_no) for j in 1:el_no]
    return (ywt0, yu_1, ywt1, yu0)
end

function inverse_fourier_transform(F, ω; W=1000, δ=0.001)
    
    x = ifftshift(fftfreq(length(ω), 1/δ) * 2 * pi)
    N = length(ω)

    f = ifftshift(ifft(F(ω)))
    f = (δ .* N .* exp.(-im .*x .*W)  ./ (2*pi)) .* f

    return (x, f)
end

function fractional_heat_fourier_solve(F, ω, timesteps)
    
    Fn = k -> F(k, 1)
    (x, f) = inverse_fourier_transform(Fn, ω)
    
    fv = [f]
    for n = 2:timesteps
        Fn = k -> F(k, n)
        (_, f) = inverse_fourier_transform(Fn, ω)
        append!(fv, [f])
    end
    return (x, fv)
end