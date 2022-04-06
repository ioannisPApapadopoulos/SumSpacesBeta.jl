using Revise
using SumSpaces, Plots
import LinearAlgebra: I

N = 5 # Truncation degree
λ = 1e2; μ = 0; Δt = 1/λ # Constants

Sp = SumSpaceP() # Primal sum space
Sd = SumSpaceD() # Dual sum space

a = [-3,-1,1.,3]
# a = [-1.,2.]
el_no = length(a)-1


M = max(N^2,5001)  # Number of collocation points in [-1,1]
Me = M ÷ 10  # Number of collocation points in [-2,-1) and (1,2].
x = collocation_points(M, Me, innergap=1e-5, endpoints=10) # Collocation points
x1 = affinetransform(a[1],a[2],x)
x2 = affinetransform(a[2],a[3],x)
x3 = affinetransform(a[3],a[4],x)


Nn = min(N,11)

A = framematrix([x1,x2,x3], Sp, Nn, M, Me) # Blocked frame matrix

# Helper functions for computing support functions
(w, yU_1, yU0, ywT0, ywT1) = supporter_functions(λ, μ, a=a)
(yU_1, yU0, ywT0, ywT1) = interpolate_supporter_functions(w, yU_1, yU0, ywT0, ywT1)
(yu_1, yu0, ywt0, ywt1) = columns_supporter_functions(A, x, yU_1, yU0, ywT0, ywT1, Nn, N)

xx = -5:0.01:5
plot(xx, yU0[1](xx))
# plot!(xx, yU0[2](xx))
eSp = ElementSumSpace{1}(a)
eSd = ElementSumSpace{2}(a)

y = eSp[xx,1:length(yu0[1])]*yu0[1]
plot!(xx, y)

# Create appended sum space
ASp = ElementAppendedSumSpace((ywT0, yU_1, ywT1, yU0), (ywt0, yu_1, ywt1, yu0), a)

# Create matrix for element 1
Id = (eSd \ ASp)[1:2N+7+(el_no-1)*(2N+6),1:2N+7+(el_no-1)*(2N+6)]

x = axes(Sp, 1); H = inv.( x .- x')
Hm = (1/π).*((eSp\(H*eSp))[1:2N+3,1:2N+3])    # Hilbert: Sp -> Sp
Cm = [(eSd\(Derivative(x)*eSp)[j])[1:2N+7,1:2N+3] for j in 1:el_no]# Derivative: Sp -> Sd
Bm = (eSd\eSp)[1:2N+7,1:2N+3]                 # Identity: Sp -> Sd


Dm =  [λ.*Bm + μ.*Bm*Hm + Cm[j]*Hm for j in 1:el_no]     # Helmholtz-like operator: Sp -> Sd   
Dm = [hcat(zeros(size(Dm[j],1), 4),Dm[j]) for j in 1:el_no] # Adding 4 columns to construct: ASp -> Sd
for j in 1:el_no
    Dm[j][2:5,1:4] = I[1:4,1:4]
    if j == 1
        Dm[j] = [Dm[j][:,5] Dm[j][:,1:4] Dm[j][:,6:end]] # permute T0 column to start
    else
        Dm[j] = [Dm[j][:,1:4] Dm[j][:,6:end]]
        Dm[j] = Dm[j][2:end,:] 
    end
end
# Initial condition, u₀ = √(1-x²)
u₀ = zeros(1+el_no*(2N+6))
u₀[6*el_no-3] = 1.
u = [u₀]

# Run solve loop for time-stepping
timesteps=50
for k = 1:timesteps
    u1 = []
    
    # Map from ASp to Sd
    v = Id * u[k]
    # Multiply RHS with λ
    v = λ.*v
    v = BlockArray(v, vcat(1,Fill(el_no,(length(v)-1)÷el_no)))
    
    # Rearrange coefficients element-wise
    v = coefficient_stack(v, N, el_no, append=true)

    for j = 1:el_no
        # Solve for each element seperately
        append!(u1, Dm[j]\v[Block.(j)])
     end

    # Rearrange coefficients back to interlaced
    u1 = coefficient_interlace(u1, N, el_no, append=true)
    
    # Append solution
    append!(u,  [u1])

    # Normalise constant so that it decays to zero
    u[k+1][1] = u[k+1][1] - ASp[1e2,1:length(u[k+1])]'*u[k+1]

end


# Plot solution
p = plot()
xx = -5:0.001:5
xlim = [xx[1],xx[end]]; ylim = [-0.1,1]
for k = 1: timesteps+1
    t = round(Δt*(k-1), digits=2)
    yy = ASp[xx,1:length(u[k])]*u[k]
    p = plot(xx,yy, label="time=$t (s)", legend=:topleft, xlim=xlim, ylim=ylim)
    sleep(0.1)
    display(p)
end  