using ClassicalOrthogonalPolynomials
using SumSpacesBeta
using HypergeometricFunctions, SpecialFunctions
using Plots

J = Jacobi(1/2,1/2)
wJ = n -> 2^(2(n-1)+1)/binomial(2(n-1)+2,n)

U = ChebyshevU()
eU = ExtendedChebyshevU()

α=1/2; β=1/2;
hJ = (x, n) -> gamma(α+1+(n-1))/(gamma(n)*gamma(α+1)) .* _₂F₁.(-(n-1), 1+α+β+(n-1), α+1, (1 .- x)/2)

xx = -4:0.01:4
n = 3;
plot(xx, eU[xx,n+2])
# plot!(xx, wJ(n).*J[xx,n])
plot!(xx, wJ(n).*real.(hJ(ComplexF64.(xx),n)))


J = Jacobi(-1/2,-1/2)
wJ = n -> 2^(2(n-1))/binomial(2(n-1),n-1)
eT = ExtendedChebyshevT()
α=-1/2; β=-1/2;
hJ = (x, n) -> gamma(α+1+(n-1))/(gamma(n)*gamma(α+1)) .* _₂F₁.(-(n-1), 1+α+β+(n-1), α+1, abs.(x))
xx = -1.:0.01:1.
n = 3;
plot(xx, eT[xx,n])
plot!(xx, wJ(n).*(hJ(xx,n)))