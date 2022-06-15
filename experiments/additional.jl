using Revise
using SumSpacesBeta
using MathLink


λ=1; η=0; μ=0; N=31;

y = 1.1;
@time a1 = weval(W`Re[1/(2*Pi)*NIntegrate[(-I)^(N+2)*Pi* BesselJ[N+2,k]/(λ + η*I*k - μ*I*Sign[k] + Abs[k]) 
* Exp[I y k], {k,-∞,∞}, MaxRecursion -> 50, WorkingPrecision -> 15, PrecisionGoal -> 12]]`;λ=λ,η=η,μ=μ,y=y,N=N)



d = x-> -(Base.MathConstants.eulergamma + log.(abs.(x)))/π
v = x-> ExtendedWeightedChebyshevT[x,1]

@time a2 = weval(W`NIntegrate[
    Which[Abs[x] < 1, ChebyshevT[N+2,x]/Sqrt[1 - x^2], True, 0]*(
        (E^(-I a (y-x)) (-I (-1 + E^(2 I a (y-x))) \[Pi] (y-x) + 
        Abs[(y-x)] (-2 (1 + E^(2 I a (y-x))) CosIntegral[a (y-x)] + 
           E^(2 I a (y-x)) Log[1/(y-x)^2] + 2 Log[(y-x)] + 2 E^(2 I a (y-x)) Log[(y-x)] - 
           2 Log[Abs[(y-x)]] + 
           2 I (-1 + E^(2 I a (y-x))) SinIntegral[a (y-x)])))/(4 \[Pi] Abs[(y-x)]))
    , {x, -Infinity, Infinity}, MaxRecursion->100, WorkingPrecision -> 15, 
    PrecisionGoal -> 12]`;a=λ, y=y,N=N)