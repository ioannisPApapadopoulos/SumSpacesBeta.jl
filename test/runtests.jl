using Test, SumSpaces


"""
Test functions in frame.jl
"""
A = Matrix([1. 0;0 2])
b = Vector([1.;1])
x = solvesvd(A,b,block=false)
@test x == Vector([1.0;0.5])

x = collocation_points(3, 3, endpoints=2, outergap=0.)
@test x ≈ Vector([-1.;0.;1.;-2;-1.5;-1;1;1.5;2])


