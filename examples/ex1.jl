using SumSpace, LinearAlgebra

A = Matrix([1 0;0 2])
b = Vector([1;1])
x = solvesvd(A,b)
@test x == Vector([1.0;0.5])