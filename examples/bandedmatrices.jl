using Revise
using SumSpaces
using BlockBandedMatrices
using LinearAlgebra
import ClassicalOrthogonalPolynomials: ∞, BlockBroadcastArray, oneto, mortar
import BlockBandedMatrices: _BandedBlockBandedMatrix,
# import BlockArrays


N = 4
rows = 2*Int32.(ones(N))
_BlockBandedMatrix(ones(10), rows, rows, (1,1))

BlockBandedMatrix(ones(8,8), rows, rows, (1,1))

BandedBlockBandedMatrix(ones(6,6), [2,2,2], [2,2,2], (2,0), (0,0))


rows = oneto(∞)
M = Ones(∞,∞)
_BandedBlockBandedMatrix(M,rows, rows, (0,0), (0,0))

rows = Fill(2,∞)
M = Ones(∞,∞)
_BandedBlockBandedMatrix(M,rows,rows, (0,0), (0,0))

_BandedBlockBandedMatrix(Ones(∞,∞), Fill(2,∞), Fill(2,∞), (2,0), (0,0))


using SumSpaces
import ClassicalOrthogonalPolynomials: ∞, mortar, Diagonal, BlockBroadcastArray