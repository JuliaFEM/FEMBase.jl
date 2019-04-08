# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase, Test, SparseArrays
import FEMBase: solve!, can_solve

mutable struct LSSolver1 <: AbstractLinearSystemSolver
    a :: Int
end

mutable struct LSSolver2 <: AbstractLinearSystemSolver
end

function LSSolver1()
    return LSSolver1(1)
end

function can_solve(ls::LinearSystem, solver::LSSolver1)
    return (false, "don't know how to do")
end

function solve!(ls::LinearSystem, solver::LSSolver2)
    fill!(ls.u, 1.0)
end

@testset "test creating new linear system solver" begin
    ls = LinearSystem(3)
    ls.u = spzeros(3)
    solver1 = LSSolver1()
    solver2 = LSSolver2()
    solve!(ls, solver1)
    @test can_solve(ls, solver1)[1] == false
    solve!(ls, [solver1, solver2])
    @test isapprox(ls.u, ones(3))
    @test_throws ErrorException solve!(ls, [solver1])
end
