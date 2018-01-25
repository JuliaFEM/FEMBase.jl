# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Base.Test

import FEMBase: solve!

type LSSolver1 <: AbstractLinearSystemSolver end
type LSSolver2 <: AbstractLinearSystemSolver end

function solve!(s::LSSolver2, ls::LinearSystem)
    fill!(ls.u, 1.0)
end

@testset "test creating new linear system solver" begin
    K = sparse(rand(3,3))
    C1 = sparse(rand(3,3))
    C2 = sparse(rand(3,3))
    D = sparse(rand(3,3))
    f = sparse(rand(3))
    g = sparse(rand(3))
    u = spzeros(3)
    la = spzeros(3)
    ls = LinearSystem(K, C1, C2, D, f, g, u, la)
    s1 = LSSolver1()
    s2 = LSSolver2()
    solve!(s1, ls)
    @test isapprox(u, spzeros(3))
    solve!(s2, ls)
    @test isapprox(u, sparse(ones(3)))
end
