# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Test

struct Dummy <: FieldProblem end

@testset "no unknown field name or assemble!-function defined" begin
    el = Element(Quad4, [1, 2, 3, 4])
    pr = Problem(Dummy, "problem", 2)
    add_elements!(pr, [el])
    assemble!(pr)
    @test true
end
