# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Test

struct Dummy <: FEMBase.FieldProblem end

@testset "add elements to problem" begin
    problem = Problem(Dummy, "test", 2)
    element = Element(Quad4, [1, 2, 3, 4])
    elements = [element]
    add_elements!(problem, elements)
    @test problem.elements[1] == element
end
