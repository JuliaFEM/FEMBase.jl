# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using FEMBase: is_field_problem, is_boundary_problem, get_formulation_type
using FEMBase: get_parent_field_name, get_global_solution
using FEMBase: get_assembly
using Base.Test

import FEMBase: get_unknown_field_dimension, get_unknown_field_name

type P1 <: FieldProblem
    A :: Bool
end

function P1()
    return P1(true)
end

type P2 <: BoundaryProblem end

get_unknown_field_name(p::P1) = "P1"
get_unknown_field_dimension(p::P1) = 1
get_unknown_field_name(p::P2) = "P2"
get_unknown_field_dimension(p::P2) = 2

@testset "test creating new problems" begin
    p1 = Problem(P1, "P1", 1)
    p2 = Problem(P2, "P2", 1, "p")
    @test get_formulation_type(p1) == :incremental
    @test get_unknown_field_name(p2) == "lambda"
    @test isempty(get_assembly(p1))
    @test is_field_problem(p1)
    @test !is_field_problem(p2)
    @test !is_boundary_problem(p1)
    @test is_boundary_problem(p2)
    @test p1.properties.A
    update!(p1.properties, "A" => "false")
    @test !p1.properties.A
    el = Element(Seg2, [1, 2])
    initialize!(p1, el, 0.0)
    initialize!(p2, el, 0.0)
    #@test haskey(el, "P1")
    #@test haskey(el, "P2")
    A = Assembly()
    u = ones(2)
    la = ones(2)
    update!(p1, A, u, la)
    @test length(A.u) == length(u)
    @test length(A.la) == length(la)
    @test size(get_global_solution(p1, A)[1]) == (2,)
    @test size(get_global_solution(p2, A)[1]) == (2,)
    add_elements!(p1, Element[el])
    @test length(p1) == 1
    update!(p1, "geometry", Dict(1 => [0.0, 0.0], 2 => [1.0, 0.0]))
    @test haskey(el, "geometry")
    @test isapprox(p1("geometry", 0.0)[1], [0.0, 0.0])
    @test get_parent_field_name(p2) == "p"
    @test get_gdofs(el, 1) == [1, 2]
    @test get_gdofs(p1, el) == [1, 2]
    empty!(p2)
end
