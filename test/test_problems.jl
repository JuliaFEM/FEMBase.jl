# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using FEMBase: is_field_problem, is_boundary_problem
using FEMBase: get_parent_field_name, get_global_solution
using FEMBase: get_assembly, get_elements, Field
using Base.Test

import FEMBase: get_unknown_field_dimension, get_unknown_field_name
import FEMBase: get_formulation_type

type P1 <: FieldProblem
    A :: Bool
end

function P1()
    return P1(true)
end

type P2 <: BoundaryProblem
    formulation :: Symbol
end

function P2()
    return P2(:incremental)
end

get_unknown_field_name(p::P1) = "P1"
get_unknown_field_name(p::P2) = "P2"
get_formulation_type(p::Problem{P2}) = p.properties.formulation

@testset "test creating new problems" begin
    p1 = Problem(P1, "P1", 1)
    p2 = Problem(P2, "P2", 2, "p")
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
    add_elements!(p1, [el])
    add_elements!(p2, [el])
    initialize!(p1, el, 0.0)
    initialize!(p2, el, 0.0)
    #@test haskey(el, "P1")
    #@test haskey(el, "P2")

    A = Assembly()
    A2 = Assembly()
    u = ones(2)
    la = ones(2)
    update!(p1, A, u, la)
    update!(p1, A, get_elements(p1), 0.0)
    @test length(A.u) == length(u)
    @test length(A.la) == length(la)
    @test size(get_global_solution(p1, A)[1]) == (2,)

    update!(p2, A2, ones(4), ones(4))
    update!(p2, A2, get_elements(p2), 0.0)
    u, la = get_global_solution(p2, A2)
    o = Vector{Float64}[ones(2), ones(2)]
    @test isapprox(u, o)
    @test isapprox(la, o)
    update!(p2, A2, ones(4), ones(4))
    u, la = get_global_solution(p2, A2)
    @test isapprox(u, 2*o)
    @test isapprox(la, o)
    p2.properties.formulation = :forwarddiff
    update!(p2, A2, ones(4), ones(4))
    u, la = get_global_solution(p2, A2)
    @test isapprox(u, 3*o)
    @test isapprox(la, 2*o)
    p2.properties.formulation = :total
    update!(p2, A2, ones(4), ones(4))
    u, la = get_global_solution(p2, A2)
    @test isapprox(u, o)
    @test isapprox(la, o)
    p2.properties.formulation = :nottell
    @test_throws Exception update!(p2, A2, ones(4), ones(4))

    @test length(p1) == 1
    update!(p1, "geometry", Dict(1 => [0.0, 0.0], 2 => [1.0, 0.0]))
    @test haskey(el, "geometry")
    @test isapprox(p1("geometry", 0.0)[1], [0.0, 0.0])
    @test get_parent_field_name(p2) == "p"
    @test get_gdofs(el, 1) == [1, 2]
    @test get_gdofs(p1, el) == [1, 2]
    empty!(p2)

    p3 = Problem(P2, "P3", 1, "p")
    as = get_assembly(p3)
    e1 = Element(Seg2, [1, 2])
    e2 = Element(Seg2, [2, 3])
    e3 = Element(Seg2, Int[])
    update!(e1, "f", [1.0, 1.0])
    update!(e2, "f", [2.0, 2.0])
    push!(p3, e1)
    push!(p3, [e2], [e3])
    initialize!(p3, e1, 0.0)
    initialize!(p3, e2, 0.0)
    p3.fields["f"] = Field(1.0)
    @test haskey(p3, "f")
    @test isapprox(p3["f"].data, 1.0)
    f = p3("f", 0.0)
    @test_throws Exception get_gdofs(e3, 1)
end
