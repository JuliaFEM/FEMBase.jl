# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using FEMBase: is_field_problem, is_boundary_problem
using FEMBase: get_parent_field_name, get_global_solution
using FEMBase: get_assembly, get_elements
using Test

import FEMBase: get_unknown_field_dimension, get_unknown_field_name
import FEMBase: get_formulation_type, assemble_elements!

mutable struct P1 <: FieldProblem
    A :: Bool
end

function P1()
    return P1(true)
end

mutable struct P2 <: BoundaryProblem
    formulation :: Symbol
end

function P2()
    return P2(:incremental)
end

mutable struct P3 <: FieldProblem
end

get_unknown_field_name(p::P1) = "P1"
get_formulation_type(p::Problem{P2}) = p.properties.formulation

@testset "test get_unknown_field_name" begin
    p = Problem(P3, "P3", 1)
    @test get_unknown_field_name(p) == "N/A"
end

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
    @debug "p1 geometry field" X = p1("geometry", 0.0)
    @debug "e1 geometry field" X = el("geometry", 0.0)
    @test isapprox(p1("geometry", 0.0)[1], [0.0, 0.0])
    @test get_parent_field_name(p2) == "p"
    @test get_gdofs(p1, el) == [1, 2]

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
    p3.fields["f"] = field(1.0)
    @test haskey(p3, "f")
    @test isapprox(p3["f"].data, 1.0)
    f = p3("f", 0.0)
    @test_throws Exception get_gdofs(problem, e3)
end

@testset "test initializing new problems" begin
    p1 = Problem(P1, "P1", 2)
    p2 = Problem(P2, "P2", 2, "P1")
    @test get_unknown_field_name(p1) == "P1"
    @test get_unknown_field_name(p2) == "lambda"
    @test get_parent_field_name(p2) == "P1"
    e1 = Element(Seg2, [1, 2])
    e2 = Element(Seg2, [1, 2])
    add_elements!(p1, [e1])
    add_elements!(p2, [e2])
    initialize!(p1, 0.0)
    initialize!(p2, 0.0)
    @debug "e1 keys" keys(e1.fields)
    @test isapprox(e1("P1", (0.0,), 0.0), [0.0, 0.0])
    @test isapprox(e2("P1", (0.0,), 0.0), [0.0, 0.0])
    @test isapprox(e2("lambda", (0.0,), 0.0), [0.0, 0.0])
end

function assemble_elements!(problem::Problem{P1},
                            assembly::Assembly,
                            elements::Vector{Element{E}},
                            time::Float64) where E

    @debug "assemble_elements!" eltype=E
    bi = BasisInfo(E)
    ndofs = length(E)
    Ke = zeros(ndofs, ndofs)
    K = assembly.K

    for element in elements
        fill!(Ke, 0.0)
        for ip in get_integration_points(element)
            J, detJ, N, dN = element_info!(bi, element, ip, time)
            c = element("coefficient", ip, time)
            Ke += ip.weight * c*dN'*dN * detJ
        end
        gdofs = get_gdofs(problem, element)
        add!(K, gdofs, gdofs, Ke)
    end

    return nothing

end

@testset "test assemble field problem" begin
    el1 = Element(Quad4, [1, 2, 3, 4])
    el2 = Element(Tri3, [3, 2, 5])
    X = Dict(1 => [0.0, 0.0],
             2 => [1.0, 0.0],
             3 => [1.0, 1.0],
             4 => [0.0, 1.0],
             5 => [2.0, 1.0])
    elements = [el1, el2]
    update!(elements, "geometry", X)
    update!(elements, "coefficient", 6.0)

    problem = Problem(P1, "test problem", 1)
    add_elements!(problem, elements)
    time = 0.0
    assemble!(problem, time)
    K_expected = [
                  4.0 -1.0 -2.0 -1.0  0.0
                 -1.0  7.0 -4.0 -2.0  0.0
                 -2.0 -4.0 10.0 -1.0 -3.0
                 -1.0 -2.0 -1.0  4.0  0.0
                  0.0  0.0 -3.0  0.0  3.0]
    @test isapprox(problem.assembly.K, K_expected)
end

mutable struct DirBC <: BoundaryProblem
end

function assemble_elements!(problem::Problem{DirBC},
                            assembly::Assembly,
                            elements::Vector{Element{E}},
                            time::Float64) where E

    name = get_parent_field_name(problem)
    dim = get_unknown_field_dimension(problem)

    data = Dict{Int,Float64}()
    for element in elements
        for i=1:dim
            haskey(element, "$name $dim") || continue
            gdofs = get_gdofs(problem, element)
            ldofs = gdofs[i:dim:end]
            xis = get_reference_element_coordinates(E)
            for (ldof, xi) in zip(ldofs, xis)
                data[ldof] = interpolate(element, "$name $dim", xi, time)
            end
        end
    end

    for (dof, val) in data
        add!(assembly.C1, dof, dof, 1.0)
        add!(assembly.C2, dof, dof, 1.0)
        add!(assembly.g, dof, val)
    end

    return nothing

end

@testset "test assemble boundary problem" begin
    el1 = Element(Seg2, [1, 2])
    el2 = Element(Seg2, [2, 3])
    X = Dict(1 => [0.0, 0.0],
             2 => [1.0, 0.0],
             3 => [1.0, 1.0],
             4 => [0.0, 1.0],
             5 => [2.0, 1.0])
    elements = [el1, el2]
    update!(elements, "geometry", X)
    T0 = Dict(1 => 0.5, 2 => 0.0, 3 => 1.0)
    update!(elements, "temperature 1", T0)
    problem = Problem(DirBC, "boundary", 1, "temperature")
    add_elements!(problem, elements)
    time = 0.0
    assemble!(problem, time)
    C1_expected = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    g_expected = [0.5, 0.0, 1.0]
    @test isapprox(problem.assembly.C1, C1_expected)
    @test isapprox(problem.assembly.g, g_expected)
end

@testset "test deprecated assembly procedure" begin
    problem = Problem(P3, "P3", 1)
    elements = [Element(Seg2, [1, 2])]
    add_elements!(problem, elements)
    assemble!(problem)
    assemble!(problem.assembly, problem, first(elements), 0.0)
    assemble!(problem.assembly, problem, problem.elements, 0.0)
    assemble_elements!(problem, problem.assembly, elements, 0.0)
end

@testset "test set and get global dofs for element" begin
    element = Element(Seg2, [1, 2])
    problem = Problem(P3, "P3", 2)
    @test get_gdofs(problem, element) == [1, 2, 3, 4]
    set_gdofs!(problem, element, [2, 3, 4, 5])
    @test get_gdofs(problem, element) == [2, 3, 4, 5]
    element2 = Element(Seg2, Int[])
    @test_throws ErrorException get_gdofs(problem, element2)
end
