# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Base.Test

using FEMBase: Field
using FEMBase: group_by_element_type
using FEMBase: get_local_coordinates, inside
using FEMBase: get_basis, get_dbasis, get_integration_order
using FEMBase: get_reference_element_coordinates, get_reference_coordinates

@testset "add time dependent field to element" begin
    el = Element(Seg2, [1, 2])
    u1 = Vector{Float64}[[0.0, 0.0], [0.0, 0.0]]
    u2 = Vector{Float64}[[1.0, 1.0], [1.0, 1.0]]
    update!(el, "displacement", 0.0 => u1)
    update!(el, "displacement", 1.0 => u2)
    @test length(el["displacement"]) == 2
    @test isapprox(el("displacement", [0.0], 0.0), [0.0, 0.0])
    @test isapprox(el("displacement", [0.0], 0.5), [0.5, 0.5])
    @test isapprox(el("displacement", [0.0], 1.0), [1.0, 1.0])
    el2 = Element(Poi1, [1])
    update!(el2, "force 1", 0.0 => 1.0)
end

@testset "add CVTV field to element" begin
    el = Element(Seg2, [1, 2])
    f(xi, time) = xi[1]*time
    update!(el, "my field", f)
    v = el("my field", [1.0], 2.0)
    @test isapprox(v, 2.0)
end

@testset "add DCTI to element" begin
    el = Element(Quad4, [1, 2, 3, 4])
    update!(el, "displacement load", DCTI([4.0, 8.0]))
    @test isa(el["displacement load"], DCTI)
    @test !isa(el["displacement load"].data, DCTI)
    update!(el, "displacement load 2", [4.0, 8.0])
    @test isa(el["displacement load 2"], DCTI)
    update!(el, "temperature", [1.0, 2.0, 3.0, 4.0])
    @test isa(el["temperature"], DVTI)
    @test isapprox(el("displacement load", [0.0, 0.0], 0.0), [4.0, 8.0])
end

@testset "interpolate DCTI from element" begin
    el = Element(Seg2, [1, 2])
    update!(el, "foobar", 1.0)
    fb = el("foobar", [0.0], 0.0)
    @test isa(fb, Float64)
    @test isapprox(fb, 1.0)
end

@testset "add elements to elements" begin
    el1 = Element(Seg2, [1, 2])
    el2 = Element(Seg2, [3, 4])
    update!(el1, "master elements", [el2])
    lst = el1("master elements", 0.0)
    @test isa(lst, Vector)
end

@testset "extend basis" begin
    el = Element(Quad4, [1, 2, 3, 4])
    expected = [
        0.25 0.00 0.25 0.00 0.25 0.00 0.25 0.00
        0.00 0.25 0.00 0.25 0.00 0.25 0.00 0.25]
    @test isapprox(el([0.0, 0.0], 0.0, 2), expected)
end

@testset "group elements" begin
    e1 = Element(Seg2, [1, 2])
    e2 = Element(Quad4, [1, 2, 3, 4])
    elements = [e1, e2]
    r = group_by_element_type(elements)
    @test length(r) == 2
    @test first(r[Element{Seg2}]) == e1
    @test first(r[Element{Quad4}]) == e2
    @test get_element_type(e1) == Seg2
    e1.id = 1
    @test get_element_id(e1) == 1
    @test is_element_type(e1, Seg2)
    @test length(filter_by_element_type(elements, Seg2)) == 1
end

@testset "inverse isoparametric mapping" begin
    el = Element(Quad4, [1, 2, 3, 4])
    X = Dict{Int64, Vector{Float64}}(
        1 => [0.0, 0.0],
        2 => [1.0, 0.0],
        3 => [1.0, 1.0],
        4 => [0.0, 1.0])
    update!(el, "geometry", X)
    time = 0.0
    X1 = el("geometry", [0.1, 0.2], time)
    xi = get_local_coordinates(el, X1, time)
    X2 = el("geometry", xi, time)
    info("X1 = $X1, X2 = $X2")
    @test isapprox(X1, X2)
end

@testset "inside of linear element" begin
    el = Element(Quad4, [1, 2, 3, 4])
    X = Dict{Int64, Vector{Float64}}(
        1 => [0.0, 0.0],
        2 => [1.0, 0.0],
        3 => [1.0, 1.0],
        4 => [0.0, 1.0])
    update!(el, "geometry", X)
    time = 0.0
    @test inside(el, [0.5, 0.5], time) == true
    @test inside(el, [1.0, 0.5], time) == true
    @test inside(el, [1.0, 1.0], time) == true
    @test inside(el, [1.01, 1.0], time) == false
    @test inside(el, [1.0, 1.01], time) == false
end

@testset "inside of quadratic element" begin
    el = Element(Tri6, [1, 2, 3, 4, 5, 6])
    X = Dict{Int64, Vector{Float64}}(
        1 => [0.0, 0.0],
        2 => [1.0, 0.0],
        3 => [0.0, 1.0],
        4 => [0.5, 0.2],
        5 => [0.8, 0.6],
        6 => [-0.2, 0.5])
    update!(el, "geometry", X)
    p = [0.94, 0.3] # visually checked to be inside
    @test inside(el, p, 0.0) == true
    p = [-0.2, 0.8] # visually checked to be outside
    @test inside(el, p, 0.0) == false
end

@testset "dict field" begin
    el = Element(Seg2, [1, 2])
    X = Dict{Int64, Vector{Float64}}(1 => [0.0, 0.0], 2 => [1.0, 0.0], 3 => [0.5, 0.5])
    f = Field(X)
    debug("field = $f")
    #update!(el, "geometry", X)
    el["geometry"] = f
    @test isapprox(el("geometry")[1], [0.0, 0.0])
    @test isapprox(el("geometry", 0.0)[1], [0.0, 0.0])
    @test isapprox(el("geometry", 0.0)[3], [0.5, 0.5])
    @test isapprox(el("geometry", [0.0], 0.0), [0.5, 0.0])
end

@testset "test nodal element" begin
    el = Element(Poi1, [1])
    @test get_basis(el, [0.0], 1.0) == [1]
    @test get_dbasis(el, [0.0], 1.0) == [0]
    @test isapprox(el([0.0], 0.0, Val{:detJ}), 1.0)
    @test get_integration_order(el.properties) == 1
    ips = get_integration_points(el.properties, 1)
    @test length(ips) == 1
    @test length(ips[1]) == 2
    @test size(el) == (0, 1)
    @test length(el) == 1
    @test get_reference_element_coordinates(Poi1) == Vector{Float64}[[0.0]]
end

@testset "get reference coordinates" begin
    el = Element(Seg2, [1, 2])
    X = get_reference_coordinates(el)
    @test isapprox(X[1][1], -1.0)
    @test isapprox(X[2][1],  1.0)
end

@testset "function field" begin
    function f(element, xi, time)
        return xi[1]*time
    end
    e1 = Element(Seg2, [1, 2])
    e1["f"] = f
    @test isapprox(e1("f", [1.0], 1.0), 1.0)
end

@testset "index" begin
    e1 = Element(Seg2, [1, 2])
    update!(e1, "f", 0.0 => 1.0)
    update!(e1, "f", 1.0 => 2.0)
    @test isapprox(last(e1, "f"), 2.0)
end
