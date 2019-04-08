# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Test

using FEMBase: group_by_element_type
using FEMBase: get_local_coordinates, inside
using FEMBase: get_basis, get_dbasis, get_integration_order
using FEMBase: get_reference_element_coordinates, get_reference_coordinates
using FEMBase: get_element_type, is_element_type, get_element_id, filter_by_element_type

struct Dummy <: FieldProblem end

@testset "add elements to problem" begin
    problem = Problem(Dummy, "test", 2)
    element = Element(Quad4, [1, 2, 3, 4])
    elements = [element]
    add_elements!(problem, elements)
    @test problem.elements[1] == element
end

@testset "no unknown field name or assemble!-function defined" begin
    el = Element(Quad4, [1, 2, 3, 4])
    pr = Problem(Dummy, "problem", 2)
    add_elements!(pr, [el])
    assemble!(pr, 0.0)
    @test true
end

@testset "add time dependent field to element" begin
    el = Element(Seg2, [1, 2])
    u1 = ([0.0, 0.0], [0.0, 0.0])
    u2 = ([1.0, 1.0], [1.0, 1.0])
    update!(el, "displacement", 0.0 => u1)
    update!(el, "displacement", 1.0 => u2)
    @test length(el["displacement"]) == 2
    xi = (0.0, )
    ui1 = el("displacement", xi, 0.0)
    ui2 = el("displacement", xi, 0.5)
    ui3 = el("displacement", xi, 1.0)
    @debug "interpolated displacement" xi=0.0 time=0.0 u=ui1
    @debug "interpolated displacement" xi=0.0 time=0.5 u=ui2
    @debug "interpolated displacement" xi=0.0 time=1.0 u=ui3
    @test isapprox(ui1, [0.0, 0.0])
    @test isapprox(ui2, [0.5, 0.5])
    @test isapprox(ui3, [1.0, 1.0])
    el2 = Element(Poi1, [1])
    update!(el2, "force 1", 0.0 => 1.0)
end

@testset "add CVTV field to element" begin
    el = Element(Seg2, [1, 2])
    el["my field"] = (xi,t) -> xi[1]*t
    v = el("my field", [1.0], 2.0)
    @test isapprox(v, 2.0)
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
    seg2_elements = filter_by_element_type(elements, Seg2)
end

@testset "inverse isoparametric mapping" begin
    el = Element(Quad4, [1, 2, 3, 4])
    X = Dict(
        1 => [0.0, 0.0],
        2 => [1.0, 0.0],
        3 => [1.0, 1.0],
        4 => [0.0, 1.0])
    update!(el, "geometry", X)
    time = 0.0
    X1 = el("geometry", [0.1, 0.2], time)
    xi = get_local_coordinates(el, X1, time)
    X2 = el("geometry", xi, time)
    @debug "inverse isoparametric mapping" time X xi X1 X2
    @test isapprox(X1, X2)
end

@testset "inside of linear element" begin
    el = Element(Quad4, [1, 2, 3, 4])
    X = Dict(
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
    X = Dict(
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
    X = Dict(1 => [0.0, 0.0], 2 => [1.0, 0.0], 3 => [0.5, 0.5])
    f = field(X)
    #update!(el, "geometry", X)
    el["geometry"] = f
    @test isapprox(el("geometry")[1], [0.0, 0.0])
    @test isapprox(el("geometry", 0.0)[1], [0.0, 0.0])
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

@testset "gradient, jacobian and determinant of jacobian" begin
    X = Dict(
             1 => [0.0, 0.0],
             2 => [2.0, 0.0],
             3 => [2.0, 2.0],
             4 => [0.0, 2.0])
    element = Element(Quad4, [1, 2, 3, 4])
    element.fields["geometry"] = field(X)
    J = element([0.0, 0.0], 0.0, Val{:Jacobian})
    detJ = element([0.0, 0.0], 0.0, Val{:detJ})
    @test isapprox(det(J), detJ)
    @test isapprox(detJ, 1.0)
    element2 = Element(Seg2, [1, 2])
    update!(element2, "geometry", X)
    detJ = element2([0.0], 0.0, Val{:detJ})
    @test isapprox(detJ, 1.0)
    X = Dict(
             1 => [0.0, 0.0, 0.0],
             2 => [2.0, 0.0, 0.0],
             3 => [2.0, 2.0, 0.0],
             4 => [0.0, 2.0, 0.0])
    element3 = Element(Quad4, [1, 2, 3, 4])
    update!(element3, "geometry", X)
    detJ = element3([0.0, 0.0], 0.0, Val{:detJ})
    @test isapprox(detJ, 1.0)
    grad = element([0.0, 0.0], 0.0, Val{:Grad})
    @test isapprox(grad, 1/4*[-1.0 1.0 1.0 -1.0; -1.0 -1.0 1.0 1.0])
    update!(element3, "youngs modulus", 1.0)
    @test isapprox(element3("youngs modulus", [0.0], 0.0), 1.0)
    update!(element3, "test", [1.0, 2.0, 3.0])
    test_val = element3("test", [0.0], 1.0)
end

@testset "update data to element" begin
    element = Element(Seg2, [1, 2])
    elements = Element[element]

#=
    # constant fields (not changing on element area)

    # time intependent

    # scalar + updating scalar
    update!(elements, "scalar 1", 0.0)
    @test isapprox(element("scalar 1", (0.0), 0.0), 0.0)
    update!(elements, "scalar 1", 1.0)
    @test isapprox(element("scalar 1", (0.0), 0.0), 1.0)
    
    # vector + updating vector
    update!(elements, "vector 1", [1.0, 2.0])
    @test isapprox(element("vector 1", (0.0), 0.0), [1.0, 2.0])
    update!(elements, "vector 1", [2.0, 1.0])
    @test isapprox(element("vector 1", (0.0), 0.0), [2.0, 1.0])

    # tensor + updating tensor
    update!(elements, "tensor 1", [1.0 2.0; 3.0 4.0])
    @test isapprox(element("tensor 1", [1.0 2.0; 3.0 4.0]))
    update!(elements, "tensor 1", [4.0 3.0; 2.0 1.0])
    @test isapprox(element("tensor 1", [4.0 3.0; 2.0 1.0]))

    # time dependent

    # scalar + interpolate
    update!(elements, "scalar 2", 0.0 => 0.0)
    @test isapprox(elements("scalar 2", (0.0), 0.5), 0.0)
    update!(elements, "scalar 2", 0.0 => 1.0)
    @test isapprox(elements("scalar 2", (0.0), 0.5), 1.0)
    update!(elements, "scalar 2", 1.0 => 0.0)
    @test isapprox(elements("scalar 2", (0.0), 0.5), 0.5)

    # vector + interpolate
    update!(elements, "vector 2", 0.0 => [0.0, 0.0])
    @test isapprox(elements("vector 2", (0.0), 0.5), [0.0, 0.0])
    update!(elements, "vector 2", 0.0 => [0.0, 1.0])
    @test isapprox(elements("vector 2", (0.0), 0.5), [0.0, 1.0])
    update!(elements, "vector 2", 1.0 => [1.0, 0.0])
    @test isapprox(elements("vector 2", (0.0), 0.5), [0.5, 0.5])

    # tensor + interpolate
    update!(elements, "tensor 2", 0.0 => [0.0 0.0; 0.0 0.0])
    @test isapprox(elements("tensor 2", (0.0), 0.5), [0.0 0.0; 0.0 0.0])
    update!(elements, "tensor 2", 0.0 => [1.0 1.0; 1.0 1.0])
    @test isapprox(elements("tensor 2", (0.0), 0.5), [1.0 1.0; 1.0 1.0])
    update!(elements, "tensor 2", 1.0 => [2.0 2.0; 2.0 2.0])
    @test isapprox(elements("tensor 2", (0.0), 0.5), [1.5 1.5; 1.5 1.5])

    # variable fields (changing on element area, needs interpolating)

    # time intependent

    # scalar changing on element area
    update!(elements, "scalar 3", (1.0, 2.0))
    @test isapprox(element("scalar 3", (0.0), 0.0), 1.5)
    update!(elements, "scalar 3", (2.0, 3.0))
    @test isapprox(element("scalar 3", (0.0), 0.0), 2.5)

    # vector changing on element area
    update!(elements, "vector 3", ([1.0, 2.0], [2.0, 3.0]))
    @test isapprox(element("vector 3", (0.0), 0.0), [1.5, 2.5])
    update!(elements, "vector 3", ([2.0, 3.0], [3.0, 4.0]))
    @test isapprox(element("vector 3", (0.0), 0.0), [2.5, 3.5])

    # tensor changing on element area
    update!(elements, "tensor 3", ([0.0 0.0; 0.0 0.0], [1.0 1.0; 1.0 1.0]))
    @test isapprox(element("tensor 3", (0.0), 0.0), [0.5 0.5; 0.5 0.5])
    update!(elements, "tensor 3", ([1.0 1.0; 1.0 1.0], [2.0 2.0; 2.0 2.0]))
    @test isapprox(element("tensor 3", (0.0), 0.0), [1.5 1.5; 1.5 1.5])

    # time dependent

    # scalar changing on element area over time
    update!(elements, "scalar 4", 0.0 => (1.0, 3.0))
    update!(elements, "scalar 4", 1.0 => (3.0, 5.0))
    @test isapprox(element("scalar 4", (0.0), 0.5), 3.0)

    # vector changing on element area over time
    update!(elements, "vector 4", 0.0 => ([1.0, 1.0], [3.0, 3.0]))
    update!(elements, "vector 4", 1.0 => ([3.0, 3.0], [5.0, 5.0]))
    @test isapprox(element("scalar 4", (0.0), 0.5), [3.0, 3.0])

    # tensor changing on element area over time
    update!(elements, "tensor 4", 0.0 => ([1.0 1.0; 1.0 1.0], [3.0 3.0; 3.0 3.0]))
    update!(elements, "tensor 4", 1.0 => ([3.0 3.0; 3.0 3.0], [5.0 5.0; 5.0 5.0]))
    @test isapprox(element("tensor 4", (0.0), 0.5), [3.0 3.0; 3.0 3.0])
=#

end

@testset "calculate displacement gradient" begin
    X = Dict(
             1 => [0.0, 0.0],
             2 => [1.0, 0.0],
             3 => [1.0, 1.0],
             4 => [0.0, 1.0])
    u = Dict(
             1 => [0.0, 0.0],
             2 => [1.0, -1.0],
             3 => [2.0, 3.0],
             4 => [0.0, 0.0])
    element = Element(Quad4, [1, 2, 3, 4])
    update!(element, "geometry", X)
    update!(element, "displacement", 0.0 => u)
    gradu = element("displacement", (0.0, 0.0), 0.0, Val{:Grad})
    gradu_expected = [1.5 0.5; 1.0 2.0]
    @test isapprox(gradu, gradu_expected)

    element = Element(Quad4, [1, 2, 3, 4])
    X = ([0.0,0.0], [1.0,0.0], [1.0,1.0], [0.0,1.0])
    u = ([0.0,0.0], [1.0,-1.0], [2.0,3.0], [0.0,0.0])
    element["geometry"] = X
    element["displacement"] = u
    gradu = element("displacement", (0.0, 0.0), 0.0, Val{:Grad})
    gradu_expected = [1.5 0.5; 1.0 2.0]
    @test isapprox(gradu, gradu_expected)
    @test isapprox(element("geometry", (0.0,0.0), 0.0), [0.5,0.5])
end

@testset "interpolate field" begin
    element = Element(Quad4, [1,2,3,4])
    update!(element, "test", 0.0 => 0.0)
    update!(element, "test", 1.0 => 1.0)
    @test isapprox(element("test", (0.0, 0.0), 0.5), 0.5)
end

@testset "get integration points for Quad4" begin
    el = Element(Quad4, [1, 2, 3, 4])
    ips = get_integration_points(el)
    @test length(ips) == 4
end

@testset "test nonconverging inverse isoparametric mapping" begin
    X = Dict(1 => [0.0, 0.0], 2 => [1.0, 0.0], 3 => [1.0, 1.0], 4 => [0.0, 1.0])
    el = Element(Quad4, [1, 2, 3, 4])
    update!(el, "geometry", X)
    @test_throws Exception get_local_coordinates(el, [0.1, 0.1], 0.0; max_iterations=0)
end

@testset "analytical functions as fields" begin
    f(xi, time) = xi[1]*time
    g(element, xi, time) = element("geometry", xi, time)*time
    X = Dict(1 => [0.0, 0.0], 2 => [1.0, 0.0], 3 => [1.0, 1.0], 4 => [0.0, 1.0])
    el = Element(Quad4, [1, 2, 3, 4])
    update!(el, "geometry", X)
    update!(el, "f", f)
    update!(el, "g", g)
    @test isapprox(el("f", (0.5, 0.5), 2.0), 1.0)
    @test isapprox(el("g", (0.0, 0.0), 2.0), [1.0, 1.0])
end

@testset "add field to element using update!" begin
    el = Element(Seg2, [1, 2])
    f = field(1)
    update!(el, "f", f)
    @test f === el.fields["f"]
end

@testset "update DVTI/DVTV with Dict" begin
    el = Element(Seg2, [1, 2])
    update!(el, "u", 0.0 => (0.0, 0.0))
    update!(el, "v", (0.0, 0.0))
    u = Dict(1 => 1.0, 2 => 1.0)
    update!(el, "u", 0.0 => u)
    update!(el, "v", u)
    @test isapprox(el("u", (0.0,), 0.0), 1.0)
    @test isapprox(el("v", (0.0,), 0.0), 1.0)
end
