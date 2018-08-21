# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using FEMBase.FEMTest

@testset "Poisson + Dirichlet" begin
    el1 = Element(Quad4, [1, 2, 3, 4])
    el2 = Element(Tri3, [3, 2, 5])
    el3 = Element(Seg2, [3, 5])
    el4 = Element(Seg2, [1, 4])
    X = Dict(1 => [0.0, 0.0],
             2 => [1.0, 0.0],
             3 => [1.0, 1.0],
             4 => [0.0, 1.0],
             5 => [2.0, 1.0])
    update!([el1, el2, el3, el4], "geometry", X)
    update!([el1, el2], "coefficient", 6.0)
    update!(el1, "source", 132.0)
    update!(el3, "flux", 264.0)
    update!(el4, "u 1", 0.0)
    problem = Problem(Poisson, "test problem", 1)
    bc = Problem(Dirichlet, "fixed", 1, "u")
    add_elements!(problem, [el1, el2, el3])
    add_elements!(bc, [el4])

    step = Analysis(Static)
    add_problems!(step, [problem, bc])
    ls, normu, normla = run!(step)

    @test isapprox(normu, 136.5979502042399)
    @test isapprox(normla, 280.52807346146307)

end

@testset "Poisson (including bc)" begin
    el1 = Element(Quad4, [1, 2, 3, 4])
    el2 = Element(Tri3, [3, 2, 5])
    el3 = Element(Seg2, [3, 5])
    el4 = Element(Seg2, [1, 4])
    X = Dict(1 => [0.0, 0.0],
             2 => [1.0, 0.0],
             3 => [1.0, 1.0],
             4 => [0.0, 1.0],
             5 => [2.0, 1.0])
    elements = [el1, el2, el3, el4]
    update!(elements, "geometry", X)
    update!([el1, el2], "coefficient", 6.0)
    update!(el1, "source", 132.0)
    update!(el3, "flux", 264.0)
    update!(el4, "fixed u", 0.0)
    problem = Problem(Poisson, "test problem", 1)
    add_elements!(problem, elements)

    step = Analysis(Static)
    add_problems!(step, [problem])
    ls, normu, normla = run!(step)

    @test isapprox(normu, 136.5979502042399)
    @test isapprox(normla, 569.2099788303083)
end

@testset "Test @test_resource" begin
    fn = @test_resource("mesh.inp")
    @test occursin("mesh.inp", fn)
end

@testset "Test read_mtx_from_file and read_mtx_from_string" begin
    data = """
    1,1, 1,1,  1.0
    """
    K = read_mtx_from_string(data)
    @test isapprox(K[1,1], 1.0)
    K2 = read_mtx_from_file(@test_resource("eb3arxs4_STIF1.mtx"))
    @test isapprox(K2[1,1], 1.8e6)
end
