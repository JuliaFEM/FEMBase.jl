# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using FEMBase: assemble_mass_matrix!
using Base.Test

type Dummy <: FieldProblem end

@testset "test tet10 mass matrix" begin
    X = Dict(
             1 => [2.0, 3.0, 4.0],
             2 => [6.0, 3.0, 2.0],
             3 => [2.0, 5.0, 1.0],
             4 => [4.0, 3.0, 6.0])
    X[5] = 1/2*(X[1] + X[2])
    X[6] = 1/2*(X[2] + X[3])
    X[7] = 1/2*(X[3] + X[1])
    X[8] = 1/2*(X[1] + X[4])
    X[9] = 1/2*(X[2] + X[4])
    X[10] = 1/2*(X[3] + X[4])

    element = Element(Tet10, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    update!(element, "geometry", X)
    update!(element, "density", 105.0)

    p = Problem(Dummy, "test", 1)
    add_elements!(p, [element])
    assemble_mass_matrix!(p, 0.0)
    M = full(p.assembly.M)

    M_expected = [
        6   1   1   1  -4  -6  -4  -4  -6  -6 
        1   6   1   1  -4  -4  -6  -6  -4  -6
        1   1   6   1  -6  -4  -4  -6  -6  -4
        1   1   1   6  -6  -6  -6  -4  -4  -4
       -4  -4  -6  -6  32  16  16  16  16   8
       -6  -4  -4  -6  16  32  16   8  16  16
       -4  -6  -4  -6  16  16  32  16   8  16
       -4  -6  -6  -4  16   8  16  32  16  16
       -6  -4  -6  -4  16  16   8  16  32  16
       -6  -6  -4  -4   8  16  16  16  16  32]

    @test isapprox(M, M_expected)
end
