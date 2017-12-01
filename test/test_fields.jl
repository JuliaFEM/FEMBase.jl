# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Base.Test

@testset "DCTI field" begin

    # scalar field
    a = DCTI(1)
    @test interpolate(a, 0.0) == 1
    @test a[1] == 1
    @test length(a) == length(1)
    @test size(a) == size(1)
    update!(a, 2)
    @test a == 2

    # vector field
    b = DCTI([1, 2])
    @test interpolate(b, 0.0) == [1, 2]
    @test b[1] == [1, 2]
    @test length(b) == length([1, 2])
    @test size(b) == size([1, 2])
    update!(b, [2, 3])
    @test b == [2, 3]

    # tensor field
    c = DCTI([1 2; 3 4])
    @test interpolate(c, 0.0) == [1 2; 3 4]
    @test c[1] == [1 2; 3 4]
    @test length(c) == length([1 2; 3 4])
    @test size(c) == size([1 2; 3 4])
    update!(c, [2 3; 4 5])
    @test c == [2 3; 4 5]

end

@testset "DVTI field" begin

    # scalar field
    a = DVTI((1, 2))
    @test a[1] == 1
    @test a[2] == 2
    @test interpolate(a, 0.0) == (1, 2)
    @test interpolate(a, 0.0, [1, 1]) == 3
    update!(a, (2, 3))
    @test a == (2, 3)
end
