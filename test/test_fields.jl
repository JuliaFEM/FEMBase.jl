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
    b = DCTI([1,2])
    @test interpolate(b, 0.0) == [1,2]
    @test b[1] == [1,2]
    @test length(b) == length([1,2])
    @test size(b) == size([1,2])
    update!(b, [2,3])
    @test b == [2,3]

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
    @test interpolate(a, 0.0, [1,1]) == 3
    update!(a, (2,3))
    @test a == (2,3)

    # vector field
    b = DVTI(([1,2], [2,3]))
    @test b[1] == [1,2]
    @test b[2] == [2,3]
    @test interpolate(b, 0.0) == ([1,2], [2,3])
    @test interpolate(b, 0.0, [1,1]) == [3,5]
    update!(b, ([2,3], [4,5]))
    @test b == ([2,3], [4,5])

    # tensor field
    c = DVTI(([1 2; 3 4], [2 3; 4 5]))
    @test c[1] == [1 2; 3 4]
    @test c[2] == [2 3; 4 5]
    @test interpolate(c, 0.0) == ([1 2; 3 4], [2 3; 4 5])
    @test interpolate(c, 0.0, [1,1]) == [3 5; 7 9]
    update!(c, ([2 3; 4 5], [5 6; 7 8]))
    @test c == ([2 3; 4 5], [5 6; 7 8])
end

@testset "DCTV field" begin

    # scalar field
    a = DCTV(0.0 => 0.0, 1.0 => 1.0)
    @test isapprox(interpolate(a, -1.0), 0.0)
    @test isapprox(interpolate(a, 0.0), 0.0)
    @test isapprox(interpolate(a, 0.5), 0.5)
    @test isapprox(interpolate(a, 1.0), 1.0)
    update!(a, 1.0 => 2.0)
    @test isapprox(interpolate(a, 0.5), 1.0)
    update!(a, 2.0 => 1.0)
    @test isapprox(interpolate(a, 1.5), 1.5)

    # vector field
    b = DCTV(0.0 => [1.0, 2.0], 1.0 => [2.0, 3.0])
    @test isapprox(interpolate(b, 0.5), [1.5, 2.5])

    # tensor field
    c = DCTV(0.0 => [1.0 2.0; 3.0 4.0], 1.0 => [2.0 3.0; 4.0 5.0])
    @test isapprox(interpolate(c, 0.5), [1.5 2.5; 3.5 4.5])
end
