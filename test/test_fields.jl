# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Test

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
    update!(a, (2,3))
    @test a == (2,3)
    @test (2,3) == a

    # vector field
    b = DVTI(([1,2], [2,3]))
    @test b[1] == [1,2]
    @test b[2] == [2,3]
    @test interpolate(b, 0.0) == ([1,2], [2,3])
    update!(b, ([2,3], [4,5]))
    @test b == ([2,3], [4,5])

    # tensor field
    c = DVTI(([1 2; 3 4], [2 3; 4 5]))
    @test c[1] == [1 2; 3 4]
    @test c[2] == [2 3; 4 5]
    @test interpolate(c, 0.0) == ([1 2; 3 4], [2 3; 4 5])
    update!(c, ([2 3; 4 5], [5 6; 7 8]))
    @test c == ([2 3; 4 5], [5 6; 7 8])

    d = DVTI(2, 3)
    @test a == d
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

@testset "DVTV field" begin
    # scalar field
    a = DVTV(0.0 => (0.0, 1.0), 1.0 => (1.0, 0.0))
    update!(a, 2.0 => (2.0, 0.0))
    r = interpolate(a, 0.5)
    @test isapprox(r[1], 0.5)
    @test isapprox(r[2], 0.5)
    update!(a, 2.0 => (4.0, 0.0))
end

@testset "CVTV field" begin
    f = CVTV((xi,t) -> xi[1]*xi[2]*t)
    @test isapprox(f([1.0,2.0],3.0), 6.0)
end

@testset "Dictionary fields" begin
    X = Dict(1=>[0.0,0.0], 1000=>[1.0,0.0], 100000=>[1.0,1.0])
    G = DVTId(X)
    @test isapprox(G[1], X[1])
    @test isapprox(G[1000], X[1000])
    @test isapprox(G[100000], X[100000])
    Y = Dict(1=>[2.0,2.0], 1000=>[3.0,2.0], 100000=>[3.0,3.0])
    F = DVTVd(0.0 => X, 1.0 => Y)
    @test isapprox(interpolate(F, 0.5)[100000], [2.0,2.0])
end

@testset "update dictionary field" begin
    f1 = Dict(1=>1.0, 2=>2.0, 3=>3.0)
    f2 = Dict(1=>2.0, 2=>3.0, 3=>4.0)
    fld = DVTVd(0.0 => f1)
    update!(fld, 1.0 => f2)
    @test isapprox(interpolate(fld, 0.5)[1], 1.5)
    update!(fld, 1.0 => f1)
    @test isapprox(interpolate(fld, 0.5)[1], 1.0)
end

@testset "use of common constructor field" begin
    @test isa(field(1.0), DCTI)
    @test isa(field(1.0 => 1.0), DCTV)
    @test isa(field((1.0,2.0)), DVTI)
    @test isa(field(1, 2), DVTI)
    @test isa(field(1.0 => (1.0,2.0)), DVTV)
    @test isa(field((xi,t) -> xi[1]*t), CVTV)
    @test isa(field(1 => [1.0, 2.0], 10 => [2.0, 3.0]), DVTId)
    @test isa(field(0.0 => (1=>1.0,10=>2.0), 1.0 => (1=>2.0,10=>3.0)), DVTVd)
    X = Dict(1 => [0.0,0.0], 2 => [1.0,0.0])
    X1 = field(X)
    X2 = field(0.0 => X)
    @test isa(X1, DVTId)
    @test isa(X2, DVTVd)
end

@testset "general interpolation" begin
    a = [1, 2, 3]
    b = (2, 3, 4)
    @test interpolate(a, b) == 2+6+12
    a = (1, 2)
    b = (2, 3, 4)
    @test interpolate(a, b) == 2+6
    @test_throws AssertionError interpolate(b, a)
end
