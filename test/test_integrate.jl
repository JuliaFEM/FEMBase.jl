# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Test

@testset "test integrate" begin
    element = Element(Seg2, [1, 2])
    ips1 = get_integration_points(element)
    ips2 = get_integration_points(element, 1)
    @test length(ips1) == 2
    @test length(ips2) == 3
    element = Element(Poi1, [1])
    ips3 = get_integration_points(element)
    @test length(ips3) == 1
end

