# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase, Test

@testset "get integration points for Element{Seg2}" begin
    element = Element(Seg2, (1, 2))
    ips = get_integration_points(element)
    @test length(ips) == 2
end
