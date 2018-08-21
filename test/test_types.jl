# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using FEMBase: IP
using Test

@testset "test integration point" begin
    a = sqrt(1.0/3.0)
    ip = IP(1, 1.0, (-a, a))
    strain = [1.0 2.0; 3.0 4.0]
    update!(ip, "strain", 0.0 => strain)
    @test isapprox(ip("strain", 0.0), strain)
    @test isapprox(ip[1], -a)
    for p in ip
        @test isapprox(abs.(p), a)
    end
    update!(ip, "strain", 1.0 => 2*strain)
    @test isapprox(ip("strain", 1.0), 2*strain)
    ip2 = IP(1, 1.0, [-a, a])
    @test isapprox(ip2[1], -a)
end

