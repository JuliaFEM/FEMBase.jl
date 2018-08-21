# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Test

# Time interpolation of fields

element = Element(Seg2, [1, 2])
update!(element, "force", 1.0 => 15.0)
update!(element, "force", 2.0 => 30.0)
@test isapprox(element("force", 1.0), 15.0)
@test isapprox(element("force", 1.2), 18.0)
@test isapprox(element("force", 1.6), 24.0)
@test isapprox(element("force", 2.0), 30.0)

update!(element, "temperature", 0.0 => (1.0, 3.0))
update!(element, "temperature", 1.0 => (2.0, 14.0))
@test isapprox(element("temperature", (0.0,), 0.0), 2.0)
@test isapprox(element("temperature", (0.0,), 0.3), 3.8)
@test isapprox(element("temperature", (0.0,), 1.0), 8.0)
