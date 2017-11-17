# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using Base.Test
using TimerOutputs
const to = TimerOutput()

test_files = String[]
push!(test_files, "test_add_elements.jl")
push!(test_files, "test_common_failures.jl")
push!(test_files, "test_elements.jl")
push!(test_files, "test_elements_2.jl")
push!(test_files, "test_fields.jl")
push!(test_files, "test_fields_time_interpolation.jl")
push!(test_files, "test_integration_points.jl")
push!(test_files, "test_sparse.jl")

@testset "FEMBase.jl" begin
    for fn in test_files
        timeit(to, fn) do
            include(fn)
        end
    end
end
println()
println("Test statistics:")
println(to)
