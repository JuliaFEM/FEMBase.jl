# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase, Test, LinearAlgebra, SparseArrays, Statistics

@testset "FEMBase.jl" begin
    @testset "test_assembly" begin include("test_assembly.jl") end
    @testset "test_elements" begin include("test_elements.jl") end
    @testset "test_fields" begin include("test_fields.jl") end
    @testset "test_integrate" begin include("test_integrate.jl") end
    @testset "test_problems" begin include("test_problems.jl") end
    @testset "test_sparse" begin include("test_sparse.jl") end
    @testset "test_solvers" begin include("test_solvers.jl") end
    @testset "test_analysis" begin include("test_analysis.jl") end
    @testset "test_test" begin include("test_test.jl") end
    @testset "test_types" begin include("test_types.jl") end
end
