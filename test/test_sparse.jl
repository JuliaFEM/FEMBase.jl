# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using FEMBase: SparseMatrixCOO, SparseVectorCOO
using FEMBase: resize_sparse, resize_sparsevec
using FEMBase: get_nonzero_rows, get_nonzero_columns
using Test
using SparseArrays

@testset "Add to SparseMatrixCOO" begin
    A = SparseMatrixCOO()
    A2 = reshape(collect(1:9), 3, 3)
    add!(A, sparse(A2))
    @test isapprox(A, A2)
end

@testset "Add to SparseVectorCOO" begin
    b = SparseVectorCOO()
    b2 = collect(1:3)
    add!(b, sparse(b2))
    @test isapprox(b, b2)
end

@testset "Failure to add data to sparse vector due dimension mismatch" begin
    b = SparseVectorCOO()
    @test_throws ErrorException add!(b, [1, 2], [1.0, 2.0, 3.0])
end

@testset "resize of sparse matrix and sparse vector" begin
    A = sparse(rand(3, 3))
    B = resize_sparse(A, 4, 4)
    @test size(B) == (4, 4)

    a = sparse(rand(3))
    b = resize_sparsevec(a, 4)
    @test size(b) == (4, )
end

@testset "type conversions" begin
    A = reshape(collect(1:9), 3, 3)
    A2 = convert(SparseMatrixCOO, A)
    A3 = convert(SparseMatrixCOO, sparse(A))
    @test isapprox(A, A2)
    @test isapprox(A, A3)
    b = sparse(collect(1:3))
    b2 = convert(SparseVectorCOO, b)
    @test isapprox(b, b2)
end

@testset "construct sparse matrix and empty it" begin
    A = SparseMatrixCOO()
    add!(A, 1, 1, 1.0)
    A_expected = sparse([1], [1], [1.0], 2, 2)
    A2 = sparse(A, 2, 2)
    @test isapprox(A2, A_expected)
    empty!(A)
    @test length(A.I) == length(A.J) == length(A.V) == 0
end

@testset "construct sparse matrix with special function" begin
    A = SparseMatrixCOO()
    add!(A, 1, 1, 1.0)
    add!(A, 1, 1, 2.0)
    B = sparse(A, 2, 2, +)
    @test isapprox(B[1,1], 3.0)
end

@testset "normal addition to sparse vector" begin
    A = SparseMatrixCOO()
    add!(A, [1, 2], [1.0 2.0])
end

@testset "get nonzero rows, columns, size" begin
    A = SparseMatrixCOO()
    add!(A, [1, 3], [1, 4], [1.0 2.0; 3.0 4.0])
    B = sparse(A)
    @test get_nonzero_rows(B) == [1, 3]
    @test get_nonzero_columns(B) == [1, 4]
    @test size(A) == (3, 4)
    @test size(A, 1) == 3
end

@testset "isapprox" begin
    A = SparseMatrixCOO()
    B = SparseMatrixCOO()
    add!(A, [1, 3], [1, 4], [1.0 2.0; 3.0 4.0])
    add!(B, [1, 3], [1, 4], [1.0 2.0; 3.0 4.0])
    @test isapprox(A, B)
end
