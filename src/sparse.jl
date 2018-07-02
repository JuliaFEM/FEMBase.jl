# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

# Sparse utils to make assembly of local and global matrices easier.
# Unoptimized but should do all necessary stuff for at start.

mutable struct SparseMatrixCOO{T<:Real}
    I :: Vector{Int}
    J :: Vector{Int}
    V :: Vector{T}
end

const SparseVectorCOO = SparseMatrixCOO

function SparseMatrixCOO()
    return SparseMatrixCOO{Float64}([], [], [])
end

function SparseVectorCOO(I::Vector, V::Vector)
    return SparseVectorCOO(I, ones(I), V)
end

function convert(::Type{SparseMatrixCOO}, A::SparseMatrixCSC)
    return SparseMatrixCOO(findnz(A)...)
end

function convert(::Type{SparseVectorCOO}, A::SparseVector)
    return SparseVectorCOO(findnz(A)...)
end

function convert(::Type{SparseMatrixCOO}, A::Matrix)
    return SparseMatrixCOO(findnz(A)...)
end

""" Convert from COO format to CSC.

Parameters
----------
tol
    used to drop near zero values less than tol.
"""
function sparse(A::SparseMatrixCOO; tol=1.0e-12)
    B = sparse(A.I, A.J, A.V)
    SparseArrays.droptol!(B, tol)
    return B
end

function sparse(A::SparseMatrixCOO, n::Int, m::Int; tol=1.0e-12)
    B = sparse(A.I, A.J, A.V, n, m)
    SparseArrays.droptol!(B, tol)
    return B
end

function sparse(A::SparseMatrixCOO, n::Int, m::Int, f::Function; tol=1.0e-12)
    B = sparse(A.I, A.J, A.V, n, m, f)
    SparseArrays.droptol!(B, tol)
    return B
end

# For backward compatibility, will be removed
function push!(A::SparseMatrixCOO, I::Int, J::Int, V::Float64)
    push!(A.I, I)
    push!(A.J, J)
    push!(A.V, V)
end

function add!(A::SparseMatrixCOO, I::Int, J::Int, V::Float64)
    push!(A.I, I)
    push!(A.J, J)
    push!(A.V, V)
end

function add!(A::SparseMatrixCOO, I::Int, V::Float64)
    push!(A.I, I)
    push!(A.J, 1)
    push!(A.V, V)
end

function empty!(A::SparseMatrixCOO)
    empty!(A.I)
    empty!(A.J)
    empty!(A.V)
end

function append!(A::SparseMatrixCOO, B::SparseMatrixCOO)
    append!(A.I, B.I)
    append!(A.J, B.J)
    append!(A.V, B.V)
end

function isempty(A::SparseMatrixCOO)
    return isempty(A.I) && isempty(A.J) && isempty(A.V)
end

function full(A::SparseMatrixCOO, args...)
    return full(sparse(A.I, A.J, A.V, args...))
end

""" Add local element matrix to sparse matrix. This basically does:

>>> A[dofs1, dofs2] = A[dofs1, dofs2] + data

Example
-------

>>> S = [3, 4]
>>> M = [6, 7, 8]
>>> data = Float64[5 6 7; 8 9 10]
>>> A = SparseMatrixCOO()
>>> add!(A, S, M, data)
>>> full(A)
4x8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0
 0.0  0.0  0.0  0.0  0.0  5.0  6.0   7.0
 0.0  0.0  0.0  0.0  0.0  8.0  9.0  10.0

"""
function add!(A::SparseMatrixCOO, dofs1::Vector{Int}, dofs2::Vector{Int}, data::Matrix)
    n, m = size(data)
    for j=1:m
        for i=1:n
            add!(A, dofs1[i], dofs2[j], data[i,j])
        end
    end
    return nothing
end

""" Add sparse matrix of CSC to COO. """
function add!(A::SparseMatrixCOO, B::SparseMatrixCSC)
    i, j, v = findnz(B)
    C = SparseMatrixCOO(i, j, v)
    append!(A, C)
end

""" Add new data to COO Sparse vector. """
function add!(A::SparseMatrixCOO, dofs::Vector{Int}, data::Array{Float64}, dim::Int=1)
    if length(dofs) != length(data)
        info("dofs = $dofs")
        info("data = $(vec(data))")
        error("when adding to sparse vector dimension mismatch!")
    end
    append!(A.I, dofs)
    append!(A.J, dim*ones(Int, length(dofs)))
    append!(A.V, vec(data))
end

""" Add SparseVector to SparseVectorCOO. """
function add!(a::SparseVectorCOO, b::SparseVector)
    i, v = findnz(b)
    c = SparseVectorCOO(i, v)
    append!(a, c)
    return
end

""" Combine (I,J,V) values if possible to reduce memory usage. """
function optimize!(A::SparseMatrixCOO)
    i, j, v = findnz(sparse(A))
    A.I = i
    A.J = j
    A.V = v
    return
end

""" Find all nonzero rows from sparse matrix.

Returns
-------

Ordered list of row indices.
"""
function get_nonzero_rows(A::SparseMatrixCSC)
    return sort(unique(rowvals(A)))
end

function get_nonzero_columns(A::SparseMatrixCSC)
    return get_nonzero_rows(transpose(A))
end

function size(A::SparseMatrixCOO)
    isempty(A) && return (0, 0)
    return maximum(A.I), maximum(A.J)
end

function size(A::SparseMatrixCOO, idx::Int)
    return size(A)[idx]
end

""" Resize sparse matrix A to (higher) dimension n x m. """
function resize_sparse(A, n, m)
    return sparse(findnz(A)..., n, m)
end

""" Resize sparse vector b to (higher) dimension n. """
function resize_sparsevec(b, n)
    return sparsevec(findnz(b)..., n)
end

""" Approximative comparison of two matricse A and B. """
function isapprox(A::SparseMatrixCOO, B::SparseMatrixCOO)
    A2 = sparse(A)
    B2 = sparse(B, size(A2)...)
    return isapprox(A2, B2)
end
