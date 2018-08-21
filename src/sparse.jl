# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using SparseArrays
import SparseArrays: sparse

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
    return SparseVectorCOO(I, fill!(similar(I), 1), V)
end

function convert(::Type{SparseMatrixCOO}, A::SparseMatrixCSC)
    return SparseMatrixCOO(findnz(A)...)
end

function convert(::Type{SparseVectorCOO}, A::SparseVector)
    return SparseVectorCOO(findnz(A)...)
end

function convert(::Type{SparseMatrixCOO}, A::Matrix)
    idx = findall(!iszero, A)
    I = getindex.(idx, 1)
    J = getindex.(idx, 2)
    V = [A[i] for i in idx]
    return SparseMatrixCOO(I, J, V)
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

"""
    add!(K, dofs1, dofs2, ke)

Add local element matrix `ke` to sparse matrix `K` for indices defined by `dofs1`
and `dofs2`. This basically does `A[dofs1, dofs2] = A[dofs1, dofs2] + data`.

# Examples

```julia
S = [3, 4]
M = [6, 7, 8]
ke = [5 6 7; 8 9 10]
K = SparseMatrixCOO()
add!(K, S, M, ke)
Matrix(A)

# output

4x8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0
 0.0  0.0  0.0  0.0  0.0  5.0  6.0   7.0
 0.0  0.0  0.0  0.0  0.0  8.0  9.0  10.0
```
"""
function add!(A::SparseMatrixCOO, dofs1::Vector{Int}, dofs2::Vector{Int}, data)
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
        @error("Dimension mismatch when adding data to sparse vector!", dofs, data)
        error("Simulation stopped.")
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

"""
    get_nonzero_rows(A)

Returns indices of all nonzero rows from a sparse matrix `A`.
"""
function get_nonzero_rows(A)
    return sort(unique(A.rowval))
end

"""
    get_nonzero_columns(A)

Returns indices of all nonzero columns from a sparse matrix `A`.
"""
function get_nonzero_columns(A)
    return get_nonzero_rows(copy(transpose(A)))
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
    idx = findall(!iszero, A)
    I = getindex.(idx, 1)
    J = getindex.(idx, 2)
    V = [A[i] for i in idx]
    return sparse(I, J, V, n, m)
end

""" Resize sparse vector b to (higher) dimension n. """
function resize_sparsevec(b, n)
    return sparsevec(b.nzind, b.nzval, n)
end

""" Approximative comparison of two matricse A and B. """
function isapprox(A::SparseMatrixCOO, B::SparseMatrixCOO)
    A2 = sparse(A)
    B2 = sparse(B, size(A2)...)
    return isapprox(A2, B2)
end

function isapprox(A::SparseMatrixCOO, B)
    return isapprox(Matrix(sparse(A)), B)
end

function isapprox(A, B::SparseMatrixCOO)
    return isapprox(B, A)
end
