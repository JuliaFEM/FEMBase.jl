# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using SparseArrays
import SparseArrays: sparse, sparsevec

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

function convert(::Type{SparseMatrixCOO}, b::Vector)
    I = findall(!iszero, b)
    J = fill(1, size(I))
    V = b[I]
    return SparseMatrixCOO(I, J, V)
end

SparseArrays.sparse(A::SparseMatrixCOO) = sparse(A.I, A.J, A.V)
SparseArrays.sparse(A::SparseMatrixCOO, n::Int, m::Int) = sparse(A.I, A.J, A.V, n, m)
SparseArrays.sparse(A::SparseMatrixCOO, n::Int, m::Int, f::Function) = sparse(A.I, A.J, A.V, n, m, f)
Base.Matrix(A::SparseMatrixCOO) = Matrix(sparse(A))
Base.Matrix(A::SparseMatrixCOO, n::Int, m::Int) = Matrix(sparse(A, n, m))

SparseArrays.sparsevec(b::SparseVectorCOO) = sparsevec(b.I, b.V)
SparseArrays.sparsevec(b::SparseVectorCOO, n::Int) = sparsevec(b.I, b.V, n)
Base.Vector(b::SparseVectorCOO) = Vector(sparsevec(b))
Base.Vector(b::SparseVectorCOO, n::Int) = Vector(sparsevec(b, n))

function add!(A::SparseMatrixCOO, I::Int, J::Int, V::Float64)
    push!(A.I, I)
    push!(A.J, J)
    push!(A.V, V)
    return nothing
end

function add!(A::SparseMatrixCOO, I::Int, V::Float64)
    push!(A.I, I)
    push!(A.J, 1)
    push!(A.V, V)
    return nothing
end

function empty!(A::SparseMatrixCOO)
    empty!(A.I)
    empty!(A.J)
    empty!(A.V)
    return nothing
end

function append!(A::SparseMatrixCOO, B::SparseMatrixCOO)
    append!(A.I, B.I)
    append!(A.J, B.J)
    append!(A.V, B.V)
    return nothing
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

""" Approximative comparison of two matrices A and B. """
function isapprox(A::SparseMatrixCOO, B::SparseMatrixCOO)
    A2 = sparse(A)
    B2 = sparse(B, size(A2)...)
    return isapprox(A2, B2)
end

isapprox(A::SparseMatrixCOO, B) = isapprox(Matrix(A), B)
isapprox(A, B::SparseMatrixCOO) = isapprox(A, Matrix(B))

# TODO: Use FEMSParse.jl
struct AssemblerSparsityPattern{Tv,Ti}
    K::SparseMatrixCSC{Tv,Ti}
    permutation::Vector{Int}
    sorteddofs::Vector{Int}
end

"""
    start_assemble(K::SparseMatrixCSC)
Create an `AssemblerSparsityPattern` from the sparsity pattern in `K`.
Also, fill the values in `K` with 0.0.
"""
function start_assemble(K::SparseMatrixCSC)
    fill!(K.nzval, 0.0)
    AssemblerSparsityPattern(K, Int[], Int[])
end


"""
    assemble!(A::AssemblerSparsityPattern, dofs2, Ke)
Assemble a local dense element matrix `Ke` into the sparse matrix
wrapped by the assembler `A`.
location given by lists of indices `dofs1` and `dofs2`.
# Example
```julia
using SparseArrays
sparsity_pattern = sparse([1. 0 1; 1 0 1; 1 1 1])
assembler = FEMSparse.start_assemble(sparsity_pattern)
dofs = [1, 3]
Ke = [1.0 2.0; 3.0 4.0]
FEMSparse.FEMSparse.assemble_local_matrix!(assembler, dofs, Ke)
Matrix(sparsity_pattern)
# output
julia> Matrix(sparsity_pattern)
3×3 Array{Float64,2}:
 1.0  0.0  2.0
 0.0  0.0  0.0
 3.0  0.0  4.0
```
"""
function assemble_local_matrix!(A::AssemblerSparsityPattern, dofs::AbstractVector{Int}, Ke::AbstractMatrix)
    permutation = A.permutation
    sorteddofs = A.sorteddofs
    K = A.K

    @boundscheck checkbounds(K, dofs, dofs)
    resize!(permutation, length(dofs))
    resize!(sorteddofs, length(dofs))
    copyto!(sorteddofs, dofs)
    sortperm2!(sorteddofs, permutation)

    current_col = 1
    @inbounds for Kcol in sorteddofs
        maxlookups = length(dofs)
        current_idx = 1
        for r in nzrange(K, Kcol)
            Kerow = permutation[current_idx]
            if K.rowval[r] == dofs[Kerow]
                Kecol = permutation[current_col]
                K.nzval[r] += Ke[Kerow, Kecol]
                current_idx += 1
            end
            current_idx > maxlookups && break
        end
        if current_idx <= maxlookups
            error("some row indices were not found")
        end
        current_col += 1
    end
end


##################
# Sort utilities #
##################
# We bundle a few sorting utilities here because the ones
# in have some unacceptable overhead.

# Sorts B and stores the permutation in `ii`
function sortperm2!(B, ii)
   @inbounds for i = 1:length(B)
      ii[i] = i
   end
   quicksort!(B, ii)
   return
end

function quicksort!(A, order, i=1,j=length(A))
    @inbounds if j > i
        if  j - i <= 12
           # Insertion sort for small groups is faster than Quicksort
           insertionsort!(A, order, i, j)
           return A
        end

        pivot = A[div(i+j,2)]
        left, right = i, j
        while left <= right
            while A[left] < pivot
                left += 1
            end
            while A[right] > pivot
                right -= 1
            end
            if left <= right
                A[left], A[right] = A[right], A[left]
                order[left], order[right] = order[right], order[left]

                left += 1
                right -= 1
            end
        end  # left <= right

        quicksort!(A,order, i,   right)
        quicksort!(A,order, left,j)
    end  # j > i

    return A
end

function insertionsort!(A, order, ii=1, jj=length(A))
    @inbounds for i = ii+1 : jj
        j = i - 1
        temp  = A[i]
        itemp = order[i]

        while true
            if j == ii-1
                break
            end
            if A[j] <= temp
                break
            end
            A[j+1] = A[j]
            order[j+1] = order[j]
            j -= 1
        end

        A[j+1] = temp
        order[j+1] = itemp
    end  # i
    return
end
