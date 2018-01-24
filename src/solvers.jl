# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

type LinearSystem{Tv, Ti<:Integer}
    K :: SparseMatrixCSC{Tv, Ti}
    C1 :: SparseMatrixCSC{Tv, Ti}
    C2 :: SparseMatrixCSC{Tv, Ti}
    D :: SparseMatrixCSC{Tv, Ti}
    f :: SparseVector{Tv, Ti}
    g :: SparseVector{Tv, Ti}
    u :: SparseVector{Tv, Ti}
    la :: SparseVector{Tv, Ti}
end 

abstract type AbstractLinearSystemSolver end

function solve!{S<:AbstractLinearSystemSolver}(::S, ::LinearSystem)
    info("This is a placeholder function for solving linear systems. To solve ",
         "linear systems, you must define a function ",
         "`solve!(solver::$S, system::LinearSystem)`")
end
