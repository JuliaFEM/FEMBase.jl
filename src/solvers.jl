# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

type LinearSystem{Tv, Ti<:Integer}
    M :: SparseMatrixCSC{Tv, Ti}
    K :: SparseMatrixCSC{Tv, Ti}
    Kg :: SparseMatrixCSC{Tv, Ti}
    C1 :: SparseMatrixCSC{Tv, Ti}
    C2 :: SparseMatrixCSC{Tv, Ti}
    D :: SparseMatrixCSC{Tv, Ti}
    f :: SparseVector{Tv, Ti}
    fg :: SparseVector{Tv, Ti}
    g :: SparseVector{Tv, Ti}
    u :: SparseVector{Tv, Ti}
    la :: SparseVector{Tv, Ti}
    dim :: Int
end 

function LinearSystem(dim::Int)
    return LinearSystem(spzeros(dim, dim), spzeros(dim, dim),
                        spzeros(dim, dim), spzeros(dim, dim),
                        spzeros(dim, dim), spzeros(dim, dim),
                        spzeros(dim), spzeros(dim), spzeros(dim),
                        spzeros(dim), spzeros(dim), dim)
end

abstract type AbstractLinearSystemSolver end

function solve!{Solver<:AbstractLinearSystemSolver}(::LinearSystem, ::Solver)
    info("This is a placeholder function for solving linear systems")
    info("To solve linear systems, you must define a function")
    info("      solve!(system::LinearSystem, solver::$Solver)")
end

function can_solve{Solver<:AbstractLinearSystemSolver}(::LinearSystem, ::Solver)
    return (true, "OK")
end

function solve!{S<:AbstractLinearSystemSolver}(ls::LinearSystem, solvers::Vector{S})
    for solver in solvers
        Solver = typeof(solver)
        cansolve, msg = can_solve(ls, solver)
        if !cansolve
            info("Solver $Solver cannot solve linear system: $msg")
            continue
        end
        timeit("solve linear system using solver $Solver") do
            solve!(ls, solver)
        end
        return
    end
    error("Failed to solve linear system.")
end
