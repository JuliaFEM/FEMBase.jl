# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

"""
    FEMBase.FEMTest

FEMBase.FEMTest contains utilities to test JuliaFEM extensions, including:

- all test macros from Base.Test (@test, @testset, ...)
- test problems `Poisson` and `Dirichlet`.
- test linear system solver `LUSolver`.
- test analysis `Static`.
- read ABAQUS matrices

"""
module FEMTest

# For convenience, one can also get all @test macros and so on just by typing
# `using FEMBase.FEMTest`, they are the very same than in `Base.Test`
using Test, LinearAlgebra, SparseArrays

export @test, @test_throws, @test_broken, @test_skip,
       @test_warn, @test_nowarn, @testset, @inferred

using ..FEMBase
import ..FEMBasis

struct Poisson <: FieldProblem end

function FEMBase.get_unknown_field_name(::Problem{Poisson})
    return "u"
end

function FEMBase.assemble_elements!(problem::Problem{Poisson},
                                    assembly::Assembly,
                                    elements::Vector{Element{E}},
                                    time::Float64) where E

    bi = FEMBasis.BasisInfo(E)
    ndofs = length(bi)
    Ke = zeros(ndofs, ndofs)
    fe = zeros(ndofs)

    for element in elements
        fill!(Ke, 0.0)
        fill!(fe, 0.0)
        for ip in get_integration_points(element)
            J, detJ, N, dN = element_info!(bi, element, ip, time)
            k = 1.0
            if haskey(element, "coefficient")
                k = element("coefficient", ip, time)
            end
            Ke += ip.weight * k*dN'*dN * detJ
            if haskey(element, "source")
                f = element("source", ip, time)
                fe += ip.weight * N'*f * detJ
            end
        end
        gdofs = get_gdofs(problem, element)
        add!(assembly.K, gdofs, gdofs, Ke)
        add!(assembly.f, gdofs, fe)
    end
end

function FEMBase.assemble_elements!(problem::Problem{Poisson},
                                    assembly::Assembly,
                                    elements::Vector{Element{E}},
                                    time::Float64) where E<:Union{FEMBasis.Seg2,FEMBasis.Seg3}

    bi = FEMBasis.BasisInfo(E)
    ndofs = length(bi)
    Ce = zeros(ndofs, ndofs)
    fe = zeros(ndofs)
    ge = zeros(ndofs)

    for element in elements
        fill!(Ce, 0.0)
        fill!(fe, 0.0)
        fill!(ge, 0.0)
        for ip in get_integration_points(element)
            J, detJ, N, dN = element_info!(bi, element, ip, time)
            if haskey(element, "flux")
                f = element("flux", ip, time)
                fe += ip.weight * N'*f * detJ
            end
            if haskey(element, "fixed u")
                f = element("fixed u", ip, time)
                Ce += ip.weight * N'*N * detJ
                ge += ip.weight * N'*f * detJ
            end
        end
        gdofs = get_gdofs(problem, element)
        add!(assembly.C1, gdofs, gdofs, Ce')
        add!(assembly.C2, gdofs, gdofs, Ce)
        add!(assembly.f, gdofs, fe)
        add!(assembly.g, gdofs, ge)
    end
end

struct Dirichlet <: BoundaryProblem end

function FEMBase.assemble_elements!(problem::Problem{Dirichlet},
                                    assembly::Assembly,
                                    elements::Vector{Element{E}},
                                    time::Float64) where E

    name = get_parent_field_name(problem)
    dim = get_unknown_field_dimension(problem)

    data = Dict{Int,Float64}()
    for element in elements
        for i=1:dim
            haskey(element, "$name $i") || continue
            gdofs = get_gdofs(problem, element)
            ldofs = gdofs[i:dim:end]
            xis = FEMBasis.get_reference_element_coordinates(E)
            for (ldof, xi) in zip(ldofs, xis)
                data[ldof] = interpolate(element, "$name $i", xi, time)
            end
        end
    end

    for (dof, val) in data
        add!(assembly.C1, dof, dof, 1.0)
        add!(assembly.C2, dof, dof, 1.0)
        add!(assembly.g, dof, val)
    end

end

struct LUSolver <: AbstractLinearSystemSolver
end

"""
    solve!(solver::LUSolver, ls::LinearSystem)

Solve linear system using LU factorization. If final system has zero rows,
add 1 to diagonal to make matrix non-singular.
"""
function FEMBase.solve!(ls::LinearSystem, ::LUSolver)
    A = [ls.K ls.C1; ls.C2 ls.D]
    b = [ls.f; ls.g]

    # add 1.0 to diagonal for any zero rows in system
    p = ones(2*ls.dim)
    p[unique(SparseArrays.rowvals(A))] .= 0.0
    A += sparse(Diagonal(p))

    # solve A*x = b using LU factorization and update solution vectors
    x = lu(A) \ Vector(b)
    ls.u = x[1:ls.dim]
    ls.la = x[ls.dim+1:end]

    return nothing
end

mutable struct Static <: AbstractAnalysis
    time :: Float64
end

function Static()
    return Static(0.0)
end

"""
    get_linear_system(problems, time)

Assemble `problems` and construct LinearSystem for test problem.
"""
function get_linear_system(problems, time)
    # assemble matrices for all problems
    N = 0 # size of resulting matrix
    for problem in problems
        assemble!(problem, time)
        N = max(N, size(problem.assembly.K, 2))
    end
    # create new LinearSystem and add assemblies to that
    ls = LinearSystem(N)
    for problem in problems
        ls.M += sparse(problem.assembly.M, N, N)
        ls.K += sparse(problem.assembly.K, N, N)
        ls.Kg += sparse(problem.assembly.K, N, N)
        ls.f += sparse(problem.assembly.f, N, 1)
        ls.fg += sparse(problem.assembly.f, N, 1)
        ls.C1 += sparse(problem.assembly.C1, N, N)
        ls.C2 += sparse(problem.assembly.C2, N, N)
        ls.D += sparse(problem.assembly.D, N, N)
        ls.g += sparse(problem.assembly.g, N, 1)
    end
    return ls
end

"""
    run!(analysis::Analysis{Static})

Run static analysis for test problem. This is for testing purposes only and
should not be used in real simulations.
"""
function FEMBase.run!(analysis::Analysis{Static})
    @info("This is from FEMBase.Test.Analysis{Static}, which is for testing " *
          "purposes only.")
    ls = get_linear_system(get_problems(analysis), analysis.properties.time)
    solve!(ls, LUSolver())
    normu = norm(ls.u)
    normla = norm(ls.la)
    @info("Solution norms are |u| = $normu, |la| = $normla")
    return ls, normu, normla
end

function read_mtx(data::IO; dim=0)
    if dim == 0
        for ln in eachline(data)
            i,idof,j,jdof,value = map(Meta.parse, split(ln, ','))
            dim = max(dim, idof, jdof)
        end
    end
    seekstart(data)
    I = Int[]
    J = Int[]
    V = Float64[]
    for ln in eachline(data)
        i,idof,j,jdof,value = map(Meta.parse, split(ln, ','))
        if i < 1 || j < 1
            continue
        end
        push!(I, (i-1)*dim+idof)
        push!(J, (j-1)*dim+jdof)
        push!(V, value)
    end
    A = Matrix(sparse(I, J, V))
    A += transpose(tril(A,-1))
    return A
end

"""
    read_mtx_from_file(filename::String; dim=0)

Read matrix data export from ABAQUS by using command:

*STEP
*MATRIX GENERATE, STIFFNESS, MASS, LOAD
*MATRIX OUTPUT, STIFFNESS, MASS, LOAD
*END STEP

If number of dofs / node (dim) is not given, it will be automatically
detemined from file.

"""
function read_mtx_from_file(filename::String; dim=0)
    open(filename) do fid
        return read_mtx(fid; dim=dim)
    end
end

function read_mtx_from_string(data::String; dim=0)
    return read_mtx(IOBuffer(data); dim=dim)
end

"""
    test_resource(resource::String) -> String

`@test_resource(resource)` expands to a string containing full path to the some
`resource` what is needed by a test. 

# Example

Say, test file name is `test_run_model.jl` inside package called `Models`.
One needs a test mesh file `mesh.inp`. Then, function

```julia
mesh_file = @test_resource("mesh.inp")
```

expands to 

    `~/.julia/v0.6/Models/test/test_run_model/mesh.inp`

Macro must be used inside the actual test file.
"""
macro test_resource(resource::String)
    #__source__.file === nothing && return nothing
    #filename = String(__source__.file)
    filename = Base.source_path()
    abs_dirname = abspath(dirname(filename))
    test_dir_name = first(splitext(basename(filename)))
    resource_ = joinpath(abs_dirname, test_dir_name, resource)
    if !isfile(resource_)
        @warn("Cannot find resource $resource_")
    end
    return resource_
end

export Poisson, Dirichlet, Static, LUSolver, get_linear_system, abaqus_read_mtx
export read_mtx_from_file, read_mtx_from_string, @test_resource

end

global const Test = FEMTest
