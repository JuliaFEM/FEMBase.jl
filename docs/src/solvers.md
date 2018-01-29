# Analyses and solvers

After a global assembly for each problem is ready, they must somehow put together,
usually construct a linear system ``\boldsymbol{A}\boldsymbol{x} = \boldsymbol{b}`` ,
solve system of equations and potentially update solution back to problems / elements.

For this task we have [`Analysis`](@ref). Just like problems are "containers" for a
set of some elements, analyses are containers for a set of problems, performing some
spesific task. This task usually (but not always) involves solving big systems of
equations.

For example, typical nonlinear quasistatic analysis is

1. initialize problems (if needed)
2. assemble problems
3. combine assemblies to construct a linear system
   ``\boldsymbol{A}\boldsymbol{x} = \boldsymbol{b}``
4. solve linear system, as effectively as possible
5. update solution back to problems or elements
6. postprocess fields
7. write requested results to file
8. check convergence, if not converged, go back to 2.

## Implementing solver for LinearSystem

When models gets big enough, step 4 is dominating in a solution process. For
that reason we have abstract type [`AbstractLinearSystemSolver`](@ref) which
can be subclassed to construct own solution strategy. This strategy can be, for
example to use Julia's build-in solvers, MUMPS, iterative solvers and so on. And
to make this as standard as possible, we have [`LinearSystem`](@ref) containing
all the relevant matrices, i.e.
```math
\begin{bmatrix}
\boldsymbol{K} & \boldsymbol{C}_{1}\\
\boldsymbol{C}_{2} & \boldsymbol{D}
\end{bmatrix}
\begin{bmatrix}
\boldsymbol{u}\\
\boldsymbol{\lambda}
\end{bmatrix} =
\begin{bmatrix}
\boldsymbol{f}\\
\boldsymbol{g}
\end{bmatrix}.
```

```@setup FEMBase
using FEMBase
using FEMBase: LinearSystem, AbstractLinearSystemSolver,
               AbstractAnalysis, Analysis, Element, Problem,
               Quad4, get_problems
import FEMBase: assemble_elements!, run!, can_solve, solve!

type Heat <: FieldProblem
end

function get_unknown_field_name(::Problem{Heat})
    return "temperature"
end

function assemble_elements!{B}(problem::Problem{Heat}, assembly::Assembly,
                               elements::Vector{Element{B}}, time::Float64)

    println("Assembling volume elements of kind $B")
    bi = BasisInfo(B)
    ndofs = length(B)
    Ke = zeros(ndofs, ndofs)
    fe = zeros(ndofs)

    for element in elements
        fill!(Ke, 0.0)
        fill!(fe, 0.0)
        for ip in get_integration_points(element)
            J, detJ, N, dN = element_info!(bi, element, ip, time)
            k = element("thermal conductivity", ip, time)
            Ke += ip.weight * k*dN'*dN * detJ
            if haskey(element, "heat source")
                f = element("heat source", ip, time)
                fe += ip.weight * N'*f * detJ
            end
        end
        gdofs = get_gdofs(problem, element)
        add!(assembly.K, gdofs, gdofs, Ke)
        add!(assembly.f, gdofs, fe)
    end

end

function assemble_elements!{B<:Union{Seg2,Seg3}}(problem::Problem{Heat},
                                                 assembly::Assembly,
                                                 elements::Vector{Element{B}},
                                                 time::Float64)

    println("Assembling boundary elements of kind $B")
    bi = BasisInfo(B)
    ndofs = length(B)
    fe = zeros(ndofs)

    for element in elements
        haskey(element, "heat flux") || continue
        fill!(fe, 0.0)
        for ip in get_integration_points(element)
            J, detJ, N, dN = element_info!(bi, element, ip, time)
            g = element("heat flux", ip, time)
            fe += ip.weight * N'*g * detJ
        end
        gdofs = get_gdofs(problem, element)
        add!(assembly.f, gdofs, fe)
    end

end

type DirichletBC <: BoundaryProblem
end

function assemble_elements!{E}(problem::Problem{DirichletBC},
                               assembly::Assembly,
                               elements::Vector{Element{E}},
                               time::Float64)

    name = get_parent_field_name(problem)
    dim = get_unknown_field_dimension(problem)

    println("Assembling Dirichlet boundary condition")
    println("Field name = $name, dofs/node = $dim")

    data = Dict{Int64,Float64}()
    for element in elements
        for i=1:dim
            haskey(element, "$name $dim") || continue
            gdofs = get_gdofs(problem, element)
            ldofs = gdofs[i:dim:end]
            xis = get_reference_element_coordinates(E)
            for (ldof, xi) in zip(ldofs, xis)
                data[ldof] = interpolate(element, "$name $dim", xi, time)
            end
        end
    end

    for (dof, val) in data
        add!(assembly.C1, dof, dof, 1.0)
        add!(assembly.C2, dof, dof, 1.0)
        add!(assembly.g, dof, val)
    end

end

```

First setup linear system, we take matrices introduced in the context of [`Problems`]:

```@example FEMBase
ls = LinearSystem(5)
ls.K = sparse([  4.0  -1.0  -2.0  -1.0   0.0
                -1.0   7.0  -4.0  -2.0   0.0
                -2.0  -4.0  10.0  -1.0  -3.0
                -1.0  -2.0  -1.0   4.0   0.0
                 0.0   0.0  -3.0   0.0   3.0])
ls.C1 = sparse([ 1.0  0.0  0.0  0.0  0.0
	         0.0  0.0  0.0  0.0  0.0
	         0.0  0.0  0.0  0.0  0.0
	         0.0  0.0  0.0  1.0  0.0
	         0.0  0.0  0.0  0.0  0.0])
ls.C2 = ls.C1
ls.f = sparsevec([33.0, 33.0, 165.0, 33.0, 132.0])
nothing # hide
```

For example, a simple solver for small problems would be to use `lufact` from
UMFPACK:

```@example FEMBase

type LUSolver <: AbstractLinearSystemSolver
    # may contain some solver-spesific settings
end

"""
    solve!(ls::LinearSystem, solver::LUSolver)

Solve linear system using LU factorization. If final system has zero rows,
add 1 to diagonal to make matrix non-singular.
"""
function solve!(ls::LinearSystem, solver::LUSolver)

    println("Solving system using LUSolver")

    A = [ls.K ls.C1; ls.C2 ls.D]
    b = [ls.f; ls.g]

    # add 1.0 to diagonal for any zero rows in system
    p = ones(2*ls.dim)
    p[unique(rowvals(A))] = 0.0
    A += spdiagm(p)

    # solve A*x = b using LU factorization and update solution vectors
    x = lufact(A) \ full(b)
    ls.u = x[1:ls.dim]
    ls.la = x[ls.dim+1:end]

    return nothing
end

nothing # hide
```

To test it:
```@repl FEMBase
solver = LUSolver()
solve!(ls, solver)
full(ls.u)
full(ls.la)
```

We try to be clever with solvers, because solution of linear system is taking
the most of the time in analysis for big models. It's possible to build a list
of solvers which are then run in order, so if first one cannot solve the system,
then we can try the next one and so on. For example, direct solvers are very
robust but they need quite lot of memory for factorization. If computer / cluster
has enough of memory, we can prefer direct solver over iterative ones. If
stiffness matrix is symmetric, and problem is positive definite, Cholesky
decomposition is preferred over LU decomposition and so on. Before starting the
actual solution, one can study matrix structure using [`can_solve`](@ref):

```@example FEMBase
type FakeSolver <: AbstractLinearSystemSolver
end

function can_solve(ls::LinearSystem, solver::FakeSolver)
    if issymmetric(ls.K)
        return (true, "OK")
    else
        return (false, "stiffness matrix not symmetric")
    end
end

nothing # hide
```

```@repl FEMBase
solvers = [FakeSolver(), LUSolver()]
ls.K[1,2] += 1
solve!(ls, solvers)
```

## Implementing new analyses

Like said, [`Analysis`](@ref) is a type containing all problems and is performing
some sort of operation with them. Typically this can be a solution of static or
dynamic equilibrium, eigenvalue analysis, linear perturbation analysis or doing
some postprocessing for already done analysis, like fatigue analysis, calculating
critical stresses and so on.

In industrial computing, a bigger workflow typically contains several different
analyses. Compared to some commerfial FEM softwares, e.g. ABAQUS, `Analysis` can
be consided as STEP, but in a more general sense. Because Julia is a real
programming language, one can then create, for example, an optimization loop,
combining several different analyses to optimize model against some design goal.

There's not so much rules about how to define new analysis. It must be a
subtype of `AbstractAnalysis` and implement `run!`-function. As a minimal
working example, implementation for static analysis would then be:

```@example FEMBase
type Static <: AbstractAnalysis
    time :: Float64
end

function Static()
    return Static(0.0)
end

function run!(analysis::Analysis{Static})

    time = analysis.properties.time
    problems = get_problems(analysis)

    # assemble matrices for all problems
    N = 0 # size of resulting matrix
    for problem in problems
        assemble!(problem, time)
        N = max(N, size(problem.assembly.K, 2))
    end

    # create new LinearSystem and add assemblies to that
    ls = LinearSystem(N)

    for problem in problems
        ls.K += sparse(problem.assembly.K, N, N)
        ls.f += sparse(problem.assembly.f, N, 1)
        ls.C1 += sparse(problem.assembly.C1, N, N)
        ls.C2 += sparse(problem.assembly.C2, N, N)
        ls.D += sparse(problem.assembly.D, N, N)
        ls.g += sparse(problem.assembly.g, N, 1)
    end

    # solve linear system using e.g. LUSolver
    solve!(ls, LUSolver())

    # update solution back to problems
    u = Dict(i=>j for (i,j) in enumerate(ls.u))
    for problem in problems
        is_field_problem(problem) || continue
        field_name = get_unknown_field_name(problem)
        elements = get_elements(problem)
        update!(elements, field_name, time => u)
    end
end

nothing # hide
```

To put everything together, typical solution process of time-intependent
variational problem looks the following:

```@example FEMBase
X = Dict(1 => [0.0, 0.0],
         2 => [1.0, 0.0],
         3 => [1.0, 1.0],
         4 => [0.0, 1.0],
         5 => [2.0, 1.0])

el1 = Element(Quad4, [1, 2, 3, 4])
el2 = Element(Tri3, [3, 2, 5])
el3 = Element(Seg2, [3, 5])
elements = [el1, el2, el3]
update!(elements, "geometry", X)
update!(elements, "thermal conductivity", 6.0)
update!(el3, "heat flux", 264.0)
update!(el1, "heat source", 132.0)
problem = Problem(Heat, "test problem", 1)
add_elements!(problem, elements)

bel1 = Element(Seg2, [1, 4])
update!(bel1, "geometry", X)
update!(bel1, "temperature 1", 0.0)
bc = Problem(DirichletBC, "fixed", 1, "temperature")
add_elements!(bc, [bel1])

analysis = Analysis(Static, "static analysis of test problem")
add_problems!(analysis, [problem, bc])
run!(analysis)
```

After analysis is ready, temperature is updated to elements:

```@repl FEMBase
el1("temperature", 0.0)
el2("temperature", 0.0)
el3("temperature", 0.0)
```
