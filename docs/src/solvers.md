# Solvers

After global assembly for each problem is ready, they must somehow put together,
construct linear system ``\boldsymbol{A}\boldsymbol{x} = \boldsymbol{b}`` , solve
it and update solution back to problems / elements.

Just like problems are "containers" for a set of elements, solvers are containers
for a set of problems. 

Typical solver structure for nonlinear quasistatic analysis is:

1. initialize problems (if needed)
2. assemble problems
3. combine assemblies to construct a linear system
   ``\boldsymbol{A}\boldsymbol{x} = \boldsymbol{b}``
4. solve linear system, as effectively as possible
5. update solution back to problems or elements
6. postprocess fields
7. write requested results to file
8. check convergence, if not converged, go back to 2.

# Developing solver for LinearSystem

When models gets big enough, step 3. is dominating in a solution process. For
that reason we have abstracit type [`AbstractLinearSystemSolver`](@ref) which
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
using FEMBase: LinearSystem, AbstractLinearSystemSolver
import FEMBase: solve! # hide
```

First setup linear system from [`Problems`]:

```@example FEMBase
K = sparse([  4.0  -1.0  -2.0  -1.0   0.0
	     -1.0   7.0  -4.0  -2.0   0.0
	     -2.0  -4.0  10.0  -1.0  -3.0
	     -1.0  -2.0  -1.0   4.0   0.0
	      0.0   0.0  -3.0   0.0   3.0])
C1 = sparse([ 1.0  0.0  0.0  0.0  0.0
	      0.0  0.0  0.0  0.0  0.0
	      0.0  0.0  0.0  0.0  0.0
	      0.0  0.0  0.0  1.0  0.0
	      0.0  0.0  0.0  0.0  0.0])
C2 = C1
D = spzeros(5, 5)
f = sparsevec([33.0, 33.0, 165.0, 33.0, 132.0])
g = spzeros(5)
u = spzeros(5)
la = spzeros(5)
ls = LinearSystem(K, C1, C2, D, f, g, u, la)
nothing # hide
```

For example, a simple solver for small problems would be to use `lufact` from
UMFPACK:

```@example FEMBase

type LUSolver <: AbstractLinearSystemSolver
    # may contain some solver-spesific settings
end

"""
    solve!(solver::LUSolver, ls::LinearSystem)

Solve linear system using LU factorization. If final system has zero rows,
add 1 to diagonal to make matrix non-singular.
"""
function solve!(solver::LUSolver, ls::LinearSystem)
    A = [ls.K ls.C1; ls.C2 ls.D]
    b = [ls.f; ls.g]

    # add 1.0 to diagonal for any zero rows in system
    p = ones(length(b))
    p[unique(rowvals(A))] = 0.0
    A += spdiagm(p)

    # solve A*x = b using LU factorization and update solution vectors
    x = lufact(A) \ full(b)
    ndofs = length(f)
    ls.u = x[1:ndofs]
    ls.la = x[ndofs+1:end]

    return nothing
end

nothing # hide
```

To test it:
```@repl FEMBase
solver = LUSolver()
solve!(solver, ls)
full(ls.u)
full(ls.la)
```


