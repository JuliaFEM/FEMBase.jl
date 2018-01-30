# Problems

The role of problems in JuliaFEM is to work as a container for a set of elements.
They contain elements and an information how the elements are assembled to the
global assembly. The key point is that thanks to multiple dispatch, each problem
defines also the physics behind discretization and all problems are assembled
using a command [`assemble_elements!`](@ref). 

As an example, a heat equation in two dimensions is discretized. Mathematically
known also as Poisson problem. Strong form of the problem is
```math
 \nabla \cdot (-\nabla u ) = f    \text{ in } \Omega, \\
                         u = u_0  \text{ on } \Gamma_\mathrm{D}, \\
      \frac{u}{\partial n} = g    \text{ on } \Gamma_\mathrm{N},
```
and corresponding variational problem is to find ``u\in\mathcal{U}`` such that
for all ``v\in\mathcal{V}``
```math
\int_{\Omega}k\nabla u\cdot\nabla v\,\mathrm{d}x=\int_{\Omega}fv\,\mathrm{d}x+\int_{\Gamma_{\mathrm{N}}}gv\,\mathrm{d}s
```

Let's call ``k`` thermal conductivity, ``f`` heat source and ``g`` heat flux.

```@setup FEMBase
using FEMBase
import FEMBase: assemble_elements!, get_unknown_field_name
```

Each new problem must be a subtype of `FieldProblem` or `BoundaryProblem`. The
main difference between these two are that `FieldProblem` is assembling local
matrices for domain ``\Omega`` whereas `BoundaryProblem` is creating (in general)
constraint matrices for boundary ``\Gamma_{\mathrm{D}}``. The general structure
to solve in JuliaFEM is currently described by four different matrices and two
vectors, i.e.
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

We believe that this construction is general enough to describe all possible
situations in future. Quite often ``\boldsymbol{C}_1 = \boldsymbol{C}_2^\mathrm{T}``
and ``\boldsymbol{D} = \boldsymbol{0}`` so that we have a typical saddle point
problem
```math
\begin{bmatrix}\boldsymbol{K} & \boldsymbol{C}^{\mathrm{T}}\\
\boldsymbol{C} & \boldsymbol{0}
\end{bmatrix}\begin{bmatrix}\boldsymbol{u}\\
\boldsymbol{\lambda}
\end{bmatrix}=\begin{bmatrix}\boldsymbol{f}\\
\boldsymbol{g}
\end{bmatrix},
```
which is equivalent to minimization problem
```math
\min_{\boldsymbol{C}\boldsymbol{u} = \boldsymbol{g}}
\frac{1}{2} \boldsymbol{u}^{\mathrm{T}} \boldsymbol{K} \boldsymbol{u} -
\boldsymbol{u}^{\mathrm{T}} \boldsymbol{f}.
```

## Discretizing field problem

So, first we must define a new type, e.g. `Heat`, which is a subclass of
[`FieldProblem`](@ref). Problem-wide parameters can be defined into struct
if needed.

```@example FEMBase
type Heat <: FieldProblem
end
```

In principle it's possible to assemble one element at a time, but way more memory
efficient is to assemble all elements with same kind (basis) at same time and
preallocate memory only one time before looping through element list. Implementation
for assembling local stiffness matrix is

```@example FEMBase
function assemble_elements!{B}(problem::Problem{Heat}, assembly::Assembly,
                               elements::Vector{Element{B}}, time::Float64)

    println("Assembling elements of kind $B")
    bi = BasisInfo(B)
    ndofs = length(B)
    Ke = zeros(ndofs, ndofs)

    for element in elements
        fill!(Ke, 0.0)
        for ip in get_integration_points(element)
            J, detJ, N, dN = element_info!(bi, element, ip, time)
            k = element("thermal conductivity", ip, time)
            Ke += ip.weight * k*dN'*dN * detJ
        end
        gdofs = get_gdofs(problem, element)
        add!(assembly.K, gdofs, gdofs, Ke)
    end

end

nothing # hide
```

Here, first some memory is allocated to calculate Jacobian, gradients etc. to
`BasisInfo`. `Ke` is used to store local stiffness matrix. Then iterate over
all elements and integration points, evaluate basis and add contribution to
local stiffness matrix. Finally, get global degrees of freedom for element by
using command `get_gdofs` and finally add contribution to global stiffness
matrix `K`.

From performance point of view, it's important that memory allocations inside
at least the innermost loop is minimized, although assembling global stiffness
matrix can be parallelized (at least almost) perfectly and is not considered
as a bottleneck when models get big enough. It's anyway a good idea to pay
attention to the memory allocations.

### Setting and getting global degrees of freedom for element

[`get_gdofs`](@ref) is returning the global degrees of freedom for element. They
can be set manually using [`set_gdofs(problem, element, dofs)`](@ref). Otherwise
they are calculated automatically based on the problem dimension using formula
`g(i,j,d) = d*(i-1)+j`, where `i` is local node number, `j` is the number of
degree of freedom and `d` is the maximum number of degrees of freedom for node.
With this formula dofs are ordered so that first comes all dofs for node 1,
then for node 2 and so on. For 3 dofs/node we get ``(u_{11}, u_{12}, u_{13},
u_{21}, u_{22}, u_{23}, \ldots, u_{N1}, u_{N2}, u_{N3})``, where the first index
is node number and second index is dof number.

Let's create some test problem and test our implementation:

```@example FEMBase
el1 = Element(Quad4, [1, 2, 3, 4])
X = Dict(1 => [0.0,0.0], 2 => [1.0,0.0], 3 => [1.0,1.0], 4 => [0.0,1.0])
update!(el1, "geometry", X)
update!(el1, "thermal conductivity", 6.0)
elements = [el1]
assembly = Assembly()
problem = Problem(Heat, "test problem", 1)
nothing # hide
```

Now the struct `Heat` is empty. If we need to store some problem-wide settings
to that struct, they can be found from `problem.properties`. When creating a
new [`Problem`](@ref), the syntax is `Problem(type, name, field_dimension)`, where
two first arguments are self descriptive. The third one is the information,
how many degrees of freedom is in this problem. As Poisson problem is scalar
equation, there is only 1 degrees of freedom per node. For example in continuum
mechanics, where the unknown field is displacement, there is usually 2-6 degrees
of freedom per node, depending on problem type. Next we do the actual assembling
of elements into global stiffness matrix:

```@example FEMBase
time = 0.0
assemble_elements!(problem, assembly, elements, time)
```

```@repl FEMBase
full(assembly.K)
```

There is actually one [`Assembly`](@ref) inside [`Problem`](@ref) and elements
are defined to problem using [`add_elements!`](@ref), so a more realistic use
case to create global assembly would be to use [`assemble!(problem, time)`](@ref)
as shown below:

```@example FEMBase
el2 = Element(Tri3, [3, 2, 5])
X[5] = [2.0, 1.0]
elements = [el1, el2]
update!(elements, "geometry", X)
update!(elements, "thermal conductivity", 6.0)
problem = Problem(Heat, "test problem", 1)
add_elements!(problem, elements)
assemble!(problem, time)
```

```@repl FEMBase
full(problem.assembly.K)
```

Now, function defined above is actually executed two times, first for elements
using `Tri3` basis and after that for `Quad4`. That is, `assemble!(problem, time)`
is grouping elements by their type and calling function for each element type
separately. It also does some initializations and gives possibility to mangle
matrices before and after assembly. These hacks may be useful if one needs to
add some discrete values to the matrices after assembly or e.g. save matrices
to disk for later diagnoses.

We also need to deal with the integrals on the right hand side. The first
integral is done over the domain and can be included to the same
`assemble_elements!`-function than stiffness matrix. Boundary term can be
handled in different ways. One option is to define it yet in same function
and search for fields like `surface traction S1`, where `S1` is one side of
the element. This is how it is done in ABAQUS. Another way is to use lower
dimensional "boundary element" in assembly and add surface term to that element.
This is how it is done in Code Aster.

`assembly_elements!`-function defined above can be overridden by restricting
the type of elements list, `elements::Vector{Element{B}}` to a some spesific
elements. This allows, for example, to optimize assembly for some certain
element type what is commonly used. Another use case is to define different
assembly function for boundary elements. In 2d setting, voluminal elements like
Tet3, Tet6, Quad4, Quad8, Quad9 are integrated over volume and they one
dimensional counterparts Seg2, Seg3 can be used to assign boundary fluxes.

The updated `assemble_elements!`-function, which can also handle volume load
from right hand side of the equation, i.e. 
```math
\int fv\,\mathrm{d}x,
```
looks like following:

```@example FEMBase
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

nothing # hide
```

At last, implement boundary elements to handle heat flux. To choose what elements
should use this assembly function, `elements::Vector{Element{B}}` must be restricted
only to group where `B<:Union{Seg2, Seg}`.

```@example FEMBase
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

nothing # hide
```

To test everything implemented, create some small test problem:

```@example FEMBase
el1 = Element(Quad4, [1, 2, 3, 4])
el2 = Element(Tri3, [3, 2, 5])
el3 = Element(Seg2, [3, 5])
X = Dict(1 => [0.0, 0.0],
         2 => [1.0, 0.0],
         3 => [1.0, 1.0],
         4 => [0.0, 1.0],
         5 => [2.0, 1.0])
elements = [el1, el2, el3]
update!(elements, "geometry", X)
update!(elements, "thermal conductivity", 6.0)
update!(el3, "heat flux", 264.0)
update!(el1, "heat source", 132.0)
problem = Problem(Heat, "test problem", 1)
add_elements!(problem, elements)
assemble!(problem, time)
```

Global stiffness matrix ``\boldsymbol{K}`` and force vector ``\boldsymbol{f}``
are

```@repl FEMBase
K = full(problem.assembly.K)
f = full(problem.assembly.f)
```

To get unique solution, some essential boundary condition must be given, e.g.
set dofs 1 and 4 fixed, ``u_1 = u_4 = 0``.

```@repl FEMBase
u = zeros(5)
all_dofs = collect(1:5)
fixed_dofs = [1, 4]
free_dofs = setdiff(all_dofs, fixed_dofs)
u[free_dofs] = K[free_dofs,free_dofs] \ f[free_dofs]
u
```

## Discretizing boundary problem

Can be e.g. Dirichlet boundary, contact / mesh tie problem between two domains,
BEM formulation, kinematic coupling (MPC) and so on.

A good question is to determine how to map element local degrees of freedom
with global ones. As the plan is to develop a multiphysical FEM platform, it's
a hard question how to create this kind of mapping in a dynamic way. Yet
another question is how to set boundary conditions for different type of
physics. For now, these questions are answered by defining the name of unknown
field using function `get_unknown_field_name` and giving the dimension of unknown
field when creating a problem, so that change of information between two problem
is possible. This solution has already identified shortcomings and can be expected
to change in future.

```@example FEMBase
function get_unknown_field_name(::Problem{Heat})
    return "temperature"
end

nothing # hide
```

```@example FEMBase
type DirichletBC <: BoundaryProblem
end
```

Our implementation to handle boundary condition ``u = u_0`` looks following:

```@example FEMBase
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

To fix dofs 1 and 4 like before:

```@example FEMBase
bel1 = Element(Seg2, [1, 4])
update!(bel1, "geometry", X)
update!(bel1, "temperature 1", 0.0)
bc = Problem(DirichletBC, "fixed", 1, "temperature")
add_elements!(bc, [bel1])
assemble!(bc, time)
```

Now we have all matrices needed:

```@repl FEMBase
C1 = full(bc.assembly.C1, 5, 5)
C2 = full(bc.assembly.C2, 5, 5)
D = full(bc.assembly.D, 5, 5)
g = full(bc.assembly.g, 5, 1)
```

Together with already calculated matrices, we can now form saddle point problem
``\boldsymbol{A}\boldsymbol{x} = \boldsymbol{b}``:

```@repl FEMBase
A = [K C1; C2 D]
b = [f; g]
```

Solution:

```@repl FEMBase
nz = [1, 2, 3, 4, 5, 6, 9]
x = zeros(10)
x[nz] = A[nz,nz] \ b[nz]
```

As a result we have found ``\boldsymbol{u}`` and ``\boldsymbol{\lambda}``:

```@repl FEMBase
u = x[1:5]
la = x[6:10]
u' * la
```

Both field problem and boundary problems can of course have all four matrices
and two vectors. For example, in finite sliding contact algorithms all four
matrices are used as algorithms are contributing to stiffness matrix also when
linearized properly.

## Assembling mass matrices
## Using problem-wide fields
## Discretizing mixed problems
## Using automatic differentiation to linearize non-linear problem
