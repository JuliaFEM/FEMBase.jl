# Elements

In JuliaFEM, elements can be considered as "containers", combining fields and
basis functions described above. Among that, element has information about
topology (connectivity) and numerical integration rule. These fundamentals
forms a finite element, the backbone of finite element method.

```@setup FEMBase
using FEMBase
```

New elements are constructed by passing basis type (e.g. Seg2, Quad4, Tet10, ...)
to `Element` and list of id numbers to where element is topologically connected.

```@repl FEMBase
el = Element(Quad4, [1, 2, 3, 4])
```

Setting fields to element is done using a command `update!`, which either
creates a new field if field with that name does not already exist, or updates
the old one. Typically, at least field called `geometry` needs to be defined to
element as it's used to calculate Jacobian of element. Fields can be discrete,
continuous, time invariant or variant, variable or constant, or anything else
that is subclassed from `AbstractField`.

```@repl FEMBase
X = Dict(1 => [0.0,0.0], 2=>[1.0,0.0], 3=>[1.0,1.0], 4=>[0.0,1.0])
update!(el, "geometry", X)
```

Internally, fields are stored in a `Dict`:

```@repl FEMBase
el.fields
```

The command `update!` is automatically inspecting field type based in input.
For example, to define temporal field varying spatially:

```@repl FEMBase
u0 = ([0.0,0.0], [0.0,0.0], [0.0,0.0], [0.0,0.0])
u1 = ([0.0,0.0], [0.0,0.0], [0.5,0.0], [0.0,0.0])
update!(el, "displacement", 0.0 => u0)
update!(el, "displacement", 1.0 => u1)
el.fields
```

Interpolating of fields can be done using syntax `Element(field_name, xi, time)`.
For example, position of material particle $X$ in initial configuration and
deformed configuration in the middle of the element at time $t=1$ can be found as

```@repl FEMBase
xi = (0.0, 0.0)
time = 1.0
X = el("geometry", xi, time)
u = el("displacement", xi, time)
x = X+u
```

Jacobian, determinant of Jacobian and gradient of field can be calculated adding
extra argument `Val{:Jacobian}`, `Val{:detJ}`, `Val{:Grad}` to the above command
and not passing field name, i.e.

```@repl FEMBase
el(xi, time, Val{:Jacobian})
el(xi, time, Val{:detJ})
el(xi, time, Val{:Grad})
```

Usually what is needed when calculating local stiffness matrices is a gradient
of some field. For example, displacement gradient and temperature gradient can
be obtained following way:

```@repl FEMBase
gradu = el("displacement", xi, time, Val{:Grad})
update!(el, "temperature", (1.0, 2.0, 3.0, 4.0))
gradT = el("temperature", xi, time, Val{:Grad})
```

Accessing integration points of element is done using function
`get_integration_points`. Combining interpolation and integration one can
already calculate local matrices of a single element or, for example area and
strain energy:

```@example FEMBase
update!(el, "lambda", 96.0)
update!(el, "mu", 48.0)

A = 0.0
W = 0.0
for ip in get_integration_points(el)
    detJ = el(ip, time, Val{:detJ})
    A += ip.weight * detJ
    ∇u = el("displacement", ip, time, Val{:Grad})
    E = 1/2*(∇u + ∇u')
    λ = el("lambda", ip, time)
    μ = el("mu", ip, time)
    W += ip.weight * ( λ/2*trace(E*E') + μ*trace(E)^2) * detJ
end

println("Area: $A")
println("Strain energy: $W")
```

To calculate local stiffness matrix for Poisson problem:

```@example FEMBase
K = zeros(4,4)
update!(el, "coefficient", 36.0)
for ip in get_integration_points(el)
    dN = el(ip, time, Val{:Grad})
    detJ = el(ip, time, Val{:detJ})
    c = el("coefficient", ip, time)
    K += ip.weight * c*dN'*dN * detJ
end
K
```

For more memory efficient code it's a good idea to use `BasisInfo` and
`element_info!` which allocates working memory to calculate all "basic stuff"
for a single integration point, like Jacobian, determinant of Jacobian, basis
and it's partial derivatives with respect to reference configuration $X$.

```@example FEMBase
bi = BasisInfo(Quad4)
fill!(K, 0.0)
for ip in get_integration_points(el)
    J, detJ, N, dN = element_info!(bi, el, ip, time)
    c = el("coefficient", ip, time)
    K += ip.weight * c*dN'*dN * detJ
end
K
```

## Using analytical fields

## Creating fields depending from other fields
