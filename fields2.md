# Fields

```@setup FEMBase
using FEMBase
```

From the beginning of a project we had a clear concept in our mind: "everything
is a field". That is, everything can vary temporally and spatially. We think
that constant is just a special case of field which does not vary in temporal
nor spatial direction. Fields can vary in spatial direction, i.e. can be either
constant or variable, and in temporal direction, i.e. can be time variant or
time invariant. From this pondering we can think that there exists four kind of
(discrete) fields:

- discrete, constant, time invariant (DCTI)
- discrete, variable, time invariant (DVTI)
- discrete, constant, time variant (DCTV)
- discrete, variable, time variant (DVTV)

Discrete, in this context, means that field is defined in point-wise in
$1 \ldots n$ locations, from where it is then interpolated to whole domain
using some interpolation polynomials, i.e.
```math
u(\xi, t) = \sum_{i} u_i[t] N_{i}(\xi,t),
```math
where
    $N_{i}(\xi, t)$
is the basis function or interpolation polymial corresponding to $i$^{th}
discrete value and 
    $u_{i}$
is the discrete value.

Then we have continuous fields, which are defined in whole domain, or at least
not point-wise. By following the already used abbreviations, we have four more
fields:

- continuous, constant, time invariant (CCTI)
- continuous, variable, time invariant (CVTI)
- continuous, constant, time variant (DCTV)
- continuous, variable, time variant (CVTV)

Continuous, again in this context, does not mean that field has to be defined
everywhere. It's enough that it's defined in function of spatial and/or temporal
coordinates, i.e. we have $u \equiv u(\xi, t)$, without a some spesific basis
needed to interpolate from discrete values. 

Field itself can be in principle anything. However, usually either scalar,
vector or tensor (matrix). Time does not to have be real, it can be for example
angle of some rotating machine or even complex value. 

From these starting points, we assume that the mentioned field system can
describe all imaginable situations.

## Creating new fields

For discrete fields that are varying in spatial direction, value for each
discrete point is defined using NTuple. The order of points is implicitly
assumed to be same than node ordering in ABAQUS. That is, first corner nodes
in anti-clockwise direction and after that middle nodes.

For example, `(1, 2, 3, 4)` is a scalar field having length of 4 and
`([1,2],[2,3],[3,4],[4,5])` is a vector field having length of 4.

For fields that are varying in temporal direction, `time => value` syntax is
used. The first item in pair is time (or similar) and second item is value 
assigned to that time. For example, `0.0 => 1.0` is a time-dependent scalar
field having value 1.0 at time 0.0.

The most simple field is a field that is constant in both time and spatial
direction. Discrete, constant, time invariant [`DCTI`](@ref):

```@repl FEMBase
a = DCTI(1.0)
```

Then we have discrete, variable, time invariant fields [`DVTI`](@ref). Notice
the use of `Tuple` when defining field.


```@repl FEMBase
b = DVTI( (1.0, 2.0) )
```

Discrete, constant, time variant field [`DCTV`](@ref) is constant in spatial
direction $\partial u/\partial x = 0$ but can vary in temporal direction,
$\partial u/\partial t\neq 0$. Here, `=>` syntax is used. New values can be
added to field using function [`update!`](@ref). If there already exists a
value for that particular time, it will be overridden. It is assumed that
content of field in time direction is monotonically increasing, i.e.
```math
t_{i-1} < t_i < t_{i+1}.
```

For the sake of clarity let's also mention that `update!` works for time
invariant fields as well if content needs to be updated.

```@repl FEMBase
c = DCTV(0.0 => 1.0, 1.0 => 2.0)
update!(c, 2.0 => 3.0)
```

Discrete, variable, time variant [`DVTV`](@ref) field is the most general
one, allowing values of field to vary in both spatial and time direction.

```@repl FEMBase
d = DVTV(0.0 => (1.0, 2.0), 1.0 => (2.0, 3.0))
update!(d, 2.0 => (3.0, 4.0))
```

In examples above, all fields defined was scalar fields. Defining vector or
tensor fields goes in the same way. The only difference is that now we define
vectors and tensors instead of a single scalar value. They can vary in spatial
and time direction in the same way than scalar fields. Here is example of
defining the abovementioned vector fields:

```@repl FEMBase
a = DCTI([1.0, 2.0])
b = DVTI(([1.0, 2.0], [2.0, 3.0]))
c = DCTV(0.0 => [1.0, 2.0], 1.0 => [2.0, 3.0])
d = DVTV(0.0 => ([1.0, 2.0], [2.0, 3.0]), 1.0 => ([2.0, 3.0], [3.0, 4.0]))
```

## Accessing fields

Accessing fields in time direction is done using a function [`interpolate`](@ref).
For example, if we have (constant) $[1,2]$ at time $t=0.0$ and $[3,4]$ at time
$t=1.0$, linear interpolation in time direction yields

```@repl FEMBase
c = DCTV(0.0 => [1.0,2.0], 1.0 => [3.0,4.0])
interpolate(c, 0.5)
```

If field is spatially varying, a `Tuple` will be returned, having one value
for each "node". This can then be interpolated in spatial direction, typically
using basis functions defined on element, i.e. $u = N_{i} u_{i}$:

```@repl FEMBase
d = DVTV(0.0 => (1.0,2.0), 1.0 => (3.0,4.0))
interpolate(d, 0.5)
```

Although the two fields above looks quite same, the key difference is that in
DCTV field we have a constant vectorial value (defined using square brackets [])
and in latter DVTV field we have a scalar value (defined using round brackets)
changing in spatial direction from 1.0 to 2.0 at time $t=0.0$ and changing from
3.0 to 4.0 at time $t=1.0$.

If a field is time invariant, interpolation in time direction returns a trivial
solution:

```@repl FEMBase
interpolate(DCTI(1.0), 0.5)
interpolate(DVTI((1.0,2.0)), 0.5)
```

For spatially varying fields, one can access to ith element using getindex:

```@repl FEMBase
a = DVTI((1.0,2.0))
getindex(a, 1)
```

For field varying both temporally and spatially, one has first to interpolate
in time direction to get iterable tuple:

```@repl FEMBase
d = DVTV(0.0 => (1.0,2.0), 1.0 => (3.0,4.0))
result = interpolate(d, 0.5)
getindex(result, 1)
```

Internally, each field is a subtype of [`AbstractField`](@ref) having a field
`data`, which be accessed directly for hacking purposes.

```@repl FEMBase
d.data
```

## Continuous fields

Continuous fields may be useful when defining analytical boundary conditions.
For that we have [`CVTV`](@ref), where "C" stands for continuous.

```@repl FEMBase
f(xi,t) = xi[1]*xi[2]*t
g = CVTV(f)
g((1.0,2.0), 3.0)
```

## Dictionary fields

Usually it is assumed that size of length of discrete field matches to the
number of basis functions of a single element, typically something small like
1-10.

However, there might be cases where it is more practical to define field in a
sense that indexing is not continuous or starting from 1. For example, we might
want to define field common for a set of elements. In that case it's natural to
think that each index in field corresponds to the certain id-number of node. For
example, if we have a triangle element connecting nodes 1, 1000 and 100000, we
still want to access that field naturally using `getindex`, e.g. `f[1]`, `f[1000]`
and `f[100000]`. In that case, more appropriate internal structure for field is
based on a dictionary, not tuple.

It only makes sense to define dictionary fields for spatially varying fields.
Two new fields are introduced: [`DVTId`](@ref) and [`DVTVd`](@ref), where the
last character "d" stands for "dictionary".

Keep on mind, that this type of field has one restriction. If and when this
field is typically defined on nodes of several elements, field must be continuous
between elements. That is, if field value in node 1000 is for example 1.0, then
it's 1.0 in all elements connecting to that node. To define jumps on field, one
must define field element-wise. 

Define eg. "geometry" for nodes 1,1000,100000:

```@repl FEMBase
X = Dict(1=>[0.0,0.0], 1000=>[1.0,0.0], 100000=>[1.0,1.0])
G = DVTId(X)
G[1], G[1000], G[100000]
```

Interpolation in time directions works in a same way than with other fields
depends from time.

```@repl FEMBase
Y = Dict(1=>[1.0,1.0], 1000=>[2.0,1.0], 100000=>[2.0,2.0])
F = DVTVd(0.0 => X, 1.0 => Y)
interpolate(F,0.5)[100000]
```

## Using common constructor `field`

Now we have introduced total of 7 fields: DCTI, DCTV, DVTI, DVTV, CVTV, DVTId,
DVTVd. A good question arises that how to remember all this stuff and is it
even necessary? Luckily not, because one can use a single constructor called
[`field`](@ref) to create all kind of fields. Type of field is inspected from
data type. It's not necessary to remember all this technical stuff, just declare
new field using more of less intuitive syntax and [`field`](@ref)-function.

```@repl FEMBase
f1 = field(1)
f2 = field(1, 2)
f3 = field(0.0 => 1)
f4 = field(0.0 => (1, 2), 1.0 => (2, 3))
f5 = field((xi,t) -> xi[1]*t)
f6 = field(1 => 1, 2 => 2)
f7 = field(0.0 => (1=>1, 10=>2), 1.0 => (1=>2,10=>3))
```

## Developing new fields

If the FEMBase ones are not enough, it's always possible to define new ones.
Minimum requirements is that field is a subtype of `AbstractField` and
`interpolate`, `getindex`, has been defined to it. Field can, for example
fetch data from random.org or market stocks, read data from hard drive or add
some stochastics behavior to it.

## Functions and types related to fields

```@meta
CurrentModule = FEMBase
DocTestSetup = quote
    using FEMBase
end
```

### Types

```@docs
AbstractField
DCTI
DVTI
DCTV
DVTV
CVTV
DVTId
DVTVd
```

### Functions (internal)

These functions needs to be defined when developing new fields:

```@docs
new_field
update_field!
interpolate_field
```


### Functions (public)

```@docs
field(x)
update!(field::F, data) where {F<:AbstractField}
interpolate(field::F, time) where {F<:AbstractField}
```
