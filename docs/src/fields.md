

```julia
using FEMBase
```

From the beginning a clear concept: "everything is a field". We think that constant is just a special case of field which does not vary in temporal and spatial direction. Fields can vary in spatial direction (constant or variable) and in temporal direction (time variant and in-variant).

From here we can deduct for types of (discrete) fields:

- discrete, constant, time invariant (DCTI)
- discrete, variable, time invariant (DVTI)
- discrete, constant, time variant (DCTV)
- discrete, variable, time variant (DVTV)

Field itself can be anything. However, usually either scalar, vector or tensor (matrix).

# Creating fields

For discrete fields that are varying in spatial direction, value for each discrete point is defined using NTuple. The order of points is implicitly assumed to be same than node ordering in ABAQUS. Usual use case for variable field is that the field is interpolated to the element volume using interpolation polynomials, i.e. $u(\xi) = u_i N_i(\xi)$, where $N_i$ is the basis function for that node and $u_i$ is the discrete value.

For example, `(1, 2, 3, 4)` is a scalar field having length of 4 and `([1,2],[2,3],[3,4],[4,5])` is a vector field having length of 4.

For fields that are varying in temporal direction, `time => value` syntax is used. The first item in pair is time and second item is value attached to that time. For example, `0.0 => 1.0` is a time-dependent scalar field having value 1.0 at time 0.0.

The most simple field is a field that is constant in both time and spatial direction. Discrete, constant, time invariant (DCTI):


```julia
a = DCTI(1.0)
```




    FEMBase.DCTI{Float64}(1.0)



Then we have discrete, variable, time invariant fields (DVTI). Note the use of `NTuple` when defining field.


```julia
b = DVTI( (1.0, 2.0) )
```




    FEMBase.DVTI{2,Float64}((1.0, 2.0))



Discrete, constant, time variant field (DCTV) is constant in spatial direction but can vary in temporal direction. Here, `=>` syntax is used. New values can be added to field using `update!`. If there already exists value for that particular time, it will be overridden. It is assumed that content of field in time direction is monotonically increasing, i.e. $t_{i-1} < t_i < t_{i+1}$. For the sake of clarity let's also mention that `update!` works for time invariant fields as well if content needs to be updated.


```julia
c = DCTV(0.0 => 1.0, 1.0 => 2.0)
update!(c, 2.0 => 3.0)
c
```




    FEMBase.DCTV{Float64}(Pair{Float64,Float64}[0.0=>1.0, 1.0=>2.0, 2.0=>3.0])



Discrete, variable, time variant (DVTV) field is the most general one, allowing values of field to vary in both spatial and time direction.


```julia
d = DVTV(0.0 => (1.0, 2.0), 1.0 => (2.0, 3.0))
update!(d, 2.0 => (3.0, 4.0))
```




    3-element Array{Pair{Float64,Tuple{Float64,Float64}},1}:
     0.0=>(1.0, 2.0)
     1.0=>(2.0, 3.0)
     2.0=>(3.0, 4.0)



In examples above, all fields was scalar fields. Defining vector or tensor fields goes in the same spirit. Only difference is that now we define vectors and tensors, not a single scalar value. They can vary in spatial and time direction in the same way than scalar fields. Here is example of defining the abovementioned vector fields:


```julia
a = DCTI([1.0, 2.0])
b = DVTI(([1.0, 2.0], [2.0, 3.0]))
c = DCTV(0.0 => [1.0, 2.0], 1.0 => [2.0, 3.0])
d = DVTV(0.0 => ([1.0, 2.0], [2.0, 3.0]), 1.0 => ([2.0, 3.0], [3.0, 4.0]))
```




    FEMBase.DVTV{2,Array{Float64,1}}(Pair{Float64,Tuple{Array{Float64,1},Array{Float64,1}}}[0.0=>([1.0, 2.0], [2.0, 3.0]), 1.0=>([2.0, 3.0], [3.0, 4.0])])



# Accessing fields

Accessing fields is done using a single command: `interpolate`. For time varying fields, one can interpolate in time direction. For example, if we have (constant) $[1,2]$ at time $t=0.0$ and $[3,4]$ at time $t=1.0$, linear interpolation yields


```julia
c = DCTV(0.0 => [1.0,2.0], 1.0 => [3.0,4.0])
interpolate(c, 0.5)
```




    2-element Array{Float64,1}:
     2.0
     3.0



If field is spatially varying, a `Tuple` will be returned, typically having one value for each "node". This can then be interpolated in spatial direction, typically using basis functions defined on element, e.g. $u = N_i u_i$:


```julia
d = DVTV(0.0 => (1.0,2.0), 1.0 => (3.0,4.0))
interpolate(d, 0.5)
```




    (2.0, 3.0)



Although the two fields above looks quite same, the key difference is that in DCTV field we have a constant vectorial value (defined using square brackets []) and in latter DVTV field we have a scalar value (defined using round brackets) changing in spatial direction from 1.0 to 2.0 at time $t=0.0$ and changing from 3.0 to 4.0 at time $t=1.0$.

One should be always able to interpolate in time direction, even if field is time invariant, to get trivial solution:


```julia
interpolate(a, 0.5), interpolate(b, 0.5),
interpolate(c, 0.5), interpolate(d, 0.5)
```




    ([1.0, 2.0], ([1.0, 2.0], [2.0, 3.0]), [2.0, 3.0], (2.0, 3.0))



For spatially varying fields, one can access to ith element using getindex:


```julia
getindex(a, 1), getindex(b, 1), getindex(b, 2)
```




    ([1.0, 2.0], [1.0, 2.0], [2.0, 3.0])



First time interpolation, then spatial lookup, i.e.


```julia
getindex(interpolate(d, 0.5), 1)
```




    2.0



By passing already evaluated basis functions to `interpolate`, field is intepolated in both temporal and spatial direction.


```julia
d = DVTV(0.0 => (1.0,2.0), 1.0 => (3.0,4.0))
xi = 0.0
N = [1/2*(1-xi) 1/2*(1+xi)] # evaluated basis functions, linear interpolation
t = 1.0 # time
interpolate(d, t, N)
```




    3.5



Internally, each field is a subtype of `AbstractField`. Internally each field has a `data` which be accessed directly for hacking.


```julia
d.data
```




    2-element Array{Pair{Float64,Tuple{Float64,Float64}},1}:
     0.0=>(1.0, 2.0)
     1.0=>(3.0, 4.0)



# Continuous fields

Then we have continuous fields which may be useful when defining analytical boundary conditions. For that we have CVTV, where "C" stands for continuous.


```julia
f(xi,t) = xi[1]*xi[2]*t
g = CVTV(f)
g([1.0,2.0],3.0)
```




    6.0



# Dictionary fields

Usually it is assumed that discrete field has values number of nodes or degrees of freedom in element, typically something small like 1-30. There might be cases where it is more practical to define field in a sense that indexing is not continuous or starting from 1. For example, we might want to define field common for a set of elements. In that case it's natural to think that each index in field corresponds to the certain id-number of node. For example, if we have triangle element connecting nodes 1, 1000 and 100000, we still want to access to that field naturally using `getindex`, e.g. `f[1]`, `f[1000]` and `f[100000]`. In that case, more appropriate internal structure for field is based on dictionary.

It only makes sense to define dictionary fields for spatially varying fields. Two new fields are introduced: DVTId and DVTVd, where last "d" stands for "dictionary".

Keep on mind, that this type of field has one restriction. If and when this field is typically defined on nodes of several elements, field must be continuous between elements. That is, if field value in node 1000 is for example 1.0, then it's 1.0 in all elements connecting to that node. To define jumps on field, one must define field element-wise. 


```julia
# Define eg. "geometry" for nodes 1,2,3,4.
X = Dict(1=>[0.0,0.0], 1000=>[1.0,0.0], 100000=>[1.0,1.0])
G = DVTId(X)
G[1], G[1000], G[100000]
```




    ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0])




```julia
Y = Dict(1=>[1.0,1.0], 1000=>[2.0,1.0], 100000=>[2.0,2.0])
F = DVTVd(0.0 => X, 1.0 => Y)
interpolate(F,0.5)[100000]
```




    2-element Array{Float64,1}:
     1.5
     1.5



# Using common constructor `field`

Now we have introduced total of 7 fields: DCTI, DCTV, DVTI, DVTV, CVTV, DVTId, DVTVd. A good question arises that how to remember all this stuff and is it even necessary? Luckily not, because one can use a single constructor called `field` to create all kind of fields. Type of field is inspected from data type. You really don't have to know about this technical stuff, just declare new field using intuition and `field` command.


```julia
f1 = field(1.0)
typeof(f1)
```




    FEMBase.DCTI{Float64}




```julia
f2 = field((1.0, 2.0))
typeof(f2)
```




    FEMBase.DVTI{2,Float64}




```julia
f3 = field(0.0 => 1.0)
typeof(f3)
```




    FEMBase.DCTV{Float64}




```julia
f4 = field(0.0 => (1.0, 2.0), 1.0 => (2.0, 3.0))
typeof(f4)
```




    FEMBase.DVTV{2,Float64}




```julia
f5 = field((xi,t) -> xi[1]*t)
typeof(f5)
```




    FEMBase.CVTV




```julia
f6 = field(1 => 1.0, 2 => 2.0)
typeof(f6)
```




    FEMBase.DVTId{Float64}




```julia
f7 = field(0.0 => (1=>1.0,10=>2.0), 1.0 => (1=>2.0,10=>3.0))
typeof(f7)
```




    FEMBase.DVTVd{Float64}



# Developing new fields

If the default ones are not enough, it's always possible to define new ones. Minimum requirements is that field is a subtype of `AbstractField` and `interpolate` has been defined to it.
