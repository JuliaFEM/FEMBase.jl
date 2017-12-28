
In JuliaFEM, elements are "containers", combining fields and basis functions described above. Among that, element has information about topology (connectivity) and integration rule. These fundamentals forms a finite element, the backbone of finite element method, as the basic idea after all is just to discretize continuous domain to smaller topological entities like tetrahedrons and perform same operations to each element.


```julia
using FEMBase
```


```julia
el = Element(Quad4, [1, 2, 3, 4])
```




    FEMBase.Element{FEMBasis.Quad4}(-1, [1, 2, 3, 4], FEMBase.Point{FEMBase.IntegrationPoint}[], Dict{String,FEMBase.AbstractField}(), FEMBasis.Quad4())



Setting fields to element is done using a command `update!`, which either creates new field if does not already exist, or updates the old one. Typically, at least field called `geometry` needs to be defined to element as it's used to calculate Jacobian of element. Fields can be discrete, continuous, time invariant or variant, variable or constant, like described earlier.


```julia
X = Dict(1 => [0.0,0.0], 2=>[1.0,0.0], 3=>[1.0,1.0], 4=>[0.0,1.0])
update!(el, "geometry", X)
```




    FEMBase.DVTId{Array{Float64,1}}(Dict(4=>[0.0, 1.0],2=>[1.0, 0.0],3=>[1.0, 1.0],1=>[0.0, 0.0]))




```julia
el.fields
```




    Dict{String,FEMBase.AbstractField} with 1 entry:
      "geometry" => FEMBase.DVTId{Array{Float64,1}}(Dict(4=>[0.0, 1.0],2=>[1.0, 0.0…




```julia
u0 = ([0.0,0.0], [0.0,0.0], [0.0,0.0], [0.0,0.0])
u1 = ([0.0,0.0], [0.0,0.0], [0.5,0.0], [0.0,0.0])
update!(el, "displacement", 0.0 => u0)
update!(el, "displacement", 1.0 => u1)
el.fields
```




    Dict{String,FEMBase.AbstractField} with 2 entries:
      "geometry"     => FEMBase.DVTId{Array{Float64,1}}(Dict(4=>[0.0, 1.0],2=>[1.0,…
      "displacement" => FEMBase.DVTV{4,Array{Float64,1}}(Pair{Float64,NTuple{4,Arra…



Interpolating of fields goes calling `Element(field_name, xi, time)`. For example, position of material particle $X$ in initial configuration and deformed configuration in the middle of the element at time $t=1$ can be found as


```julia
xi = (0.0, 0.0)
time = 1.0
X = el("geometry", xi, time)
u = el("displacement", xi, time)
x = X+u
println("X = $X, x = $x")
```

    X = [0.5, 0.5], x = [0.625, 0.5]
    

Jacobian, determinant of Jacobian and gradient of field can be calculated adding extra argument `Val{:Jacobian}`, `Val{:detJ}`, `Val{:Grad}` to the above command and not passing field name, i.e.


```julia
el(xi, time, Val{:Jacobian})
```




    2×2 Array{Float64,2}:
     0.5  0.0
     0.0  0.5




```julia
el(xi, time, Val{:detJ})
```




    0.25




```julia
el(xi, time, Val{:Grad})
```




    2×4 Array{Float64,2}:
     -0.5   0.5  0.5  -0.5
     -0.5  -0.5  0.5   0.5



Usually what the user wants is still a gradient of some field. For example, displacement gradient:


```julia
gradu = el("displacement", xi, time, Val{:Grad})
gradu
```




    2×2 Array{Float64,2}:
     0.25  0.25
     0.0   0.0 



Or temperature gradient:


```julia
update!(el, "temperature", (1.0, 2.0, 3.0, 4.0))
gradT = el("temperature", xi, time, Val{:Grad})
```




    1×2 RowVector{Float64,Array{Float64,1}}:
     0.0  2.0



Accessing integration points of element is done using command `get_integration_points`. Combining interpolation and integration one can already calculate local matrices of a single element or, for example area and strain energy:


```julia
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

    Area: 1.0
    Strain energy: 10.0
    

Local stiffness matrix for Poisson problem:


```julia
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




    4×4 Array{Float64,2}:
      24.0   -6.0  -12.0   -6.0
      -6.0   24.0   -6.0  -12.0
     -12.0   -6.0   24.0   -6.0
      -6.0  -12.0   -6.0   24.0


