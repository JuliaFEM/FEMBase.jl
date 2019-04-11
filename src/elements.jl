# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

abstract type AbstractElement{T<:FEMBasis.AbstractBasis} end

mutable struct Element{T} <: AbstractElement{T}
    id :: Int
    connectivity :: Vector{Int}
    integration_points :: Vector{IP}
    dfields :: Dict{Symbol, AbstractField}
    properties :: T
end

"""
    Element(topology, connectivity)

Construct a new element where `topology` is the topological type of the element
and connectivity contains node numbers where element is connected.

# Topological types

## 1d elements
- `Seg2`
- `Seg3`

## 2d elements
- `Tri3`
- `Tri6`
- `Tri7`
- `Quad4`
- `Quad8`
- `Quad9`

## 3d elements
- `Tet4`
- `Tet10`
- `Hex8`
- `Hex20`
- `Hex27`
- `Pyr5`
- `Wedge6`
- `Wedge15`

# Examples

```julia
element = Element(Tri3, (1, 2, 3))
```
"""
function Element(::Type{T}, connectivity::NTuple{N, Int}) where {N, T<:FEMBasis.AbstractBasis}
    element_id = -1
    topology = T()
    integration_points = []
    fields = Dict()
    element = Element{T}(element_id, collect(connectivity), integration_points, fields, topology)
    return element
end

function Element(::Type{T}, connectivity::Vector{Int}) where T<:FEMBasis.AbstractBasis
    return Element(T, (connectivity...,))
end

function get_element_type(::AbstractElement{E}) where E
    return E
end

function get_element_id(element::AbstractElement{E}) where E
    return element.id
end

function is_element_type(::AbstractElement{E}, element_type) where E
    return E === element_type
end

function filter_by_element_type(element_type, elements)
    return Iterators.filter(element -> is_element_type(element, element_type), elements)
end

function get_connectivity(element::AbstractElement)
    return element.connectivity
end

"""
    group_by_element_type(elements)

Given a vector of elements, group elements by element type to several vectors.
Returns a dictionary, where key is the element type and value is a vector
containing all elements of type `element_type`.
"""
function group_by_element_type(elements)
    results = Dict{DataType, Any}()
    basis_types = map(element -> typeof(element.properties), elements)
    for basis in unique(basis_types)
        element_type = Element{basis}
        subset = filter(element -> isa(element, element_type), elements)
        results[element_type] = convert(Vector{element_type}, subset)
    end
    return results
end

### dfields - dynamically defined fields

# This is the "old" field system, where fields are defined to dictionary.
# It is known that this approach is having a performance issue caused by
# type instability.

function has_dfield(element::AbstractElement, field_name)
    return haskey(element.dfields, field_name)
end

function get_dfield(element::AbstractElement, field_name)
    return getindex(element.dfields, field_name)
end

function create_dfield!(element::AbstractElement, field_name, field_::AbstractField)
    T = typeof(field_)
    if has_dfield(element, field_name)
        @debug("Replacing the content of a field $field_name with a new field of type $T.")
    else
        @debug("Creating a new dfield $field_name of type $T")
    end
    element.dfields[field_name] = field_
    return
end

function create_dfield!(element::AbstractElement, field_name, field_data)
    create_dfield!(element, field_name, field(field_data))
end

function update_dfield!(element::AbstractElement, field_name, field_data)
    if has_dfield(element, field_name)
        field = get_dfield(element, field_name)
        @debug("Update $field_name with data $field_data")
        update_field!(field, field_data)
    else
        create_dfield!(element, field_name, field_data)
    end
end

# A helper function to pick element data from dictionary
function pick_data_(element, field_data)
    connectivity = get_connectivity(element)
    N = length(connectivity)
    picked_data = ntuple(i -> getindex(field_data, connectivity[i]), N)
    return picked_data
end

function update_dfield!(element::AbstractElement, field_name, (time, field_data)::Pair{Float64, Dict{Int,V}}) where V
    update_dfield!(element, field_name, time => pick_data_(element, field_data))
end

function update_dfield!(element::AbstractElement, field_name, field_data::Dict{Int,V}) where V
    update_dfield!(element, field_name, pick_data_(element, field_data))
end

function update_dfield!(element::AbstractElement, field_name, field_data::Function)
    if hasmethod(field_data, Tuple{Element, Any, Any})
        element.dfields[field_name] = field((ip, time) -> field_data(element, ip, time))
    else
        element.dfields[field_name] = field(field_data)
    end
end

function interpolate_dfield(element::AbstractElement, field_name, time)
    field = get_dfield(element, field_name)
    return interpolate(field, time)
end

### dfield & sfield -- common routines

function has_field(element, field_name)
    return has_dfield(element, field_name)
end

function get_field(element, field_name)
    return get_dfield(element, field_name)
end

function update_field!(element::AbstractElement, field_name, field_data)
    update_dfield!(element, field_name, field_data)
end

function interpolate_field(element::AbstractElement, field_name, time)
    interpolate_dfield(element, field_name, time)
end

function interpolate(element::AbstractElement, field_name, time)
    return interpolate_field(element, field_name, time)
end

#function update_field!(element::AbstractElement, field_name, field_data::Dict)
#    update_dfield!(element, field_name, field_data)
#end

function update_field!(elements::Vector{Element}, field_name, field_data)
    for element in elements
        update_field!(element, field_name, field_data)
    end
end

# Update fields when given a dictionary or time => dictionary:
# pick data from dictionary diven by the connectivity information of element
#=
function update_field!(element::AbstractElement, field::F,
                       data::Dict{T,V}) where {F<:DVTI,T,V}
    connectivity = get_connectivity(element)
    N = length(connectivity)
    picked_data = ntuple(i -> data[connectivity[i]], N)
    update_field!(field, picked_data)
end

function update_field!(element::AbstractElement, field::F,
                       ddata::Pair{Float64, Dict{T,V}}) where {F<:DVTV,T,V}
    time, data = ddata
    connectivity = get_connectivity(element)
    N = length(connectivity)
    picked_data = ntuple(i -> data[connectivity[i]], N)
    update_field!(field, time => picked_data)
end
=#

"""
    interpolate(element, field_name, time)

Interpolate field `field_name` from element at given `time`.

# Example
```
element = Element(Seg2, [1, 2])
data1 = Dict(1 => 1.0, 2 => 2.0)
data2 = Dict(1 => 2.0, 2 => 3.0)
update!(element, "my field", 0.0 => data1)
update!(element, "my field", 1.0 => data2)
interpolate(element, "my field", 0.5)

# output

(1.5, 2.5)

```
"""
function interpolate(element::AbstractElement, field_name::String, time::Float64)
    field = element[field_name]
    result = interpolate(field, time)
    if isa(result, Dict)
        connectivity = get_connectivity(element)
        return tuple((result[i] for i in connectivity)...)
    else
        return result
    end
end

function info_update_field(elements, field_name, data)
    nelements = length(elements)
    @info("Updating field `$field_name` for $nelements elements.")
end

function info_update_field(elements, field_name, data::Float64)
    nelements = length(elements)
    @info("Updating field `$field_name` => $data for $nelements elements.")
end

"""
    update!(elements, field_name, data)

Given a list of elements, field name and data, update field to elements. Data
is passed directly to the `field`-function.

# Examples

Create two elements with topology `Seg2`, one is connecting to nodes (1, 2) and
the other is connecting to (2, 3). Some examples of updating fields:

```julia
elements = [Element(Seg2, [1, 2]), Element(Seg2, [2, 3])]
X = Dict(1 => 0.0, 2 => 1.0, 3 => 2.0)
u = Dict(1 => 0.0, 2 => 0.0, 3 => 0.0)
update!(elements, "geometry", X)
update!(elements, "displacement", 0.0 => u)
update!(elements, "youngs modulus", 210.0e9)
update!(elements, "time-dependent force", 0.0 => 0.0)
update!(elements, "time-dependent force", 1.0 => 100.0)
```

When using dictionaries in definition of fields, key of dictionary corresponds
to node id, that is, updating field `geometry` in the example above is updating
values `(0.0, 1.0)` for the first elements and values `(1.0, 2.0)` to the second
element. For time dependent field, syntax `time => data` is used. If field is
initialized without time-dependency, it cannot be changed to be time-dependent
afterwards. If unsure, it's better to initialize field with time dependency.

"""
function update!(elements, field_name, data)
    info_update_field(elements, field_name, data)
    for element in elements
        update!(element, field_name, data)
    end
end


## Interpolate fields in spatial direction

const ConstantField = Union{DCTI, DCTV}
const VariableFields = Union{DVTV, DVTI}
const DictionaryFields = Union{DVTVd, DVTId}

function interpolate_field(::AbstractElement, field::ConstantField, ip, time)
    return interpolate_field(field, time)
end

function interpolate_field(element::AbstractElement, field::VariableFields, ip, time)
    data = interpolate_field(field, time)
    basis = get_basis(element, ip, time)
    N = length(basis)
    return sum(data[i]*basis[i] for i=1:N)
end

function interpolate_field(element::AbstractElement, field::DictionaryFields, ip, time)
    data = interpolate_field(field, time)
    basis = element(ip, time)
    N = length(element)
    c = get_connectivity(element)
    return sum(data[c[i]]*basis[i] for i=1:N)
end

function interpolate_field(::AbstractElement, field::CVTV, ip, time)
    return field(ip, time)
end

function interpolate(element::AbstractElement, field_name, ip, time)
    field = get_field(element, field_name)
    interpolate_field(element, field, ip, time)
end


## Other stuff

function get_basis(element::AbstractElement{B}, ip, ::Any) where B
    T = typeof(first(ip))
    N = zeros(T, 1, length(element))
    FEMBasis.eval_basis!(B, N, tuple(ip...))
    return N
end

function get_dbasis(element::AbstractElement{B}, ip, ::Any) where B
    T = typeof(first(ip))
    dN = zeros(T, size(element)...)
    FEMBasis.eval_dbasis!(B, dN, tuple(ip...))
    return dN
end

function (element::Element)(ip, time::Float64=0.0)
    return get_basis(element, ip, time)
end

#"""
#Examples
#julia> el = Element(Quad4, [1, 2, 3, 4]);
#julia> el([0.0, 0.0], 0.0, 1)
#1x4 Array{Float64,2}:
# 0.25  0.25  0.25  0.25
#julia> el([0.0, 0.0], 0.0, 2)
#2x8 Array{Float64,2}:
# 0.25  0.0   0.25  0.0   0.25  0.0   0.25  0.0
# 0.0   0.25  0.0   0.25  0.0   0.25  0.0   0.25
#"""
function (element::Element)(ip, time::Float64, dim::Int)
    dim == 1 && return get_basis(element, ip, time)
    Ni = vec(get_basis(element, ip, time))
    N = zeros(dim, length(element)*dim)
    for i=1:dim
        N[i,i:dim:end] += Ni
    end
    return N
end

function (element::Element)(ip, time, ::Type{Val{:Jacobian}})
    X = element("geometry", time)
    J = FEMBasis.jacobian(element.properties, X, ip)
    return J
end

function (element::Element)(ip, time::Float64, ::Type{Val{:detJ}})
    J = element(ip, time, Val{:Jacobian})
    n, m = size(J)
    if n == m  # volume element
        return det(J)
    end
    JT = transpose(J)
    if size(JT, 2) == 1  # boundary of 2d problem, || ∂X/∂ξ ||
        return norm(JT)
    else # manifold on 3d problem, || ∂X/∂ξ₁ × ∂X/∂ξ₂ ||
        return norm(cross(JT[:,1], JT[:,2]))
    end
end

function (element::Element)(ip, time::Float64, ::Type{Val{:Grad}})
    J = element(ip, time, Val{:Jacobian})
    return inv(J)*get_dbasis(element, ip, time)
end

function (element::Element)(field_name::String, ip, time::Float64, ::Type{Val{:Grad}})
    X = element("geometry", time)
    u = element(field_name, time)
    return grad(element.properties, u, X, ip)
end


function get_integration_points(element::AbstractElement{E}) where E
    # first time initialize default integration points
    if length(element.integration_points) == 0
        ips = get_integration_points(element.properties)
        element.integration_points = [IP(i, w, xi) for (i, (w, xi)) in enumerate(ips)]
    end
    return element.integration_points
end

""" This is a special case, temporarily change order
of integration scheme mainly for mass matrix.
"""
function get_integration_points(element::AbstractElement{E}, change_order::Int) where E
    ips = get_integration_points(element.properties, Val{change_order})
    return [IP(i, w, xi) for (i, (w, xi)) in enumerate(ips)]
end

""" Find inverse isoparametric mapping of element. """
function get_local_coordinates(element::AbstractElement, X::Vector, time::Float64; max_iterations=10, tolerance=1.0e-6)
    haskey(element, "geometry") || error("element geometry not defined, cannot calculate inverse isoparametric mapping")
    dim = size(element, 1)
    dim == length(X) || error("manifolds not supported.")
    xi = zeros(dim)
    dX = element("geometry", xi, time) - X
    for i=1:max_iterations
        J = element(xi, time, Val{:Jacobian})'
        xi -= J \ dX
        dX = element("geometry", xi, time) - X
        norm(dX) < tolerance && return xi
    end
    debug("get_local_coordinates", X, dX, xi)
    error("Unable to find inverse isoparametric mapping for element $element for X = $X")
end

""" Test is X inside element. """
function inside(element::AbstractElement{E}, X, time) where E
    xi = get_local_coordinates(element, X, time)
    return inside(E, xi)
end

## Convenience functions

# element("displacement", 0.0)
function (element::Element)(field_name::String, time::Float64)
    return interpolate(element, field_name, time)
end

# element("displacement", (0.0, 0.0), 0.0)
function (element::Element)(field_name::String, ip, time::Float64)
    return interpolate(element, field_name, ip, time)
end

function element_info!(bi::FEMBasis.BasisInfo{E,T}, element::AbstractElement{E}, ip, time) where {E,T}
    X = interpolate(element, "geometry", time)
    eval_basis!(bi, X, ip)
    return bi.J, bi.detJ, bi.N, bi.grad
end
