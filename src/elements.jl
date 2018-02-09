# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

type Element{E<:AbstractBasis}
    id :: Int
    connectivity :: Vector{Int}
    integration_points :: Vector{IP}
    fields :: Dict{String, AbstractField}
    properties :: E
end

"""
    Element(element_type, connectivity_vector)

Construct a new element where element_type is the type of the element
and connectivity_vector is the vector of nodes that the element is connected to.

Examples
--------
In the example a new element (E in the figure below) of type Tri3 is created.
This spesific element connects to nodes 89, 43, 12 in the finite element mesh.

```@example
element = Element(Tri3, [89, 43, 12])
```
![img](figs/mesh.png)
"""
function Element{E<:AbstractBasis}(::Type{E}, connectivity::Vector{Int})
    return Element{E}(-1, connectivity, [], Dict(), E())
end

"""
    length(element)

Return the length of basis (number of nodes).
"""
function length(element::Element)
    return length(element.properties)
end

"""
    size(element)

Return the size of basis (dim, nnodes).
"""
function size(element::Element)
    return size(element.properties)
end

function getindex(element::Element, field_name::String)
    return element.fields[field_name]
end

function setindex!{T<:AbstractField}(element::Element, data::T, field_name)
    element.fields[field_name] = data
end

function get_element_type{E}(::Element{E})
    return E
end

function get_element_id{E}(element::Element{E})
    return element.id
end

function is_element_type{E}(::Element{E}, element_type)
    return E === element_type
end

function filter_by_element_type(element_type, elements)
    return Iterators.filter(element -> is_element_type(element, element_type), elements)
end

function get_connectivity(element::Element)
    return element.connectivity
end

"""
    group_by_element_type(elements::Vector{Element})

Given a vector of elements, group elements by element type to several vectors.
Returns a dictionary, where key is the element type and value is a vector
containing all elements of type `element_type`.
"""
function group_by_element_type(elements::Vector{Element})
    results = Dict{DataType, Any}()
    basis_types = map(element -> typeof(element.properties), elements)
    for basis in unique(basis_types)
        element_type = Element{basis}
        subset = filter(element -> isa(element, element_type), elements)
        results[element_type] = convert(Vector{element_type}, subset)
    end
    return results
end

function setindex!(element::Element, data::Function, field_name)
    if method_exists(data, Tuple{Element, Vector, Float64})
        # create enclosure to pass element as argument
        element.fields[field_name] = field((ip,time) -> data(element,ip,time))
    else
        element.fields[field_name] = field(data)
    end
end

function setindex!(element::Element, data, field_name)
    element.fields[field_name] = field(data)
end

#""" Return a Field object from element.
#Examples
#--------
#>>> element = Element(Seg2, [1, 2])
#>>> data = Dict(1 => 1.0, 2 => 2.0)
#>>> update!(element, "my field", data)
#>>> element("my field")
#"""
function (element::Element)(field_name::String)
    return element[field_name]
end

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
function interpolate(element::Element, field_name::String, time::Float64)
    field = element[field_name]
    result = interpolate(field, time)
    if isa(result, Dict)
        connectivity = get_connectivity(element)
        return tuple((result[i] for i in connectivity)...)
    else
        return result
    end
end

function get_basis{B}(element::Element{B}, ip, ::Any)
    T = typeof(first(ip))
    N = zeros(T, 1, length(element))
    eval_basis!(B, N, tuple(ip...))
    return N
end

function get_dbasis{B}(element::Element{B}, ip, ::Any)
    T = typeof(first(ip))
    dN = zeros(T, size(element)...)
    eval_dbasis!(B, dN, tuple(ip...))
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
    J = jacobian(element.properties, X, ip)
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

function interpolate(element::Element, field_name, ip, time)
    field = element[field_name]
    return interpolate_field(element, field, ip, time)
end

function interpolate_field{T<:Union{DCTI,DCTV}}(::Element, field::T, ::Any, time)
    return interpolate_field(field, time)
end

function interpolate_field{T<:Union{DVTV,DVTI}}(element::Element, field::T, ip, time)
    data = interpolate_field(field, time)
    basis = get_basis(element, ip, time)
    N = length(basis)
    return sum(data[i]*basis[i] for i=1:N)
end

function interpolate_field{T<:Union{DVTVd,DVTId}}(element::Element, field::T, ip, time)
    data = interpolate_field(field, time)
    basis = element(ip, time)
    N = length(element)
    c = get_connectivity(element)
    return sum(data[c[i]]*basis[i] for i=1:N)
end

function interpolate_field(::Element, field::CVTV, ip, time)
    return field(ip, time)
end

function size(element::Element, dim)
    return size(element)[dim]
end

function update!{F<:AbstractField}(::Element, field::F, data)
    update!(field, data)
end

function update!{F<:AbstractField,T,V}(element::Element, field::F,
                                       data_::Dict{T,V})
    data = (collect(data_[i] for i in get_connectivity(element))...)
    update!(field, data)
end

function update!{F<:AbstractField,T,V}(element::Element, field::F,
                                       data__::Pair{Float64, Dict{T,V}})
    time_, data_ = data__
    data = (collect(data_[i] for i in get_connectivity(element))...)
    update!(field, time_ => data)
end

function update!(element::Element, field_name, data)
    if haskey(element.fields, field_name)
        update!(element, element.fields[field_name], data)
    else
        element.fields[field_name] = field(data)
    end
end

function update!{F<:AbstractField}(element::Element, field_name, field_::F)
    element.fields[field_name] = field_
end

function update!(element::Element, field_name, data::Function)
    if method_exists(data, Tuple{Element, Any, Any})
        element.fields[field_name] = field((ip, time) -> data(element, ip, time))
    else
        element.fields[field_name] = field(data)
    end
end

function update!(elements::Vector, field_name, data)
    for element in elements
        update!(element, field_name, data)
    end
end

""" Check existence of field. """
function haskey(element::Element, field_name)
    return haskey(element.fields, field_name)
end

function get_integration_points{E}(element::Element{E})
    # first time initialize default integration points
    if length(element.integration_points) == 0
        ips = get_integration_points(element.properties)
        if E in (Seg2, Seg3, NSeg)
            element.integration_points = [IP(i, w, (xi,)) for (i, (w, xi)) in enumerate(ips)]
        else
            element.integration_points = [IP(i, w, xi) for (i, (w, xi)) in enumerate(ips)]
        end
    end
    return element.integration_points
end

""" This is a special case, temporarily change order
of integration scheme mainly for mass matrix.
"""
function get_integration_points{E}(element::Element{E}, change_order::Int)
    ips = get_integration_points(element.properties, Val{change_order})
    if E in (Seg2, Seg3, NSeg)
        return [IP(i, w, (xi,)) for (i, (w, xi)) in enumerate(ips)]
    else
        return [IP(i, w, xi) for (i, (w, xi)) in enumerate(ips)]
    end
end

""" Find inverse isoparametric mapping of element. """
function get_local_coordinates(element::Element, X::Vector, time::Float64; max_iterations=10, tolerance=1.0e-6)
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
    info("X = $X, dX = $dX, xi = $xi")
    error("Unable to find inverse isoparametric mapping for element $element for X = $X")
end

""" Test is X inside element. """
function inside{E}(element::Element{E}, X, time)
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

function element_info!{E,T}(bi::BasisInfo{E,T}, element::Element{E}, ip, time)
    X = interpolate(element, "geometry", time)
    eval_basis!(bi, X, ip)
    return bi.J, bi.detJ, bi.N, bi.grad
end
