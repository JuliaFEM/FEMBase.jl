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
    length(element::Element)

Return the number of nodes in element.
"""
function length{B}(element::Element{B})
    return length(B)
end

function size{B}(element::Element{B})
    return size(B)
end

function getindex(element::Element, field_name::AbstractString)
    return element.fields[field_name]
end

function setindex!{T<:AbstractField}(element::Element, data::T, field_name)
    element.fields[field_name] = data
end

function get_element_type{E}(element::Element{E})
    return E
end

function get_element_id{E}(element::Element{E})
    return element.id
end

function is_element_type{E}(element::Element{E}, element_type)
    return E === element_type
end

function filter_by_element_type(element_type, elements)
    return Iterators.filter(element -> is_element_type(element, element_type), elements)
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

#""" Return a Field object from element and interpolate in time direction.
#Examples
#--------
#>>> element = Element(Seg2, [1, 2])
#>>> data1 = Dict(1 => 1.0, 2 => 2.0)
#>>> data2 = Dict(1 => 2.0, 2 => 3.0)
#>>> update!(element, "my field", 0.0 => data1, 1.0 => data2)
#>>> element("my field", 0.5)
#"""
function (element::Element)(field_name::String, time::Float64)
    field = element[field_name]
    result = interpolate(field, time)
    if isa(result, Dict)
        return tuple((result[i] for i in get_connectivity(element))...)
    else
        return result
    end
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

function (element::Element)(field_name::String, ip, time::Float64)
    field = element[field_name]
    return element(field, ip, time)
end

function (element::Element)(field::DCTI, ip, time::Float64)
    return field.data
end

function (element::Element)(field::DCTV, ip, time::Float64)
    return interpolate(field, time)
end

function (element::Element)(field::CVTV, ip, time::Float64)
    return field(ip, time)
end

function (element::Element){F<:AbstractField}(field::F, ip, time::Float64)
    field_ = interpolate(field, time)
    basis = element(ip, time)
    n = length(element)
    if isa(field_, Tuple)
        m = length(field_)
        if n != m
            error("Error when trying to interpolate field $field at coords $ip and time $time: element length is $n and field length is $m, f = Nᵢfᵢ makes no sense!")
        end
        return sum(field_[i]*basis[i] for i=1:n)
    else
        c = get_connectivity(element)
        return sum(field_[c[i]]*basis[i] for i=1:n)
    end
end

function size(element::Element, dim)
    return size(element)[dim]
end

function update!(element::Element, field_name, data)
    if haskey(element.fields, field_name)
        update!(element.fields[field_name], data)
    else
        element.fields[field_name] = field(data)
    end
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

function get_connectivity(element::Element)
    return element.connectivity
end

function get_integration_points{E}(element::Element{E})
    # first time initialize default integration points
    if length(element.integration_points) == 0
        ips = get_integration_points(element.properties)
        if E in (Poi1, Seg2, Seg3, NSeg)
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
    if E in (Poi1, Seg2, Seg3, NSeg)
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
