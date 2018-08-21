# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBasis

struct Poi1 <: AbstractBasis end

function get_basis(::Element{Poi1}, ::Any, ::Any)
    return [1]
end

function get_dbasis(::Element{Poi1}, ::Any, ::Any)
    return [0]
end

function (::Element{Poi1})(::Any, ::Float64, ::Type{Val{:detJ}})
    return 1.0
end

function get_integration_order(::Poi1)
    return 1
end

function get_integration_points(::Poi1, ::Int64)
    return [ (1.0, (0.0, )) ]
end

function size(::Type{Poi1})
    return (0, 1)
end

function length(::Type{Poi1})
    return 1
end

function FEMBasis.get_reference_element_coordinates(::Type{Poi1})
    Vector{Float64}[[0.0]]
end

function inside(::Union{Type{Seg2}, Type{Seg3}, Type{Quad4}, Type{Quad8},
                        Type{Quad9}, Type{Pyr5}, Type{Hex8}, Type{Hex20},
                        Type{Hex27}}, xi)
    return all(-1.0 .<= xi .<= 1.0)
end

function inside(::Union{Type{Tri3}, Type{Tri6}, Type{Tri7}, Type{Tet4}, Type{Tet10}}, xi)
    return all(xi .>= 0.0) && (sum(xi) <= 1.0)
end

function get_reference_coordinates(::Element{B}) where B
    return get_reference_element_coordinates(B)
end

