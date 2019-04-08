# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

struct Poi1 <: FEMBasis.AbstractBasis end

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

function get_integration_points(::Poi1, ::Int)
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

function inside(::Union{Type{FEMBasis.Seg2}, Type{FEMBasis.Seg3}, Type{FEMBasis.Quad4},
                        Type{FEMBasis.Quad8}, Type{FEMBasis.Quad9}, Type{FEMBasis.Pyr5},
                        Type{FEMBasis.Hex8}, Type{FEMBasis.Hex20},
                        Type{FEMBasis.Hex27}}, xi)
    return all(-1.0 .<= xi .<= 1.0)
end

function inside(::Union{Type{FEMBasis.Tri3}, Type{FEMBasis.Tri6}, Type{FEMBasis.Tri7},
                        Type{FEMBasis.Tet4}, Type{FEMBasis.Tet10}}, xi)
    return all(xi .>= 0.0) && (sum(xi) <= 1.0)
end

function get_reference_coordinates(::Element{B}) where B
    return get_reference_element_coordinates(B)
end

