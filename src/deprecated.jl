# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

"""
    length(element)

Return the length of basis (number of nodes).
"""
function length(element::AbstractElement)
    return length(element.properties)
end

"""
    size(element)

Return the size of basis (dim, nnodes).
"""
function size(element::AbstractElement)
    return size(element.properties)
end

function getindex(element::AbstractElement, field_name::String)
    return element.fields[field_name]
end

function setindex!(element::AbstractElement, data::T, field_name) where T<:AbstractField
    element.fields[field_name] = data
end

function setindex!(element::AbstractElement, data::Function, field_name)
    if hasmethod(data, Tuple{AbstractElement, Vector, Float64})
        # create enclosure to pass element as argument
        element.fields[field_name] = field((ip,time) -> data(element,ip,time))
    else
        element.fields[field_name] = field(data)
    end
end

function setindex!(element::AbstractElement, data, field_name)
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

function size(element::AbstractElement, dim)
    return size(element)[dim]
end

""" Check existence of field. """
function haskey(element::AbstractElement, field_name)
    return haskey(element.fields, field_name)
end

# will be deprecated
function assemble!(::Assembly, ::Problem{P}, ::AbstractElement, ::Any) where P
    @warn("One must define assemble! function for problem of type $P. " *
          "Not doing anything.")
    return nothing
end

# will be deprecated
function assemble!(assembly::Assembly, problem::Problem{P},
                   elements::Vector{Element}, time) where P
    @warn("This is default assemble! function. Decreased performance can be " *
          "expected without preallocation of memory. One should implement " *
          "`assemble_elements!(problem, assembly, elements, time)` function.")
    for element in elements
        assemble!(assembly, problem, element, time)
    end
    return nothing
end

# generally a bad idea to have functions like update! not explicitly giving targe?
function update!(field::AbstractField, data)
    update_field!(field, data)
end
