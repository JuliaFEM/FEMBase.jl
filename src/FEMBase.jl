# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

module FEMBase

import Base: getindex, setindex!, convert, length, size, isapprox,
             similar, start, first, next, done, last, endof, vec,
             ==, +, -, *, /, haskey, copy, push!, isempty, empty!,
             append!, sparse, full, read

using TimerOutputs
export @timeit, print_timer

using Logging
Logging.configure(level=INFO)
if haskey(ENV, "JULIAFEM_LOGLEVEL")
    Logging.configure(level=LogLevel(ENV["JULIAFEM_LOGLEVEL"]))
end
export info, debug

using FEMBasis
using FEMBasis: AbstractBasis
using FEMQuad: get_quadrature_points
include("fields.jl")
export interpolate, update!, DCTI, DCTV, DVTI, DVTV, CVTV, DVTId, DVTVd, field
include("types.jl")
include("sparse.jl")
include("elements.jl")
include("elements_lagrange.jl")
include("integrate.jl")
include("problems.jl")
include("assembly.jl")

export FieldProblem, BoundaryProblem, Problem, Element, Assembly
export Poi1, Seg2, Seg3, Tri3, Tri6, Tri7, Quad4, Quad8, Quad9,
       Tet4, Tet10, Pyr5, Wedge6, Wedge15, Hex8, Hex20, Hex27
export update!, add_elements!, add!, get_gdofs, group_by_element_type,
       get_unknown_field_name, get_unknown_field_dimension,
       get_integration_points, initialize!, assemble!
export DCTI, DVTI, DCTV, DVTV, CCTI, CVTI, CCTV, CVTV

end
