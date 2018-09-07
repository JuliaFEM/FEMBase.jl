# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

module FEMBase

using LinearAlgebra

macro lintpragma(s)
end

import Base: getindex, setindex!, convert, length, size, isapprox,
             first, last, vec,
             ==, +, -, *, /, haskey, copy, push!, isempty, empty!,
             append!, read

using TimerOutputs
export @timeit, print_timer

import FEMBasis
import FEMQuad

include("fields.jl")
export interpolate, update!, DCTI, DCTV, DVTI, DVTV, CVTV, DVTId, DVTVd, field
include("types.jl")
include("sparse.jl")
include("elements.jl")
export element_info!
include("elements_lagrange.jl")
include("integrate.jl")

include("problems.jl")
export set_gdofs!, get_gdofs, get_formulation_type, add_element!, add_elements!

include("assembly.jl")
export assemble_prehook!, assemble_posthook!, assemble_elements!

include("solvers.jl")
export LinearSystem, AbstractLinearSystemSolver, solve!

include("analysis.jl")
export AbstractAnalysis, Analysis, add_problems!, get_problems, run!,
       AbstractResultsWriter, add_results_writer!, get_results_writers,
       write_results!, get_problem

export FieldProblem, BoundaryProblem, Problem, Element, Assembly
export add!, group_by_element_type,
       get_unknown_field_name, get_unknown_field_dimension,
       get_integration_points, initialize!, assemble!
export SparseMatrixCOO, SparseVectorCOO, Node, BasisInfo,
       IP, AbstractProblem, IntegrationPoint, AbstractField
export is_field_problem, is_boundary_problem, get_elements,
       get_connectivity,
       get_parent_field_name, get_reference_coordinates,
       get_reference_element_coordinates,
       get_assembly, get_nonzero_rows, get_nonzero_columns,
       eval_basis!, get_basis, get_dbasis, grad!,
       assemble_mass_matrix!, get_local_coordinates, inside,
       get_element_type, filter_by_element_type, get_element_id,
       resize_sparse, resize_sparsevec

include("test.jl")

function eliminate_boundary_conditions! end
export eliminate_boundary_conditions!

function add_slave_elements! end
export add_slave_elements!

function get_slave_elements end
export get_slave_elements

function add_master_elements! end
export add_master_elements!

function get_master_elements end
export get_master_elements

using FEMBasis
export Poi1, Seg2, Seg3, Tri3, Tri6, Tri7, Quad4, Quad8, Quad9,
       Tet4, Tet10, Pyr5, Wedge6, Wedge15, Hex8, Hex20, Hex27

end
