
# Developing new input/output interfaces

E.g. mesh reader or results writer.


```julia
mutable struct MeshReader <: AbstractMeshReader
end

"""
    read_mesh(mesh::MeshReader, data)

Reads mesh from disk/memory/cloud/sql/etc. and returns Mesh object
"""
function read_mesh(mesh::MeshReader, filename::String)
    # do something ...
end
```


```julia
mutable struct ResultsWriter <: AbstractResultsWriter
end

"""
    write_results!(results::ResultsWriter, data)

Given data, write calculation results back to disk/memory/cloud/sql/etc.
"""
function write_results!(results::ResultsWriter, data)
    # write results ...
end
```

# Developing new physical models

Starting point is a weak formulation, which is then discretized to elements.


```julia
mutable struct Heat <: AbstractProblem
end

"""
    assemble!(settings::Heat, problem::Problem, assembly::Assembly, elements::Vector{Element}, time::Float64)

Given `problem` spesification, assemble `elements` to global stiffness matrix and force
vector defined in `assembly` for some given `time`.
"""
function assemble!(settings::Heat, problem::Problem, assembly::Assembly,
                   elements::Vector{Element}, time::Float64)
    # integrate local matrices and add them to assembly
end
```

# Developing new material model

Aim is to define material response given data.


```julia
mutable struct LinearIsotropic <: AbstractMaterial
end

"""
    calculate_material_response!(material::LinearIsotropic, data)

Given strain tensor and some other quantities, calculate material response
"""
function calculate_response!(material::LinearIsotropic, data)
    # given strain tensor and some other quantities, calculate stress
end
```

# Developing new interpolation functions

E.g. basis functions


```julia
type LinQuad4Basis <: AbstractBasis
end

"""
    evaluate_basis!(basis::LinQuad4Basis, element::Element, xi, time, N::Matrix{Float64})

Evaluate basis functions at some point $\xi$ and store results to `N`.
"""
function evaluate_basis!(basis::LinQuad4Basis, element::Element, xi, time, N::Matrix{Float64})
    # populate N with new basis
end
```

# Developing new integration rules


```julia
mutable struct Quad4PointGaussLegendre <: AbstractIntegrationRule
end

"""
    get_integration_points(q::Quad4PointQaussLegendre)

Return integration point locations and weights.
"""
function get_integration_points(q::Quad4PointGaussLegendre)
    # return integration points
end
```

# Developing new solver


```julia
mutable struct ImplicitTimeSolver <: AbstractSolver
end

"""
    solve!(settings::ImplicitTimeSolver, solver::Solver, time)

Assemble problems, solve problem, update problems, write results and so on.
"""
function solve!(settings::ImplicitTimeSolver, solver::Solver, time)
    # do solution
end
```
