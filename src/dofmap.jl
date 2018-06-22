# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

abstract type AbstractDOFMap end

"""
    SimpleDOFMap

DOF maps makes conversions between node id/dof pair and global dofs. The basic
question to anwwer is that given some node `j` and its local dof index `i`, what
is the corresponding row number in global matrix assembly.

In the simplest case, each node has same number of degrees of freedom and nodes
starts from 1. For example, given continuum mechanics model having ``d=3``
displacement degrees of freedom per dof, we have ``(u_{11}, u_{12}, u_{13},
u_{21}, u_{22}, u_{23}, \ldots, u_{N1}, u_{N2}, u_{N3})``, where, in general,
``u_{ij}`` maps to row d*(i-1)+j in matrix assembly level.
"""
type SimpleDOFMap <: AbstractDOFMap
    local_dof_indices :: Dict{Symbol, Int64}
end

function SimpleDOFMap()
    return SimpleDOFMap(Dict())
end

"""
    set_local_dof_indices!(dofmap, local_dof_indices)

Set local dof indices for dofmap, connecting together dof abbreviation and its
order in node.

# Example

```julia
local_dof_indices = Dict(1 => :u1, 2 => :u2, 3 => :ur3, 4 => :T)
set_local_dof_indices!(dofmap, local_dof_indices)

# output

nothing
```
"""
function set_local_dof_indices!(dofmap, local_dof_indices)
    dofmap.local_dof_indices = local_dof_indices
    return nothing
end

"""
    get_gdofs(dofmap, nodes, dof_names)

Return global dofs for some nodes and dofs. With `SimpleDOFMap`, this is a
trivial calculation operation `gdof(i,j) = d(j-1)+i`, where ``j`` is node id,
``i`` is the local dof index and ``d`` is number of dofs per node.
"""
function get_gdofs(dofmap::SimpleDOFMap, nodes, dof_names)
    ldi = dofmap.local_dof_indices
    max_dim = length(ldi)
    return (max_dim*(j-1)+ldi[n] for j in nodes for n in dof_names)
end

function get_gdofs!(gdofs, dofmap::SimpleDOFMap, nodes, dof_names)
    for (k,j) in enumerate(get_gdofs(dofmap, nodes, dof_names))
        gdofs[k] = j
    end
    return gdofs
end

function DOFMap()
    return SimpleDOFMap()
end
