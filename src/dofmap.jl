# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

"""
    DOFMap

This type makes conversions between node id/dof pair and global dofs of freedom.

In the simplest case, each node has same number of degrees of freedom and nodes
starts from 1. For example, given continuum mechanics model having three
displacement degrees of freedom per dof, we have (u_11, u_12, u_13, u_21, u_22,
u_23, ..., u_N1, u_N2, u_N3), where in general u_ij maps to row 3*(i-1)+j in
matrix assembly level. `DOFMap` also supports variable number of degrees of
freedom per node, so some nodes can have, for instance, temperature dofs or
rotation dofs.

`map` is a dictionary, containing node id as key value, and its value is another
dictionary, containing name of the degree of freedom and its location in matrix
level. Name can be e.g. `:u1`, `:u2`, `:T` or similar get from `Problem` definition.
Then, `dofmap.map[j][:u1]` is the global index of displacement of node ``j`` in
direction of ``u_1``.

`local_dof_indices` is another dictionary used when number of degrees of freedom
is assumed to be constant. Typical setting could be then e.g. `Dict(:u1=>1, :u2=>2)`.
"""
type DOFMap
    map :: Dict{Int64, Dict{Symbol, Int64}}
    local_dof_indices :: Dict{Symbol, Int64}
end

function DOFMap()
    return DOFMap(Dict(), Dict())
end

"""
    set_gdofs!(dm, nodes, dof_names, gdofs)

Create mapping of `nodes` and `dof_names` to `gdofs`.

# Example

``jltest
dm = DOFMap()
set_gdofs!(dm, (1, 3), (:u1, :u2), (1, 2, 5, 6))
``
"""
function set_gdofs!(dm::DOFMap, nodes, dof_names, gdofs)
    k = 0
    for j in nodes
        if !haskey(dm.map, j)
            dm.map[j] = Dict()
        end
        for n in dof_names
            k = k + 1
            dm.map[j][n] = gdofs[k]
        end
    end
    return nothing
end

function get_gdofs(dofmap::DOFMap, nodes, dof_names)
    ldi = dofmap.local_dof_indices
    max_dim = length(ldi)
    if max_dim != 0
        return (max_dim*(j-1)+ldi[n] for j in nodes for n in dof_names)
    else
        return (dofmap.map[j][n] for j in nodes for n in dof_names)
    end
end

function get_gdofs!(gdofs, dofmap::DOFMap, nodes, dof_names)
    for (k,j) in enumerate(get_gdofs(dofmap, nodes, dof_names))
        gdofs[k] = j
    end
    return gdofs
end

function get_gdofs(nodes, max_dim, dof_indices)
    return (max_dim*(j-1)+k for j in nodes for k in dof_indices)
end

function get_gdofs!(gdofs, nodes, max_dim, dof_indices)
    for (k,j) in enumerate(get_gdofs(nodes, max_dim, dof_indices))
        gdofs[k] = j
    end
    return gdofs
end
