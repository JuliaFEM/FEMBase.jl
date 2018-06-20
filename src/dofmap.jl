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
"""
type DOFMap
    map :: Dict{Int64, Dict{Symbol, Int64}}
    local_dof_indices :: Dict{Symbol, Int64}
end

function DOFMap()
    return DOFMap(Dict(), Dict())
end
