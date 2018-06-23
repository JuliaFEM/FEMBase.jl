# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

# # Mapping nodal degrees of freedom to matrix rows (global dofs)

using FEMBase
using FEMBase.Test
using FEMBase: DOFMap, set_local_dof_indices!

# With `SimpleDOFMap`, matrix row indices, what we call to global dofs, are
# calculated based on the number of node ``j`` and maximum number of dofs per
# node ``d``. So for example, if each node has two dofs ``u_1``, ``u_2``, we
# have ``(u_{11}, u_{12}, u_{21}, u_{22}, \ldots, u_{N1}, u_{N2})``. In general,
# this formula is then `gdof(i,j) = d*(1-j)+i`. `SimpleDOFMap` mode does not
# need any storage of mapping because matrix row is explicitly determined based
# on node id number and local dof index.

dm = DOFMap()
set_local_dof_indices!(dm, Dict(:u1=>1, :u2=>2, :T=>3))

# Accessing of data is done using `get_gdofs` or in-place version `get_gdofs!`:

@test collect(get_gdofs(dm, (1, 3), (:u1, :u2))) == [1, 2, 7, 8]

# In general, `get_gdofs` takes two iterables, one describing nodes and another
# describing name of the dofs. Order does matter:

@test collect(get_gdofs(dm, (3, 1), (:u1, :u2))) == [7, 8, 1, 2]
@test collect(get_gdofs(dm, (1, 3), (:u2, :u1))) == [2, 1, 8, 7]

# There is also in-place version to get dofs:

gdofs = zeros(Int64, 4)
@test collect(get_gdofs!(dm, (1, 2), (:u1, :u2))) == [1, 2, 3, 4]
