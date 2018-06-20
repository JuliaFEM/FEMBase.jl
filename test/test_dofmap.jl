# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

# # Mapping nodal degrees of freedom to matrix rows

using FEMBase
using FEMBase.Test
using FEMBase: DOFMap
import FEMBase: set_gdofs!, get_gdofs, get_gdofs!

# Initialize a new `DOFMap`:

dm = DOFMap()

@test isempty(dm.map)
@test isempty(dm.local_dof_indices)

# ## Basic usage

# First we assign that nodes 1-4, each node having dof ``u_1`` and ``u_2``, are
# connected to matrix rows 1-8:

set_gdofs!(dm, (1, 2, 3, 4), (:u1, :u2), (1, 2, 3, 4, 5, 6, 7, 8))

# Accessing of data is done using `get_gdofs`

gdofs = zeros(Int64, 4)
@test get_gdofs!(gdofs, dm, (1, 2, 3, 4), (:u1,)) == [1, 3, 5, 7]
@test get_gdofs!(gdofs, dm, (1, 2), (:u1, :u2)) == [1, 2, 3, 4]

# Order does matter:

@test get_gdofs!(gdofs, dm, (2, 1), (:u1, :u2)) == [3, 4, 1, 2]
@test get_gdofs!(gdofs, dm, (1, 2), (:u2, :u1)) == [2, 1, 4, 3]

# For fast access, we have functions to return dofs, given maximum dimension
# and dof ids. For exapmle, to get dofs u1 and u2 (in indices 1,2) for nodes
# 2,3, and maximum dimension is 3, is done using the following command. Here
# dofs are (u11, u12, u13, u21, u22, u23, u31, u32, u33), so requesting dofs
# 1,2 for nodes 2,3 is returning 4, 5, 7, 8.

@test get_gdofs!(gdofs, (2,3), 3, (1,2)) == [4, 5, 7, 8]
