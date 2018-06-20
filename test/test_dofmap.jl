# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using FEMBase.Test
using FEMBase: DOFMap

# Initialize a new `DOFMap`:

dm = DOFMap()

@test isempty(dm.map)
@test isempty(dm.local_dof_indices)
