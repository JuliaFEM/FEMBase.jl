# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase, Test, Tensors

element = Element(Seg2, (1, 2))
data = Vec(1,2,3)
update!(element, 0.0 => data)

