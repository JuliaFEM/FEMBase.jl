# FEMBase.jl

[![Build status][travis-img]][travis-url]
[![Coverage status][coveralls-img]][coveralls-url]
[![Stable docs][docs-stable-img]][docs-stable-url]
[![Latest docs][docs-latest-img]][docs-latest-url]

![somebody](docs/src/figs/somebody.svg)

Package author: Jukka Aho (@ahojukka5, ahojukka5@gmail.com)

FEMBase.jl is a base package for FEM codes. It implements all the basic data
structures in an efficient way and aims to be the rock-solid backbone which can
be used to build your very own FEM solver. The basis principle is that FEMBase
is a tiny and well-optimized, efficient subset of larger finite element method
software with pre- and postprocessors. FEMBase.jl is used in JuliaFEM.jl.

In order to make a stable package, all operations that are not closely related
to the basic operations are left out. This encourages to develop a highly
modular FEM where features are implemented in their own packages. The stragegy
chosen is intentional. Pre- and posprocessing, physical models, file formats and
things like that evolve in time, while the basic principles stay. We want that
those evolving things stay in their own packages as much as possible.

## What FEMBase.jl does

FEMBase implements all kind of basic data types like `Problem`, `Element`,
`Node`, `IntegrationPoint`, `Analysis` and so on so that user can define the
discretized geometry and integration scheme. Each element does have a topology
like `Poi1`, `Seg2`, `Tri3`, `Quad4`, `Tet4`, `Hex8` for isoparametric linear
elements having only corner nodes or `Seg3`, `Tri6`, `Quad8`, `Tet10`, `Hex20`
for isoparametric quadratic elements having nodes in middle of the edges. Then
there is also elements having middle nodes like `Tri7` and `Quad9`.

Data of model is repsented in fields. Data type can be scalar, vector, tensor or
something user-defined for spesific needs. Data can be defined in analysis,
problem, element, node or integration point. Data can depend on time and it can
be interpolated with respect to time using different interpolation schemes.

Data can also be interpolated in spatial direction using basis functions, which
are Lagrange basis functions at the moment.

## What FEMBase.jl doesn't do

FEMBase.jl doesn't solve the system of equations. You can define model and data,
integrate and interpolate, but you need to write your own solver. And as a
consequence it doesn't write the results to anywhere, so you need to implement
your own results writer or use some of the existing solutions like
`XdmfWriter.jl`.

FEMBase.jl doesn't have any physical models. FEMBase doesn't know anything about
displacements, temperatures, velocities or pressures. FEMBase.jl does assemble a
problem, but it doesn't define it. So you need to build your physical model in
some other package or use some of the existing solutions like `HeatTransfer.jl`.
Things related to geometry, like calculating Jacobians, are supported.

FEMBase does not parse mesh from other format. For that, you need to implement
your own parser or use some existing solutions, like `AbaqusReader.jl`,
`AsterReader.jl` or `Gmsh.jl`.

## What FEMBase.jl is planning to do

Given data and basis functions, we should be able to naturally interpolate any
kind of data. In current implementation, we have an assumption that data is
defined in nodes like it is done standard FE methods. Then we have things like
hierarchical elements, NURBS, material particle method and so on, which should
also be possible to implement using FEMBase.

## Contributing

Contributions in code and idea level are welcome and highly appreciated. If you
have some great ideas how to improve package, raise an issue or drop me email.

[travis-img]: https://travis-ci.org/JuliaFEM/FEMBase.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaFEM/FEMBase.jl
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://juliafem.github.io/FEMBase.jl/stable
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://juliafem.github.io/FEMBase.jl/latest
[coveralls-img]: https://coveralls.io/repos/github/JuliaFEM/FEMBase.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JuliaFEM/FEMBase.jl?branch=master
