# FEMBase.jl

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

## Nomenclature

In this package, the roof concept is analysis, which can mean any kind of action
done given a set of problems. We analyze a problem. Analysis can be, for
example, integrating the volume of some area, calculating the pressure,
temperature or displacement of some volume, calculating it's eigenvalues or
things like that. In even highler level, a simulation is a sequence of analyses
performed in order to achieve of final deliverable. In commercial softwares,
analysis would resemble maybe like "STEP" in ABAQUS.

The subject of analysis doesn't have to be anything what we physically observe,
like a solid volume. In abstract sense, it's just a set of problems to where we
apply a series of actions. At this level of abstraction, we don't care about the
internal structure of the problem. So what are then the problems?

Each analysis is having a set of problems. There can be only one or there can be
several. One could think "problem" something equivalent to "body" in some other
FEM codes, or something similar, byt it's not physically restricted in any
manner. Problem describes the physics, what is considered to be solved. The
reason why problems are not called to bodies in FEMBase is the following: not
all problems relates to the volume.

For example, when elastic bodies are having surface interactions, it's called to
contact problem, and the problem can be fully described just by knowing the
boundary / surface of the bodies. In the final discretized system, we do not
need any information about the internal structure of the volume. Thus it would
be unnecessary complication to carry all that data with the problem.

In problem level, we need to start thinking about how we describe it. In the
context of finite elements, we are describing the problem with the help of
finite elements. And in this level, we need to take care of our solution
strategy. Discretized problem has to be somehow transferred to something what
linear system solver understands. So, the problem is kind of "element set".
Problem have a matrix level representation and elements describing the problem.
This is the level where new physical rules are applied. Given a set of elements,
populate the matrices ready for solver.

For this, we need to have some mapping, which connect each element local degrees
of freedom to the global ones.

The element itself is an individual entity. It doesn't need anything surrounding
complexity to fully describe it's own state. It can be connected to problem, but
it's not mandatory. Element is having it's own topology, basis functions and
fields. Element doesn't even need a connection information to other elements.
Element needs to know only it's own geometry. Given basis and data, element can
interpolate its results in the geometrical domain occupied by that element.
Element can form it own local matrices. Your own one element FEM should be
possible to do just by using `Element`, without any complexities.

When doing FEM, we also need to integrate, which is done numerically using some
discrete integration scheme. Thus we need integation points. They are defined
inside elements.

    analysis -> problems -> elements -> integration points

In this description, we did go in top-down-approach. This same setup should also
make sense when going in other direction.

## Outline

This documentation aims to describe the structures of the package in "bottom to
up" manner. starting from the smaller data types proceeding to the big scope.
