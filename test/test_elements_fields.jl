# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Test

## Elements - field sets

# Each element have two kind of fields, dynamically defined and statically
# defined.

### Dynamically defined fields (dfields)

using FEMBase: has_dfield, get_dfield, create_dfield!, update_dfield!,
               interpolate_dfield

# Like name suggests, dynamically defined fields can be defined on fly.
# For example, to define a fields called *geometry* and *displacement*:

element = Element(Quad4, (1, 2, 3, 4))

@test has_dfield(element, :geometry) == false

X = ([0.0,0.0], [1.0,0.0], [1.0,1.0], [0.0,1.0])
u0 = ([0.0,0.0], [0.0,0.0], [0.0,0.0], [0.0,0.0])
u1 = ([0.0,0.0], [0.0,0.0], [1.0,0.0], [0.0,0.0])

create_dfield!(element, :geometry, X)
@test has_dfield(element, :geometry) == true

create_dfield!(element, :displacement, 0.0 => u0)
update_dfield!(element, :displacement, 1.0 => u1)

# The main property of the dynamic fields are that they can be defined
# abolutely in what state of simulation ever, without any predefined
# initialization steps. `update_dfield!` is creating field on fly, if
# it doesn't exists before hand.

update_dfield!(element, :youngs_modulus, 288.0)
@test has_dfield(element, :youngs_modulus)

# If field depends from time (DVTV, DCTV, DVTVd, ...), interpolation
# works as expected:

um = interpolate_dfield(element, :displacement, 0.5)
@test um[3] == [0.5, 0.0]

### Statically defined fields (sfields)

# The problem of the dynamically defined fields are that they have a poor
# performance. Julia is JIT compiled language and from performance perspective,
# knowing concrete types in compile time leads to fast implementation. For that
# reason, we also have another fieldset, which must be given when a new element
# is created, allowing high performance.

using FEMBase: has_sfield, get_sfield, update_sfield!, interpolate_sfield, AbstractFieldSet

struct MyFieldSet{N} <: AbstractFieldSet{N}
    geometry :: DVTI{N, Vector{Float64}}
    displacement :: DVTV{N, Vector{Float64}}
end

function MyFieldSet{N}() where N
    # create tuple of length N with zeros(3): ([0.0,0.0,0.0], ..., )
    data = ntuple(i -> zeros(3), N)
    geometry = DVTI(data)
    displacement = DVTV(0.0 => data)
    return MyFieldSet{N}(geometry, displacement)
end

# When constructing new element, one must pass now the field set:

element = Element(Quad4, MyFieldSet, (1, 2, 3, 4))

@test has_sfield(element, :geometry) == true

X = ([0.0,0.0], [1.0,0.0], [1.0,1.0], [0.0,1.0])
u0 = ([0.0,0.0], [0.0,0.0], [0.0,0.0], [0.0,0.0])
u1 = ([0.0,0.0], [0.0,0.0], [1.0,0.0], [0.0,0.0])

# Indeed using sfield does not mean that they cannot be changed. The only
# difference is that they cannot be created on the fly, thus there is no
# `create_dfield!` at all. But `update_dfield!` works as expected.

@test has_sfield(element, :geometry) == true

update_sfield!(element, :geometry, X)
update_sfield!(element, :displacement, 0.0 => u0)
update_sfield!(element, :displacement, 1.0 => u1)

# Interpolating fields etc. works as expected:

um = interpolate_sfield(element, :displacement, 0.5)
@test um[3] == [0.5, 0.0]

### Common routines to dfields and sfields


using FEMBase: has_field, get_field, update_field!, interpolate_field

# There is common routines `has_field`, `update_field!`, `interpolate_field`
# and so on, which tries to be clever. We prefer statically defined fields
# as they outperform dynamically defined ones in access time, but then e.g.
# if field is not statically defined, `update_field!` is creating a new dfield.

element = Element(Quad4, MyFieldSet, (1, 2, 3, 4))

update_field!(element, :geometry, X)
update_field!(element, :displacement, 0.0 => u0)
update_field!(element, :youngs_modulus, 288.0)
update_field!(element, :poissons_ratio, 1/3)

# For example, geometry and displacement are here sfields because they are
# predefined using fieldset, but youngs modulus is dfield and created on-the-fly:

@test has_field(element, :geometry) == true
@test has_dfield(element, :geometry) == false
@test has_sfield(element, :geometry) == true
@test has_dfield(element, :youngs_modulus) == true
@test has_sfield(element, :youngs_modulus) == false
@test get_field(element, :youngs_modulus) === get_dfield(element, :youngs_modulus)
@test interpolate_field(element, :poissons_ratio, 0.0) == 1/3
@test interpolate_field(element, :geometry, 0.0) == X

### Creating and updating DVTI and DVTV fields using dictionaries

# It's a very convenient way to update fields using dictionaries.

element = Element(Seg2, (1, 2))
update!(element, :temperature, 0.0 => (0.0, 0.0))

# Let say that one wants to update the following temperature data:

T = Dict(1 => 1.0, 2 => 2.0, 3 => 3.0)

# to a element, can be done using that dictionary:

update!(element, :temperature, 1.0 => T)

# Internally, we pick the data from dictionary based on the node numbers and
# update only that data to element.

@test interpolate(element, :temperature, 1.0) == (T[1], T[2])

### Performance difference between fields

# The reason, like explained, for having two kinds of fields comes from the
# performance. We need certain fields, like geometry, displacement, connectivity
# and so on, work with maximum performance, because they are accessed all the
# time. To make a short timing, we create geometry and geometry2, where the
# first one is statically defined and latter one dynamically.

#=
X = ([0.0,0.0], [1.0,0.0], [1.0,1.0], [0.0,1.0])
element = Element(Quad4, MyFieldSet, (1, 2, 3, 4))
update!(element, :geometry, X)
update!(element, :geometry2, X)
using BenchmarkTools
@btime interpolate($element, :geometry, 0.0)
@btime interpolate($element, :geometry2, 0.0)
=#

# The following code will give in first case 1.466 ns and 0 memory allocations,
# where the latter case is 37.755 ns with 1 memory allocation:
#
#  1.466 ns (0 allocations: 0 bytes)
#  37.755 ns (1 allocation: 16 bytes)
