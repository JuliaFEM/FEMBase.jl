# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using Test

## Elements - field sets

# Each element have two kind of fields, dynamically defined and statically
# defined.

element = Element(Quad4, (1, 2, 3, 4))

### Dynamically defined fields

using FEMBase: has_dfield, get_dfield, create_dfield!, update_dfield!,
               interpolate_dfield

# Like name suggests, dynamically defined fields can be defined on fly.
# For example, to define a fields called *geometry* and *displacement*:

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


### Common routines to dfields and sfields

using FEMBase: has_field, get_field, update_field!, interpolate_field

# There is common routines `has_field`, `update_field!`, `interpolate_field`
# and so on, which tries to be clever. We prefer statically defined fields
# as they outperform dynamically defined ones in access time, but then e.g.
# if field is not statically defined, `update_field!` is creating a new dfield.

@test has_field(element, :geometry) == true
@test get_field(element, :youngs_modulus) === get_dfield(element, :youngs_modulus)
update_field!(element, :poissons_ratio, 1/3)
@test interpolate_field(element, :poissons_ratio, 0.0) == 1/3

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
