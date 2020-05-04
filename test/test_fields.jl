using FEMBase, Test

# From the beginning of a project we had a clear concept in our mind: "everything
# is a field". That is, everything can vary temporally and spatially. We think
# that constant is just a special case of field which does not vary in temporal
# nor spatial direction. Fields can vary in spatial direction, i.e. can be either
# constant or variable, and in temporal direction, i.e. can be time variant or
# time invariant. From this pondering we can think that there exists four kind of
# (discrete) fields:

# - discrete, constant, time invariant (DCTI)
# - discrete, variable, time invariant (DVTI)
# - discrete, constant, time variant (DCTV)
# - discrete, variable, time variant (DVTV)

# Discrete, in this context, means that field is defined in point-wise in
# $1 \ldots n$ locations, from where it is then interpolated to whole domain
# using some interpolation polynomials, i.e.
# ```math
# u(\xi, t) = \sum_{i} u_i[t] N_{i}(\xi,t),
# ```math
# where
#     $N_{i}(\xi, t)$
# is the basis function or interpolation polymial corresponding to $i$^{th}
# discrete value and 
#     $u_{i}$
# is the discrete value.

# Then we have continuous fields, which are defined in whole domain, or at least
# not point-wise. By following the already used abbreviations, we have four more
# fields:

# - continuous, constant, time invariant (CCTI)
# - continuous, variable, time invariant (CVTI)
# - continuous, constant, time variant (DCTV)
# - continuous, variable, time variant (CVTV)

# Continuous, again in this context, does not mean that field has to be defined
# everywhere. It's enough that it's defined in function of spatial and/or temporal
# coordinates, i.e. we have $u \equiv u(\xi, t)$, without a some spesific basis
# needed to interpolate from discrete values. 

# Field itself can be in principle anything. However, usually either scalar,
# vector or tensor (matrix). Time does not to have be real, it can be for example
# angle of some rotating machine or even complex value. 

# From these starting points, we assume that the mentioned field system can
# describe all imaginable situations.

# ## Creating new fields

# For discrete fields that are varying in spatial direction, value for each
# discrete point is defined using NTuple. The order of points is implicitly
# assumed to be same than node ordering in ABAQUS. That is, first corner nodes
# in anti-clockwise direction and after that middle nodes.

# For example, `(1, 2, 3, 4)` is a scalar field having length of 4 and
# `([1,2],[2,3],[3,4],[4,5])` is a vector field having length of 4.

# For fields that are varying in temporal direction, `time => value` syntax is
# used. The first item in pair is time (or similar) and second item is value 
# assigned to that time. For example, `0.0 => 1.0` is a time-dependent scalar
# field having value 1.0 at time 0.0.

# ## Dicrete, constant, time invariant field (DCTI)

# The most simple field is a field that is constant in both time and spatial
# direction. Discrete, constant, time invariant field. For example, youngs
# modulus could be this kind of field.

a = DCTI(1)

# Accessing data is done using `interpolate`. In FEM codes, we try to hide the
# actual type of the field, so for example interpolating constant field works,
# but the result is quite unsuprising.

@test interpolate(a, 0.0) == 1

# Field value value can be updated with `update!` function:

update!(a, 2)
@test a == 2

# Constant field of course doesn't have to be scalar field. It can be e.g.
# vector field. I use here packate Tensors.jl because of its excellent
# performance and other features, but normal `Vector` would work just fine
# also:

using Tensors

b = DCTI(Vec(1, 2))

# Interpolation, again, returns just the original data:

@test interpolate(b, 0.0) == [1, 2]

# Updating field is done using `update!`-function:

update!(b, Vec(2, 3))
@test interpolate(b, 0.0) == [2, 3]

# Constant tensor field:

c = DCTI(Tensor{2,2}((1.0, 2.0, 3.0, 4.0)))

# Data can be accessed also using `getindex`. Also things like `length` and
# `size` are defined.

@test interpolate(c, 0.0) == [1 3; 2 4]
@test c[1] == [1 3; 2 4]
@test length(c) == 4
@test size(c) == (2, 2)

# For now everything might look like extra complexity, but later on we see how
# to combine field with some basis functions in order to interpolate in element
# domain. Another nice feature is that we can interpolate fields in time. In this
# particular case of time invariant fields it of course doesn't give anything
# extra.

# ## Dicrete, variable, time invariant fields (DVTI)

@testset "DVTI field" begin

    # scalar field
    a = DVTI((1, 2))
    @test a[1] == 1
    @test a[2] == 2
    @test interpolate(a, 0.0) == (1, 2)
    update!(a, (2, 3))
    @test a == (2, 3)
    @test (2, 3) == a

    # vector field
    b = DVTI(([1, 2], [2, 3]))
    @test b[1] == [1, 2]
    @test b[2] == [2, 3]
    @test interpolate(b, 0.0) == ([1, 2], [2, 3])
    update!(b, ([2, 3], [4, 5]))
    @test b == ([2, 3], [4, 5])

    # tensor field
    c = DVTI(([1 2; 3 4], [2 3; 4 5]))
    @test c[1] == [1 2; 3 4]
    @test c[2] == [2 3; 4 5]
    @test interpolate(c, 0.0) == ([1 2; 3 4], [2 3; 4 5])
    update!(c, ([2 3; 4 5], [5 6; 7 8]))
    @test c == ([2 3; 4 5], [5 6; 7 8])

    d = DVTI(2, 3)
    @test a == d
end

@testset "DCTV field" begin

    # scalar field
    a = DCTV(0.0 => 0.0, 1.0 => 1.0)
    @test isapprox(interpolate(a, -1.0), 0.0)
    @test isapprox(interpolate(a, 0.0), 0.0)
    @test isapprox(interpolate(a, 0.5), 0.5)
    @test isapprox(interpolate(a, 1.0), 1.0)
    update!(a, 1.0 => 2.0)
    @test isapprox(interpolate(a, 0.5), 1.0)
    update!(a, 2.0 => 1.0)
    @test isapprox(interpolate(a, 1.5), 1.5)

    # vector field
    b = DCTV(0.0 => [1.0, 2.0], 1.0 => [2.0, 3.0])
    @test isapprox(interpolate(b, 0.5), [1.5, 2.5])

    # tensor field
    c = DCTV(0.0 => [1.0 2.0; 3.0 4.0], 1.0 => [2.0 3.0; 4.0 5.0])
    @test isapprox(interpolate(c, 0.5), [1.5 2.5; 3.5 4.5])
end

@testset "DVTV field" begin
    # scalar field
    a = DVTV(0.0 => (0.0, 1.0), 1.0 => (1.0, 0.0))
    update!(a, 2.0 => (2.0, 0.0))
    r = interpolate(a, 0.5)
    @test isapprox(r[1], 0.5)
    @test isapprox(r[2], 0.5)
    update!(a, 2.0 => (4.0, 0.0))
end

@testset "CVTV field" begin
    f = CVTV((xi, t) -> xi[1] * xi[2] * t)
    @test isapprox(f([1.0, 2.0], 3.0), 6.0)
end

@testset "Dictionary fields" begin
    X = Dict(1 => [0.0, 0.0], 1000 => [1.0, 0.0], 100000 => [1.0, 1.0])
    G = DVTId(X)
    @test isapprox(G[1], X[1])
    @test isapprox(G[1000], X[1000])
    @test isapprox(G[100000], X[100000])
    Y = Dict(1 => [2.0, 2.0], 1000 => [3.0, 2.0], 100000 => [3.0, 3.0])
    F = DVTVd(0.0 => X, 1.0 => Y)
    @test isapprox(interpolate(F, 0.5)[100000], [2.0, 2.0])
end

@testset "update dictionary field" begin
    f1 = Dict(1 => 1.0, 2 => 2.0, 3 => 3.0)
    f2 = Dict(1 => 2.0, 2 => 3.0, 3 => 4.0)
    fld = DVTVd(0.0 => f1)
    update!(fld, 1.0 => f2)
    @test isapprox(interpolate(fld, 0.5)[1], 1.5)
    update!(fld, 1.0 => f1)
    @test isapprox(interpolate(fld, 0.5)[1], 1.0)
end

@testset "use of common constructor field" begin
    @test isa(field(1.0), DCTI)
    @test isa(field(1.0 => 1.0), DCTV)
    @test isa(field((1.0, 2.0)), DVTI)
    @test isa(field(1, 2), DVTI)
    @test isa(field(1.0 => (1.0, 2.0)), DVTV)
    @test isa(field((xi, t) -> xi[1] * t), CVTV)
    @test isa(field(1 => [1.0, 2.0], 10 => [2.0, 3.0]), DVTId)
    @test isa(field(0.0 => (1 => 1.0, 10 => 2.0), 1.0 => (1 => 2.0, 10 => 3.0)), DVTVd)
    X = Dict(1 => [0.0, 0.0], 2 => [1.0, 0.0])
    X1 = field(X)
    X2 = field(0.0 => X)
    @test isa(X1, DVTId)
    @test isa(X2, DVTVd)
end

@testset "general interpolation" begin
    a = [1, 2, 3]
    b = (2, 3, 4)
    @test interpolate(a, b) == 2 + 6 + 12
    a = (1, 2)
    b = (2, 3, 4)
    @test interpolate(a, b) == 2 + 6
    @test_throws AssertionError interpolate(b, a)
end
