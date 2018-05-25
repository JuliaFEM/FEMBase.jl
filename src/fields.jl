# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

"""
    AbstractField

Abstract supertype for all fields in JuliaFEM.
"""
abstract type AbstractField end

function length{F<:AbstractField}(f::F)
    return length(f.data)
end

function size{F<:AbstractField}(f::F)
    return size(f.data)
end

function =={F<:AbstractField}(x::F, y)
    return ==(x.data, y)
end

function =={F<:AbstractField}(x, y::F)
    return ==(x, y.data)
end

function =={F<:AbstractField}(x::F, y::F)
    return ==(x.data, y.data)
end

function getindex{F<:AbstractField}(f::F, i::Int64)
    return getindex(f.data, i)
end

function interpolate_field(field, ::Any)
    return field.data
end

function update_field!(field, data)
    field.data = data
end

"""
    DCTI{T} <: AbstractField

Discrete, constant, time-invariant field.

This field is constant in both spatial direction and time direction,
i.e. df/dX = 0 and df/dt = 0.

# Example

```jldoctest
julia> DCTI(1)
FEMBase.DCTI{Int64}(1)
```
"""
type DCTI{T} <: AbstractField
    data :: T
end

function getindex(field::DCTI, ::Int64)
    return field.data
end

"""
    DVTI{N,T} <: AbstractField

Discrete, variable, time-invariant field.

This is constant in time direction, but not in spatial direction, i.e. df/dt = 0
but df/dX != 0. The basic structure of data is `Tuple`, and it is implicitly
assumed that length of field matches to the number of shape functions, so that
interpolation in spatial direction works.

# Example

```jldoctest
julia> DVTI(1, 2, 3)
FEMBase.DVTI{3,Int64}((1, 2, 3))
```
"""
type DVTI{N,T} <: AbstractField
    data :: NTuple{N,T}
end

function DVTI(data...)
    return DVTI(data)
end

"""
    DCTV{T} <: AbstractField

Discrete, constant, time variant field. This type of field can change in time
direction but not in spatial direction.

# Example

Field having value 5 at time 0.0 and value 10 at time 1.0:

```jldoctest
julia> DCTV(0.0 => 5, 1.0 => 10)
FEMBase.DCTV{Int64}(Pair{Float64,Int64}[0.0=>5, 1.0=>10])
```

"""
type DCTV{T} <: AbstractField
    data :: Vector{Pair{Float64,T}}
end

function DCTV{T}(data::Pair{Float64,T}...)
    return DCTV(collect(data))
end

function update_field!{T}(f::DCTV, data::Pair{Float64, T})
    if isapprox(last(f.data).first, data.first)
        f.data[end] = data
    else
        push!(f.data, data)
    end
end

function interpolate_field(field::DCTV, time)
    time < first(field.data).first && return first(field.data).second
    time > last(field.data).first && return last(field.data).second
    for i=reverse(1:length(field))
        isapprox(field.data[i].first, time) && return field.data[i].second
    end
    for i=length(field.data):-1:2
        t0 = field.data[i-1].first
        t1 = field.data[i].first
        if t0 < time < t1
            y0 = field.data[i-1].second
            y1 = field.data[i].second
            dy = y1-y0
            dt = t1-t0
            return y0 + (time-t0)*dy/dt
        end
    end
end

"""
    DVTV{N,T} <: AbstractField

Discrete, variable, time variant field. The most general discrete field can
change in both temporal and spatial direction.

# Example

```jldoctest
julia> DVTV(0.0 => (1, 2), 1.0 => (2, 3))
FEMBase.DVTV{2,Int64}(Pair{Float64,Tuple{Int64,Int64}}[0.0=>(1, 2), 1.0=>(2, 3)])
```
"""
type DVTV{N,T} <: AbstractField
    data :: Vector{Pair{Float64,NTuple{N,T}}}
end

function DVTV{N,T}(data::Pair{Float64,NTuple{N,T}}...)
    return DVTV(collect(data))
end

function update_field!{N,T}(f::DVTV, data::Pair{Float64, NTuple{N,T}})
    if isapprox(last(f.data).first, data.first)
        f.data[end] = data
    else
        push!(f.data, data)
    end
end

function interpolate_field{N,T}(field::DVTV{N,T}, time)
    time < first(field.data).first && return first(field.data).second
    time > last(field.data).first && return last(field.data).second
    for i=reverse(1:length(field))
        isapprox(field.data[i].first, time) && return field.data[i].second
    end
    for i=length(field.data):-1:2
        t0 = field.data[i-1].first
        t1 = field.data[i].first
        if t0 < time < t1
            y0 = field.data[i-1].second
            y1 = field.data[i].second
            dt = t1-t0
            return map((a,b) -> a + (time-t0)*(b-a)/dt, y0, y1)
        end
    end
end

"""
    CVTV <: AbstractField

Continuous, variable, time variant field.

# Example

```jldoctest
julia> f = CVTV((xi,t) -> xi*t)
FEMBase.CVTV(#1)
```
"""
type CVTV <: AbstractField
    data :: Function
end

function (f::CVTV)(xi, time)
    return f.data(xi, time)
end

"""
    DVTId(X::Dict)

Discrete, variable, time invariant dictionary field.
"""
type DVTId{T} <: AbstractField
    data :: Dict{Int64, T}
end

function update_field!{T}(field::DVTId{T}, data::Dict{Int64, T})
    merge!(field.data, data)
end

"""
    DVTVd(time => data::Dict)

Discrete, variable, time variant dictionary field.
"""
type DVTVd{T} <: AbstractField
    data :: Vector{Pair{Float64,Dict{Int64,T}}}
end

function DVTVd{T}(data::Pair{Float64,Dict{Int64,T}}...)
    return DVTVd(collect(data))
end

function interpolate_field{T}(field::DVTVd{T}, time)
    time < first(field.data).first && return first(field.data).second
    time > last(field.data).first && return last(field.data).second
    for i=reverse(1:length(field))
        isapprox(field.data[i].first, time) && return field.data[i].second
    end
    for i=length(field.data):-1:2
        t0 = field.data[i-1].first
        t1 = field.data[i].first
        if t0 < time < t1
            y0 = field.data[i-1].second
            y1 = field.data[i].second
            f = (time-t0)/(t1-t0)
            new_data = similar(y0)
            for i in keys(y0)
                new_data[i] = f*y0[i] + (1-f)*y1[i]
            end
            return new_data
        end
    end
end

function update_field!{T}(f::DVTVd, data::Pair{Float64,Dict{Int64,T}})
    if isapprox(last(f.data).first, data.first)
        f.data[end] = data
    else
        push!(f.data, data)
    end
end

function new_field(data)
    return DCTI(data)
end

function new_field(data...)
    return DVTI(data)
end

function new_field{N,T}(data::NTuple{N,T})
    return DVTI(data)
end

function new_field{T}(data::Pair{Float64,T}...)
    return DCTV(collect(data))
end

function new_field{N,T}(data::Pair{Float64,NTuple{N,T}}...)
    return DVTV(collect(data))
end

function new_field(data::Function)
    return CVTV(data)
end

function new_field{T}(data::Pair{Int64, T}...)
    return DVTId(Dict(data))
end

function new_field{N,T}(data::Pair{Float64, NTuple{N, Pair{Int64, T}}}...)
    return DVTVd(collect(t => Dict(d) for (t, d) in data))
end

function new_field{T}(data::Dict{Int64,T})
    return DVTId(data)
end

function new_field{T}(data::Pair{Float64, Dict{Int64, T}}...)
    return DVTVd(collect(data))
end

"""
    field(x)

Create new field. Field type is deduced from data type.
"""
function field(data...)
    return new_field(data...)
end

"""
    interpolate(field, time)

Interpolate field in time direction.

# Examples

For time invariant fields [`DCTI`](@ref), [`DVTI`](@ref), [`DVTId`](@ref)
solution is trivially the data inside field as fields does not depend from
the time:

```jldoctest
julia> a = field(1.0)
FEMBase.DCTI{Float64}(1.0)

julia> interpolate(a, 0.0)
1.0
```

```jldoctest
julia> a = field((1.0, 2.0))
FEMBase.DVTI{2,Float64}((1.0, 2.0))

julia> interpolate(a, 0.0)
(1.0, 2.0)
```

```jldoctest
julia> a = field(1=>1.0, 2=>2.0)
FEMBase.DVTId{Float64}(Dict(2=>2.0,1=>1.0))

julia> interpolate(a, 0.0)
Dict{Int64,Float64} with 2 entries:
  2 => 2.0
  1 => 1.0
```

DVTId trivial solution is returned. For time variant fields DCTV, DVTV, DVTVd
linear interpolation is performed.

# Other notes

First algorithm checks that is time out of range, i.e. time is smaller than
time of first frame or larger than last frame. If that is the case, return
first or last frame. Secondly algorithm finds is given time exact match to
time of some frame and return that frame. At last, we find correct bin so
that t0 < time < t1 and use linear interpolation.

"""
function interpolate{F<:AbstractField}(field::F, time)
    return interpolate_field(field, time)
end

"""
    interpolate(a, b)

A helper function for interpolate routines. Given iterables `a` and `b`,
calculate c = aᵢbᵢ. Length of `a` can be less than `b`, but not vice versa.
"""
function interpolate(a, b)
    @assert length(a) <= length(b)
    return sum(a[i]*b[i] for i=1:length(a))
end

"""
    update!(field, data)

Update new value to field.
"""
function update!{F<:AbstractField}(field::F, data)
    update_field!(field, data)
end
