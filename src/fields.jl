# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

abstract type AbstractField end

function length{F<:AbstractField}(f::F)
    return length(f.data)
end

function size{F<:AbstractField}(f::F)
    return size(f.data)
end

function =={F<:AbstractField}(f::F, y)
    return ==(f.data, y)
end

function getindex{F<:AbstractField}(f::F, i::Int64)
    return getindex(f.data, i)
end

function first{F<:AbstractField}(f::F)
    return first(f.data)
end

function last{F<:AbstractField}(f::F)
    return last(f.data)
end

function endof{F<:AbstractField}(f::F)
    return endof(f.data)
end

# DISCRETE, CONSTANT, TIME INVARIANT

"""
    DCTI(T)

Discrete, constant, time-invariant field.
This field is constant in both spatial direction and time direction,
i.e. df/dX = 0 and df/dt = 0.
"""
type DCTI{T} <: AbstractField
    data :: T
end

"""
    update!(f::DCTI, data)

Update new value to field.
"""
function update!(f::DCTI, data)
    f.data = data
end

"""
    interpolate(f::DCTI, args...)

Interpolate constant, time-invariant DCTI field in time direction. That is
trivially only the data itself.
"""
function interpolate(f::DCTI, args...)
    return f.data
end

function getindex(f::DCTI, i::Int64)
    return f.data
end

## DISCRETE, VARIABLE, TIME INVARIANT

""" Discrete, variable, time-invariant field. This is constant in time direction,
but not in spatial direction, i.e. df/dt = 0 but df/dX != 0. The basic structure
of data is Vector, and it is implicitly assumed that length of field matches to
the number of shape functions, so that interpolation in spatial direction works.
"""
type DVTI{N,T} <: AbstractField
    data :: NTuple{N,T}
end

"""
    update!(f::DVTI, data)

Update new value to field.
"""
function update!(f::DVTI, data)
    f.data = data
end

"""
    interpolate(f::DVTI, t)

Interpolate variable, time-invariant DVTI field in time direction.
"""
function interpolate{N,T}(f::DVTI{N,T}, t)
    return f.data
end

"""
    interpolate(f::DVTI, t, B)

Interpolate variable, time-invariant DVTI field in time and spatial direction.
"""
function interpolate{N,T}(f::DVTI{N,T}, t, B)
    return sum(f.data[i]*B[i] for i=1:N)
end

## DISCRETE, CONSTANT, TIME VARIANT

type DCTV{T} <: AbstractField
    data :: Vector{Pair{Float64,T}}
end

"""
    update!(f::DCTV, time => data)

Update new value to field.
"""
function update!{T}(f::DCTV, data::Pair{Float64, T})
    if isapprox(last(f.data).first, data.first)
        f.data[end] = data
    else
        push!(f.data, data)
    end
end

function DCTV{T}(a::Pair{Float64,T}, b::Pair{Float64,T})
    return DCTV([a, b])
end

"""
    interpolate(f::DCTV, time)

Interpolate constant time variant DCTV field in time direction.

# Notes
First check that is outside of range -> extrapolate
Secondly check is "exact match" in time
At last, find the correct bin and use linear interpolation

"""
function interpolate(field::DCTV, time)
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
            return f*y0 + (1-f)*y1
        end
    end
end

## DISCRETE, VARIABLE, TIME VARIANT

type DVTV{N,T} <: AbstractField
    data :: Vector{Pair{Float64,NTuple{N,T}}}
end

"""
    update!(f::DVTV, time => data)

Update new value to field.
"""
function update!{N,T}(f::DVTV, data::Pair{Float64, NTuple{N,T}})
    if isapprox(last(f.data).first, data.first)
        f.data[end] = data
    else
        push!(f.data, data)
    end
end

function DVTV{N,T}(a::Pair{Float64,NTuple{N,T}}, b::Pair{Float64,NTuple{N,T}})
    return DVTV([a, b])
end

"""
    interpolate(f::DVTV, time)

Interpolate variable, time variant DVTV field in time direction.

# Notes
First check that is outside of range -> extrapolate
Secondly check is "exact match" in time
At last, find the correct bin and use linear interpolation

"""
function interpolate{N,T}(field::DVTV{N,T}, time)
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
            new_data = tuple((f*y0[i] + (1-f)*y1[i] for i=1:N)...)
            return new_data
        end
    end
end

"""
    interpolate(f::DVTV, time, basis)

Interpolate variable, time variant DVTV field in both time and spatial direction.

# Notes
First check that is outside of range -> extrapolate
Secondly check is "exact match" in time
At last, find the correct bin and use linear interpolation

"""
function interpolate{N,T}(field::DVTV{N,T}, time, basis)
    data = interpolate(field, time)
    return sum(data[i]*basis[i] for i=1:N)
end

"""
    CVTV(f)

Continuous, variable, time variant field.
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

"""
    interpolate(f::DVTId, time)

Interpolate DVTId, returns trivially the content as this is time invariant field.
"""
function interpolate(f::DVTId, time)
    return f.data
end

"""
    update!(field::DVTId, data::Dict)

Update data to field.
"""
function update!{T}(field::DVTId{T}, data::Dict{Int64, T})
    merge!(field.data, data)
end

"""
    DVTVd(time => data::Dict)

Discrete, variable, time variant dictionary field.
"""
type DVTVd{T} <: AbstractField
    data :: Vector{Pair{Float64,Dict{Int64,T}}}
end

function DVTVd{T}(a::Pair{Float64,Dict{Int64,T}}, b::Pair{Float64,Dict{Int64,T}})
    return DVTVd([a, b])
end

"""
    interpolate(f::DVTVd, time)

Interpolate variable, time variant DVTVd dictionary field in time direction.

# Notes
First check that is outside of range -> extrapolate
Secondly check is "exact match" in time
At last, find the correct bin and use linear interpolation

"""
function interpolate{T}(field::DVTVd{T}, time)
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

"""
    field(x)

Create new discrete, constant, time invariant field from value `x`.

# Example

```julia
f = field(1.0)
```
"""
function field(x)
    return DCTI(x)
end

"""
    field(x::NTuple{N,T})

Create new discrete, variable, time invariant field from tuple `x`.

# Example

```julia
f = field( (1.0, 2.0) )
```
"""
function field{N,T}(data::NTuple{N,T})
    return DVTI(data)
end

"""
    field(x::Pair{Float64,T})

Create new discrete, constant, time variant field from pair `x`.

# Example

```julia
f = field(1.0=>1.0, 2.0=>2.0)
```
"""
function field{T}(data::Pair{Float64,T}...)
    return DCTV(collect(data))
end

"""
    field(x::Pair{Float64,NTuple{N,T})

Create new discrete, variable, time variant field from pair `x`.

# Example

```julia
f = field(1.0=>(1.0,2.0), 2.0=>(2.0,3.0))
```
"""
function field{N,T}(data::Pair{Float64,NTuple{N,T}}...)
    return DVTV(collect(data))
end

"""
    field(x::Function)

Create new, continuous, variable, time variant field from function `x`.

# Example

```julia
f = field( (xi,t) -> xi[1]*t )
```
"""
function field(data::Function)
    return CVTV(data)
end

"""
    field(x::Pair{Int64,T})

Create new discrete, variable, time invariant dictionary field from `x`.

# Example

```julia
f = field(1 => 1.0, 2 => 2.0)
```
"""
function field{T}(data::Pair{Int64, T}...)
    return DVTId(Dict(data))
end

"""
    field(x::Pair{Float64, NTuple{N, Pair{Int64, T}}})

Create new discrete, variable, time variant dictionary field from `x`.

# Example

```julia
X1 = (1 => 1.0, 2 => 2.0)
X2 = (1 => 2.0, 2 => 3.0)
f = field(0.0 => X1, 1.0 => X2)
```
"""
function field{N,T}(data::Pair{Float64, NTuple{N, Pair{Int64, T}}}...)
    return DVTVd(collect(t => Dict(d) for (t, d) in data))
end

"""
    field(x::Dict)

Create new discrete, variable, time invariant dictionary field from dictionary `x`.

# Example

```julia
X = Dict(1 => [0.0,0.0], 2 => [1.0,0.0])
f = field(X)
```
"""
function field{T}(data::Dict{Int64,T})
    return DVTId(data)
end

"""
    field(t::Float64 => x::Dict, ...)

Create new discrete, variable, time variant dictionary field from pair of time
and dictionary.
"""
function field{T}(data::Pair{Float64, Dict{Int64, T}}...)
    return DVTVd(collect(data))
end
