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
            dt = t1-t0
            new_data = y0*(1-(time-t0)/dt) + y1*(1-(t1-time)/dt)
            return new_data
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
    interpolate(f::DVTV, time, basis)

Interpolate variable, time variant DVTV field in both time and spatial direction.

# Notes
First check that is outside of range -> extrapolate
Secondly check is "exact match" in time
At last, find the correct bin and use linear interpolation

"""
function interpolate{N,T}(field::DVTV{N,T}, time, basis)
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
            Y0 = sum(y0[i]*basis[i] for i=1:N)
            Y1 = sum(y1[i]*basis[i] for i=1:N)
            new_data = Y0*(1-(time-t0)/dt) + Y1*(1-(t1-time)/dt)
            return new_data
        end
    end
end
###

abstract type Discrete<:AbstractField end
abstract type Continuous<:AbstractField end
abstract type Constant<:AbstractField end
abstract type Variable<:AbstractField end
abstract type TimeVariant<:AbstractField end
abstract type TimeInvariant<:AbstractField end

type Field{A<:Union{Discrete, Continuous},
           B<:Union{Constant, Variable},
           C<:Union{TimeVariant, TimeInvariant},
          T}
    data :: T
end

### Different field combinations and other typealiases

const CCTI{T} = Field{Continuous, Constant, TimeInvariant, T}
const CVTI{T} = Field{Continuous, Variable, TimeInvariant, T}
const CCTV{T} = Field{Continuous, Constant, TimeVariant, T}
const CVTV{T} = Field{Continuous, Variable, TimeVariant, T}

function CCTI{T}(a::T)
    return CCTI{T}(a)
end

function CVTI{T}(a::T)
    return CVTI{T}(a)
end

function CCTV{T}(a::T)
    return CCTV{T}(a)
end

function CVTV{T}(a::T)
    return CVTV{T}(a)
end

function Field()
    return DCTI()
end

function Field{T}(data::T)
    return DCTI{T}(data)
end

""" For vector data, DVTI is automatically created.

julia> DVTI([1.0, 2.0]) == Field([1.0, 2.0])
true

"""
function Field(data::Vector)
    return DVTI(data)
end

""" For dictionary data, DVTI is automatically created.

Define e.g. nodal coordinates in dictionary
julia> X = Dict(1 => [1.0, 2.0], 2 => [3.0, 4.0])
julia> Field(X) == DVTI(X)

"""
function Field(data::Dict)
    return DVTI(data)
end

""" Multi-slicing of field.

julia> f = DVTI([1.0, 2.0, 3.0])
julia> f[[1, 3]]
[1.0, 3.0]

"""
function getindex(field::DVTI, I::Array{Int64, 1})
    return [field.data[i] for i in I]
end

function length(field::DVTI)
    return length(field.data)
end

function start(field::DVTI)
    return 1
end

function +(f1::DVTI, f2::DVTI)
    return DVTI(f1.data + f2.data)
end

function -(f1::DVTI, f2::DVTI)
    return DVTI(f1.data - f2.data)
end

function update!(field::DVTI, data::Union{Vector, Dict})
    field.data = data
end

""" Take scalar product of DVTI and constant T. """
function *(T::Number, field::DVTI)
    return DVTI(T*field.data)
end

""" Take dot product of DVTI field and vector T. Vector length must match to the
field length and this can be used mainly for interpolation purposes, i.e., u = ∑ Nᵢuᵢ.
"""
function *(T::Union{Vector, RowVector}, f::DVTI)
    @assert length(T) <= length(f)
    return sum([T[i]*f[i] for i=1:length(T)])
end

""" Take outer product of DVTI field and matrix T. """
function *(T::Matrix, f::DVTI)
    n, m = size(T)
    return sum([kron(T[:,i], f[i]') for i=1:m])'
end

function vec(field::DVTI)
    return [field.data...;]
end

# """ Interpolate discrete, variable, time-variant field in time direction. """
function (field::DVTV)(time::Float64)
    time < first(field).time && return DVTI(first(field).data)
    time > last(field).time && return DVTI(last(field).data)
    for i=reverse(1:length(field))
        isapprox(field[i].time, time) && return DVTI(field[i].data)
    end
    for i=reverse(2:length(field))
        t0 = field[i-1].time
        t1 = field[i].time
        if t0 < time < t1
            y0 = field[i-1].data
            y1 = field[i].data
            dt = t1-t0
            new_data = y0*(1-(time-t0)/dt) + y1*(1-(t1-time)/dt)
            return DVTI(new_data)
        end
    end
end

### Basic data structure for continuous field

type Basis
    basis :: Function
    dbasis :: Function
end

### Convenient functions to create fields

function Field(func::Function)
    if method_exists(func, Tuple{})
        return CCTI(func)
    elseif method_exists(func, Tuple{Float64})
        return CCTV(func)
    elseif method_exists(func, Tuple{Vector})
        return CVTI(func)
    elseif method_exists(func, Tuple{Vector, Float64})
        return CVTV(func)
    else
        error("no proper definition found for function: check methods.")
    end
end

### Accessing continuous fields

function (field::CCTI)(xi::Vector, time::Number)
    return field.data()
end

function (field::CVTI)(xi::Vector, time::Number)
    return field.data(xi)
end

function (field::CCTV)(xi::Vector, time::Number)
    return field.data(time)
end

function (field::CVTV)(xi, time)
    return field.data(xi, time)
end
