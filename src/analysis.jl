# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

abstract type AbstractAnalysis end
abstract type AbstractResultsWriter end

type Analysis{A<:AbstractAnalysis}
    name :: String
    problems :: Vector{Problem}
    fields :: Dict{String, AbstractField}
    results_writers :: Vector{AbstractResultsWriter}
    properties :: A
end

function Analysis{A<:AbstractAnalysis}(::Type{A}, name::String="$A Analysis")
    analysis = Analysis{A}(name, [], Dict(), [], A())
    return analysis
end

function add_problems!(analysis::Analysis, problems::Vector)
    append!(analysis.problems, problems)
end

function get_problems(analysis::Analysis)
    return analysis.problems
end

function add_results_writer!{W<:AbstractResultsWriter}(analysis::Analysis, writer::W)
    push!(analysis.results_writers, writer)
end

function get_results_writers(analysis::Analysis)
    return analysis.results_writers
end

function run!{A<:AbstractAnalysis}(::Analysis{A})
    info("This is a placeholder function for running an analysis $A for a set of problems.")
end
