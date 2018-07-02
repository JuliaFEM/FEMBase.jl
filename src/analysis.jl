# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

abstract type AbstractAnalysis end
abstract type AbstractResultsWriter end

mutable struct Analysis{A<:AbstractAnalysis}
    name :: String
    problems :: Vector{Problem}
    fields :: Dict{String, AbstractField}
    results_writers :: Vector{AbstractResultsWriter}
    properties :: A
end

function Analysis(::Type{A}, name::String="$A Analysis") where A<:AbstractAnalysis
    analysis = Analysis{A}(name, [], Dict(), [], A())
    return analysis
end

function add_problems!(analysis::Analysis, problems::Vector)
    append!(analysis.problems, problems)
end

function get_problems(analysis::Analysis)
    return analysis.problems
end

function add_results_writer!(analysis::Analysis, writer::W) where W<:AbstractResultsWriter
    push!(analysis.results_writers, writer)
    return nothing
end

function get_results_writers(analysis::Analysis)
    return analysis.results_writers
end

function run!(::Analysis{A}) where A<:AbstractAnalysis
    info("This is a placeholder function for running an analysis $A for a set of problems.")
end

function write_results!(::Analysis{A}, ::W) where {A<:AbstractAnalysis, W<:AbstractResultsWriter}
    info("Writing the results of analysis $A is not supported by a results writer $W")
    return nothing
end

function write_results!(analysis)
    results_writers = get_results_writers(analysis)
    if isempty(results_writers)
        info("No result writers attached to the analysis $(analysis.name). ",
             "In order to get results of the analysis stored to the disk, one ",
             "must attach some results writer to the analysis using ",
             "add_results_writer!, e.g. xdmf_writer = Xdmf(\"results\"); ",
             "add_results_writer!(analysis, xdmf_writer)")
        return nothing
    end
    for results_writer in results_writers
        write_results!(analysis, results_writer)
    end
    return nothing
end
