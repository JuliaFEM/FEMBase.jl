# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using FEMBase: get_problems
using Test

import FEMBase: run!

struct MyAnalysis1 <: AbstractAnalysis end

struct MyAnalysis2 <: AbstractAnalysis end

mutable struct MyProblem <: FieldProblem
    value :: Bool
end

MyProblem() = MyProblem(false)

function run!(analysis::Analysis{MyAnalysis2})
    for problem in get_problems(analysis)
        problem.properties.value = true
    end
end

@testset "test creating new analysis" begin
    analysis1 = Analysis(MyAnalysis1, "my test analysis 1")
    analysis2 = Analysis(MyAnalysis2, "my test analysis 2")
    problem = Problem(MyProblem, "my test problem", 1)
    add_problems!(analysis1, [problem])
    add_problems!(analysis2, [problem])
    run!(analysis1)
    @test problem.properties.value == false
    run!(analysis2)
    @test problem.properties.value == true
end

## Result writers

# Writing results of some analysis is done, in abstract level, by using
# results writers. They must be subtypes of `AbstractResultsWriter`:

mutable struct TestResultsWriter <: AbstractResultsWriter
end

# Results writer is added to analysis using `add_results_writer!`:

writer = TestResultsWriter()
analysis1 = Analysis(MyAnalysis1)
add_results_writer!(analysis1, writer)
@test first(get_results_writers(analysis1)) === writer

# Results are written with function `write_results!`-function, taking analysis
# as input argument. Function applies all results writers (there can be many!)
# to the analysis in the order they are attached to analysis using
# `add_results_writer!`. Internally the function runs
# `write_results!(analysis, results_writer)`, which must be implemented for each
# analysis / results_writer pair.

write_results!(analysis1)

# This of course doesn't do anything useful but inform that writing the results
# of MyAnalysis1 using TestResultsWriter is not implemented as the actual
# implementation is missing. In pseudo-level, the actual implementation would
# be something like:

function FEMBase.write_results!(analysis::MyAnalysis1, writer::TestResultsWriter)
    # fid = open(writer.filename)
    # for problem in get_problems(analysis)
    #     field = problem("displacement", analysis.time)
    #     write(fid, field)
    # end
    return nothing!
end

# Lastly, it is not mandatory to have a results writer in the analysis at all
# or the results writer can write directly to stdout or do some other fancy
# stuff. If the results writer is not defined, user gets only info message
# about missing results writers.

empty!(analysis1.results_writers)
write_results!(analysis1)
