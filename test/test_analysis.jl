# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using FEMBase
using FEMBase.get_problems
using Base.Test

import FEMBase: run!

type MyAnalysis1 <: AbstractAnalysis end

type MyAnalysis2 <: AbstractAnalysis end

type TestResultsWriter <: AbstractResultsWriter
end

type MyProblem <: FieldProblem
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

@testset "test results writer" begin
    writer = TestResultsWriter()
    analysis1 = Analysis(MyAnalysis1)
    add_results_writer!(analysis1, writer)
    @test first(get_results_writers(analysis1)) === writer
end
