# Testing extensions

Own extensions to JuliaFEM can be done to own separate packages which are then
used in JuliaFEM. The main idea is that `FEMBase.jl` is giving all supporting
functions and types for all kind of extensions. Extensions lives in their own
modules. Extensions can be tested using `FEMBase.Test`, which itself is a
extension to FEMBase, containing types introduced in this manual.

Some guidelines:
- use Coverage.jl to check coverage, should be 100 %
- use Lint.jl to check syntax
- no use of tabulators in files allowed
- no use of fancy utf-8 in code
- licence header should match in every source file to one defined in main file
- keep version history clean and understandable
- unit tests should test only implemented functions, nothing else
