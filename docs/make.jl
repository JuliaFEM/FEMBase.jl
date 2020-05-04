# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using Documenter, FEMBase
using Literate

Literate.markdown("test/test_fields.jl", "docs/src"; name="fields", documenter=true)

DEVELOPER_GUIDE = [
          "mesh.md",
          "fields.md",
          "basis.md",
          "integration.md",
          "elements.md",
          "problems.md",
          "solvers.md",
          "postprocessing.md",
          "results.md",
          "materials.md",
          "testing.md"]

LIBRARY = ["api.md"]

PAGES = [
         "Introduction" => "index.md",
        #  "Developing JuliaFEM" => DEVELOPER_GUIDE,
        #  "Library" => LIBRARY,
         "Fields" => "fields.md",
        ]

makedocs(modules=[FEMBase],
         format = Documenter.HTML(analytics = "UA-83590644-1"),
         checkdocs = :all,
         sitename = "FEMBase.jl",
         authors = "Jukka Aho",
         pages = PAGES)
