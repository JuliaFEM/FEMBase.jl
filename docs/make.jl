# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBase.jl/blob/master/LICENSE

using Documenter, FEMBase

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
         "Home" => "index.md",
         "Developing JuliaFEM" => DEVELOPER_GUIDE,
         "Library" => LIBRARY
        ]

makedocs(modules=[FEMBase],
         format = Documenter.HTML(analytics = "UA-83590644-1"),
         checkdocs = :all,
         sitename = "FEMBase.jl",
         authors = "Jukka Aho",
         pages = PAGES)

include("deploy.jl")
