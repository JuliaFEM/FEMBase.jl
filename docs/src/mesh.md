# Mesh

```@meta
CurrentModule = FEMBase
DocTestSetup = quote
    using FEMBase
end
```

## Mesh structure

## Modifying mesh

## Defining new mesh parsers


Whould work somehow like this:

```julia
mutable struct DemoReader <: AbstractMeshReader
    handle :: String
end

function read_mesh!(m::Mesh, r::DemoReader)
    # parse file and insert results to `m`.
end

mesh = Mesh()
reader = DemoReader("file.inp")
read_mesh!(mesh, reader)
```
