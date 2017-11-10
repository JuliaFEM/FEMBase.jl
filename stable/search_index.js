var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#FEMBase.jl-documentation-1",
    "page": "Introduction",
    "title": "FEMBase.jl documentation",
    "category": "section",
    "text": "Pages = [\"index.md\", \"theory.md\", \"api.md\"]"
},

{
    "location": "theory.html#",
    "page": "Theory",
    "title": "Theory",
    "category": "page",
    "text": ""
},

{
    "location": "theory.html#Theory-1",
    "page": "Theory",
    "title": "Theory",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#API-documentation-1",
    "page": "API",
    "title": "API documentation",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#FEMBase.Assembly",
    "page": "API",
    "title": "FEMBase.Assembly",
    "category": "Type",
    "text": "General linearized problem to solve     (K₁+K₂)Δu  +   C1'Δλ = f₁+f₂          C2Δu  +     DΔλ = g\n\n\n\n"
},

{
    "location": "api.html#FEMBase.Element-Union{Tuple{E}, Tuple{Type{E},Array{Int64,1}}} where E<:FEMBasis.AbstractBasis",
    "page": "API",
    "title": "FEMBase.Element",
    "category": "Method",
    "text": "Element(element_type, connectivity_vector)\n\nConstruct a new element where element_type is the type of the element and connectivity_vector is the vector of nodes that the element is connected to.\n\nExamples\n\nIn the example a new element (E in the figure below) of type Tri3 is created. This spesific element connects to nodes 89, 43, 12 in the finite element mesh.\n\nelement = Element(Tri3, [89, 43, 12])\n\n(Image: img)\n\n\n\n"
},

{
    "location": "api.html#FEMBase.Problem",
    "page": "API",
    "title": "FEMBase.Problem",
    "category": "Type",
    "text": "Defines types for Problem variables.\n\nExamples\n\nThe type of 'elements' is Vector{Element}\n\nAdd elements into the Problem element list.\n\na = [1, 2, 3]\nProblem.elements = a\n\n\n\n"
},

{
    "location": "api.html#FEMBase.Problem-Union{Tuple{P}, Tuple{Type{P},AbstractString,Int64}} where P<:FEMBase.FieldProblem",
    "page": "API",
    "title": "FEMBase.Problem",
    "category": "Method",
    "text": "Problem(problem_type, problem_name::String, problem_dimension)\n\nConstruct a new field problem where problem_type is the type of the problem (Elasticity, Dirichlet, etc.), problem_name is the name of the problem and problem_dimension is the number of DOF:s in one node (2 in a 2D problem, 3 in an elastic 3D problem, 6 in a 3D beam problem, etc.).\n\nExamples\n\nCreate a vector-valued (dim=3) elasticity problem:\n\nprob1 = Problem(Elasticity, \"this is my problem\", 3)\n\n\n\n"
},

{
    "location": "api.html#FEMBase.Problem-Union{Tuple{P}, Tuple{Type{P},Any,Any,Any}} where P<:FEMBase.BoundaryProblem",
    "page": "API",
    "title": "FEMBase.Problem",
    "category": "Method",
    "text": "Construct a new boundary problem.\n\nExamples\n\nCreate a Dirichlet boundary problem for a vector-valued (dim=3) elasticity problem.\n\njulia> bc1 = Problem(Dirichlet, \"support\", 3, \"displacement\") solver.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.add!",
    "page": "API",
    "title": "FEMBase.add!",
    "category": "Function",
    "text": "Add new data to COO Sparse vector. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.add!-Tuple{FEMBase.SparseMatrixCOO,Array{Int64,1},Array{Int64,1},Array{T,2} where T}",
    "page": "API",
    "title": "FEMBase.add!",
    "category": "Method",
    "text": "Add local element matrix to sparse matrix. This basically does:\n\nA[dofs1, dofs2] = A[dofs1, dofs2] + data\n\nExample\n\nS = [3, 4] M = [6, 7, 8] data = Float64[5 6 7; 8 9 10] A = SparseMatrixCOO() add!(A, S, M, data) full(A)\n\n4x8 Array{Float64,2}:  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0  5.0  6.0   7.0  0.0  0.0  0.0  0.0  0.0  8.0  9.0  10.0\n\n\n\n"
},

{
    "location": "api.html#FEMBase.add!-Tuple{FEMBase.SparseMatrixCOO,SparseMatrixCSC}",
    "page": "API",
    "title": "FEMBase.add!",
    "category": "Method",
    "text": "Add sparse matrix of CSC to COO. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.add!-Tuple{FEMBase.SparseMatrixCOO,SparseVector}",
    "page": "API",
    "title": "FEMBase.add!",
    "category": "Method",
    "text": "Add SparseVector to SparseVectorCOO. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.add_elements!-Tuple{FEMBase.Problem,Any}",
    "page": "API",
    "title": "FEMBase.add_elements!",
    "category": "Method",
    "text": "add_elements!(problem::Problem, elements)\n\nAdd new elements into the problem.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_gdofs-Tuple{FEMBase.Problem,FEMBase.Element}",
    "page": "API",
    "title": "FEMBase.get_gdofs",
    "category": "Method",
    "text": "Return global degrees of freedom for element.\n\nNotes\n\nFirst look dofs from problem.dofmap, it not found, update dofmap from element.element connectivity using formula gdofs = [dim*(nid-1)+j for j=1:dim]\n\nlook element dofs from problem.dofmap\nif not found, use element.connectivity to update dofmap and 1.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_integration_points-Union{Tuple{E}, Tuple{FEMBase.Element{E},Int64}} where E",
    "page": "API",
    "title": "FEMBase.get_integration_points",
    "category": "Method",
    "text": "This is a special case, temporarily change order of integration scheme mainly for mass matrix.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_unknown_field_dimension-Tuple{FEMBase.Problem}",
    "page": "API",
    "title": "FEMBase.get_unknown_field_dimension",
    "category": "Method",
    "text": "Return the dimension of the unknown field of this problem. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_unknown_field_name-Union{Tuple{FEMBase.Problem{P}}, Tuple{P}} where P",
    "page": "API",
    "title": "FEMBase.get_unknown_field_name",
    "category": "Method",
    "text": "Return the name of the unknown field of this problem. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.group_by_element_type-Tuple{Array{FEMBase.Element,1}}",
    "page": "API",
    "title": "FEMBase.group_by_element_type",
    "category": "Method",
    "text": "group_by_element_type(elements::Vector{Element})\n\nGiven a vector of elements, group elements by element type to several vectors. Returns a dictionary, where key is the element type and value is a vector containing all elements of type element_type.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.initialize!-Tuple{FEMBase.Problem,FEMBase.Element,Float64}",
    "page": "API",
    "title": "FEMBase.initialize!",
    "category": "Method",
    "text": "function initialize!(problem_type, element_name, time)\n\nInitialize the element ready for calculation, where problem_type is the type of the problem (Elasticity, Dirichlet, etc.), element_name is the name of a constructed element (see Element(element_type, connectivity_vector)) and time is the starting time of the initializing process.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Tuple{FEMBase.Problem,FEMBase.Assembly,Array{T,1} where T,Array{T,1} where T}",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "update!(problem, assembly, u, la)\n\nUpdate the problem solution vector for assembly.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Union{Tuple{E}, Tuple{FEMBase.Element{E},AbstractString,Dict}} where E",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "Update element field based on a dictionary of nodal data and connectivity information.\n\nExamples\n\njulia> data = Dict(1 => [0.0, 0.0], 2 => [1.0, 2.0]) julia> element = Seg2([1, 2]) julia> update!(element, \"geometry\", data)\n\nAs a result element now have time invariant (variable) vector field \"geometry\" with data ([0.0, 0.0], [1.0, 2.0]).\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Union{Tuple{FEMBase.Problem{P},FEMBase.Assembly,Array{FEMBase.Element,1},Float64}, Tuple{P}} where P<:FEMBase.FieldProblem",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "update!(problem, assembly, elements, time)\n\nUpdate a solution from the assebly to elements.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Union{Tuple{P,Vararg{Pair{String,String},N} where N}, Tuple{P}} where P<:FEMBase.AbstractProblem",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "update!(problem.properties, attr...)\n\nUpdate properties for a problem.\n\nExample\n\nupdate!(body.properties, \"finite_strain\" => \"false\")\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Union{Tuple{T}, Tuple{Union{FEMBase.Field{FEMBase.Discrete,FEMBase.Constant,FEMBase.TimeVariant,T} where T, FEMBase.Field{FEMBase.Discrete,FEMBase.Variable,FEMBase.TimeVariant,T} where T},Pair{Float64,T}}} where T",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "Update time-dependent fields with new values.\n\nExamples\n\njulia> f = Field(0.0 => 1.0) julia> update!(f, 1.0 => 2.0)\n\nNow field has two (time, value) pairs: (0.0, 1.0) and (1.0, 2.0)\n\nNotes\n\nTime vector is assumed to be ordered t_i-1 < t_i < t_i+1. If updating field with already existing time the old value is replaced with new one.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.DCTI-Tuple{}",
    "page": "API",
    "title": "FEMBase.DCTI",
    "category": "Method",
    "text": "Discrete, constant, time-invariant field. This is constant in both spatial direction and time direction, i.e. df/dX = 0 and df/dt = 0.\n\nThis is the most basic type of field having no anything special functionality.\n\nExamples\n\njulia> f = DCTI() julia> update!(f, 1.0)\n\nMultiplying by constant works:\n\njulia> 2*f 2.0\n\nInterpolation in time direction gives the same constant:\n\njulia> f(1.0) 1.0\n\nBy default, when calling Field with scalar, DCTI is assumed, i.e.\n\njulia> Field(0.0) == DCTI(0.0) true\n\n\n\n"
},

{
    "location": "api.html#FEMBase.DCTV-Tuple{Vararg{Pair,N} where N}",
    "page": "API",
    "title": "FEMBase.DCTV",
    "category": "Method",
    "text": "Discrete, constant, time variant field. This is constant in spatial direction but non-constant in time direction, i.e. df/dX = 0 but df/dt != 0.\n\nExamples\n\njulia> t0 = 0.0; t1=1.0; y0 = 0.0; y1 = 1.0 julia> f = DCTV(t0 => y0, t1 => y1)\n\n\n\n"
},

{
    "location": "api.html#FEMBase.DVTI-Tuple{}",
    "page": "API",
    "title": "FEMBase.DVTI",
    "category": "Method",
    "text": "Discrete, variable, time-invariant field. This is constant in time direction, but not in spatial direction, i.e. df/dt = 0 but df/dX != 0. The basic structure of data is Vector, and it is implicitly assumed that length of field matches to the number of shape functions, so that interpolation in spatial direction works.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.DVTV-Tuple{}",
    "page": "API",
    "title": "FEMBase.DVTV",
    "category": "Method",
    "text": "Discrete, variable, time variant fields. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.Field-Tuple{Array{T,1} where T}",
    "page": "API",
    "title": "FEMBase.Field",
    "category": "Method",
    "text": "For vector data, DVTI is automatically created.\n\njulia> DVTI([1.0, 2.0]) == Field([1.0, 2.0]) true\n\n\n\n"
},

{
    "location": "api.html#FEMBase.Field-Tuple{Dict}",
    "page": "API",
    "title": "FEMBase.Field",
    "category": "Method",
    "text": "For dictionary data, DVTI is automatically created.\n\nDefine e.g. nodal coordinates in dictionary julia> X = Dict(1 => [1.0, 2.0], 2 => [3.0, 4.0]) julia> Field(X) == DVTI(X)\n\n\n\n"
},

{
    "location": "api.html#FEMBase.Increment",
    "page": "API",
    "title": "FEMBase.Increment",
    "category": "Type",
    "text": "Simple time frame / increment to contain both time and data. \n\n\n\n"
},

{
    "location": "api.html#Base.:*-Tuple{Array{T,2} where T,FEMBase.Field{FEMBase.Discrete,FEMBase.Constant,FEMBase.TimeInvariant,T} where T}",
    "page": "API",
    "title": "Base.:*",
    "category": "Method",
    "text": "Kind of spatial interpolation of DCTI. \n\n\n\n"
},

{
    "location": "api.html#Base.:*-Tuple{Array{T,2} where T,FEMBase.Field{FEMBase.Discrete,FEMBase.Variable,FEMBase.TimeInvariant,T} where T}",
    "page": "API",
    "title": "Base.:*",
    "category": "Method",
    "text": "Take outer product of DVTI field and matrix T. \n\n\n\n"
},

{
    "location": "api.html#Base.:*-Tuple{Number,FEMBase.Field{FEMBase.Discrete,FEMBase.Variable,FEMBase.TimeInvariant,T} where T}",
    "page": "API",
    "title": "Base.:*",
    "category": "Method",
    "text": "Take scalar product of DVTI and constant T. \n\n\n\n"
},

{
    "location": "api.html#Base.:*-Tuple{Union{Array{T,1} where T, RowVector},FEMBase.Field{FEMBase.Discrete,FEMBase.Variable,FEMBase.TimeInvariant,T} where T}",
    "page": "API",
    "title": "Base.:*",
    "category": "Method",
    "text": "Take dot product of DVTI field and vector T. Vector length must match to the field length and this can be used mainly for interpolation purposes, i.e., u = ∑ Nᵢuᵢ.\n\n\n\n"
},

{
    "location": "api.html#Base.SparseArrays.sparse-Tuple{FEMBase.SparseMatrixCOO}",
    "page": "API",
    "title": "Base.SparseArrays.sparse",
    "category": "Method",
    "text": "Convert from COO format to CSC.\n\nParameters\n\ntol     used to drop near zero values less than tol.\n\n\n\n"
},

{
    "location": "api.html#Base.getindex-Tuple{FEMBase.Field{FEMBase.Discrete,FEMBase.Variable,FEMBase.TimeInvariant,T} where T,Array{Int64,1}}",
    "page": "API",
    "title": "Base.getindex",
    "category": "Method",
    "text": "Multi-slicing of field.\n\njulia> f = DVTI([1.0, 2.0, 3.0]) julia> f[[1, 3]] [1.0, 3.0]\n\n\n\n"
},

{
    "location": "api.html#Base.getindex-Tuple{FEMBase.Field{FEMBase.Discrete,FEMBase.Variable,FEMBase.TimeInvariant,T} where T,Int64}",
    "page": "API",
    "title": "Base.getindex",
    "category": "Method",
    "text": "Default slicing of field.\n\njulia> f = DVTI([1.0, 2.0]) julia> f[1] 1.0\n\n\n\n"
},

{
    "location": "api.html#Base.haskey-Tuple{FEMBase.Element,String}",
    "page": "API",
    "title": "Base.haskey",
    "category": "Method",
    "text": "Check existence of field. \n\n\n\n"
},

{
    "location": "api.html#Base.isapprox-Tuple{FEMBase.SparseMatrixCOO,FEMBase.SparseMatrixCOO}",
    "page": "API",
    "title": "Base.isapprox",
    "category": "Method",
    "text": "Approximative comparison of two matricse A and B. \n\n\n\n"
},

{
    "location": "api.html#Base.length-Union{Tuple{B}, Tuple{FEMBase.Element{B}}} where B",
    "page": "API",
    "title": "Base.length",
    "category": "Method",
    "text": "length(element::Element)\n\nReturn the number of nodes in element.\n\n\n\n"
},

{
    "location": "api.html#Base.similar-Tuple{FEMBase.Field{FEMBase.Discrete,FEMBase.Variable,FEMBase.TimeInvariant,T} where T,Array{T,1} where T}",
    "page": "API",
    "title": "Base.similar",
    "category": "Method",
    "text": "Create a similar DVTI field from vector data. \n\njulia> f1 = DVTI(Vector[[1.0, 2.0], [3.0, 4.0]]) julia> f2 = similar(f1, [2.0, 3.0, 4.0, 5.0]) julia> f2 == DVTI(Vector[[2.0, 3.0], [4.0, 5.0]]) true\n\n\n\n"
},

{
    "location": "api.html#FEMBase.assemble_mass_matrix!-Tuple{FEMBase.Problem,Array{FEMBase.Element{FEMBasis.Tet10},1},Any}",
    "page": "API",
    "title": "FEMBase.assemble_mass_matrix!",
    "category": "Method",
    "text": "assemble_mass_matrix!(problem, elements::Vector{Element{Tet10}}, time)\n\nAssemble Tet10 mass matrices using special method. If Tet10 has constant metric if can be integrated analytically to gain performance.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_dualbasis",
    "page": "API",
    "title": "FEMBase.get_dualbasis",
    "category": "Function",
    "text": "Return dual basis transformation matrix Ae. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_global_solution-Tuple{FEMBase.Problem,FEMBase.Assembly}",
    "page": "API",
    "title": "FEMBase.get_global_solution",
    "category": "Method",
    "text": "get_global_solution(problem, assembly)\n\nReturn a global solution (u, la) for a problem.\n\nNotes\n\nIf the length of solution vector != number of nodes, i.e. the field dimension is something else than 1, reshape vectors so that their length matches to the number of nodes. This helps to get nodal results easily.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_local_coordinates-Tuple{FEMBase.Element,Array{T,1} where T,Float64}",
    "page": "API",
    "title": "FEMBase.get_local_coordinates",
    "category": "Method",
    "text": "Find inverse isoparametric mapping of element. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_nonzero_rows-Tuple{SparseMatrixCSC}",
    "page": "API",
    "title": "FEMBase.get_nonzero_rows",
    "category": "Method",
    "text": "Find all nonzero rows from sparse matrix.\n\nReturns\n\nOrdered list of row indices.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_parent_field_name-Union{Tuple{FEMBase.Problem{P}}, Tuple{P}} where P<:FEMBase.BoundaryProblem",
    "page": "API",
    "title": "FEMBase.get_parent_field_name",
    "category": "Method",
    "text": "Return the name of the parent field of this (boundary) problem. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.inside-Union{Tuple{E}, Tuple{FEMBase.Element{E},Any,Any}} where E",
    "page": "API",
    "title": "FEMBase.inside",
    "category": "Method",
    "text": "Test is X inside element. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.optimize!-Tuple{FEMBase.SparseMatrixCOO}",
    "page": "API",
    "title": "FEMBase.optimize!",
    "category": "Method",
    "text": "Combine (I,J,V) values if possible to reduce memory usage. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.resize_sparse-Tuple{Any,Any,Any}",
    "page": "API",
    "title": "FEMBase.resize_sparse",
    "category": "Method",
    "text": "Resize sparse matrix A to (higher) dimension n x m. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.resize_sparsevec-Tuple{Any,Any}",
    "page": "API",
    "title": "FEMBase.resize_sparsevec",
    "category": "Method",
    "text": "Resize sparse vector b to (higher) dimension n. \n\n\n\n"
},

{
    "location": "api.html#FEMBasis.BasisInfo",
    "page": "API",
    "title": "FEMBasis.BasisInfo",
    "category": "Type",
    "text": "Data type for fast FEM.\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.BasisInfo-Union{Tuple{B}, Tuple{Type{B},Any}, Tuple{Type{B}}} where B<:FEMBasis.AbstractBasis",
    "page": "API",
    "title": "FEMBasis.BasisInfo",
    "category": "Method",
    "text": "Initialization of data type BasisInfo.\n\nExamples\n\n\nBasisInfo(Tri3)\n\n# output\n\nFEMBasis.BasisInfo{FEMBasis.Tri3,Float64}([0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0], 0.0)\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.NSeg",
    "page": "API",
    "title": "FEMBasis.NSeg",
    "category": "Type",
    "text": "NURBS segment. \n\n\n\n"
},

{
    "location": "api.html#FEMBasis.eval_basis!-Union{Tuple{B}, Tuple{FEMBasis.BasisInfo{B,T} where T,Any,Any}} where B",
    "page": "API",
    "title": "FEMBasis.eval_basis!",
    "category": "Method",
    "text": "Evaluate basis, gradient and so on for some point xi.\n\nExamples\n\n\nb = BasisInfo(Quad4)\nX = ((0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0))\nxi = (0.0, 0.0)\neval_basis!(b, X, xi)\n\n# output\n\nFEMBasis.BasisInfo{FEMBasis.Quad4,Float64}([0.25 0.25 0.25 0.25], [-0.25 0.25 0.25 -0.25; -0.25 -0.25 0.25 0.25], [-0.5 0.5 0.5 -0.5; -0.5 -0.5 0.5 0.5], [0.5 0.0; 0.0 0.5], [2.0 -0.0; -0.0 2.0], 0.25)\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.grad!-Union{Tuple{B}, Tuple{FEMBasis.BasisInfo{B,T} where T,Any,Any}} where B",
    "page": "API",
    "title": "FEMBasis.grad!",
    "category": "Method",
    "text": "grad!(bi, gradu, u)\n\nEvalute gradient ∂u/∂X and store result to matrix gradu. It is assumed that eval_basis! has been already run to bi so it already contains all necessary matrices evaluated with some X and xi.\n\nExample\n\nFirst setup and evaluate basis using eval_basis!:\n\nB = BasisInfo(Quad4)\nX = ((0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0))\nxi = (0.0, 0.0)\neval_basis!(B, X, xi)\n\n# output\n\nFEMBasis.BasisInfo{FEMBasis.Quad4,Float64}([0.25 0.25 0.25 0.25], [-0.25 0.25 0.25 -0.25; -0.25 -0.25 0.25 0.25], [-0.5 0.5 0.5 -0.5; -0.5 -0.5 0.5 0.5], [0.5 0.0; 0.0 0.5], [2.0 -0.0; -0.0 2.0], 0.25)\n\n\nNext, calculate gradient of u:\n\nu = ((0.0, 0.0), (1.0, -1.0), (2.0, 3.0), (0.0, 0.0))\ngradu = zeros(2, 2)\ngrad!(B, gradu, u)\n\n# output\n\n2×2 Array{Float64,2}:\n 1.5  0.5\n 1.0  2.0\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.grad-NTuple{4,Any}",
    "page": "API",
    "title": "FEMBasis.grad",
    "category": "Method",
    "text": "grad(B, T, X, xi)\n\nCalculate gradient of T with respect to X in point xi using basis B.\n\nExample\n\nB = Quad4()\nX = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])\nu = ([0.0, 0.0], [1.0, -1.0], [2.0, 3.0], [0.0, 0.0])\ngrad(B, u, X, (0.0, 0.0))\n\n# output\n\n2×2 Array{Float64,2}:\n 1.5  0.5\n 1.0  2.0\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.grad-Tuple{Any,Any,Any}",
    "page": "API",
    "title": "FEMBasis.grad",
    "category": "Method",
    "text": "grad(B, X, xi)\n\nGiven basis B, calculate gradient dB/dX at xi.\n\nExample\n\nB = Quad4()\nX = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])\ngrad(B, X, (0.0, 0.0))\n\n# output\n\n2×4 Array{Float64,2}:\n -0.5   0.5  0.5  -0.5\n -0.5  -0.5  0.5   0.5\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.interpolate-Tuple{Any,Any,Any}",
    "page": "API",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(B, T, xi)\n\nGiven basis B, interpolate T at xi.\n\nExample\n\nB = Quad4()\nX = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))\nT = (1.0, 2.0, 3.0, 4.0)\ninterpolate(B, T, (0.0, 0.0))\n\n# output\n\n2.5\n\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.calculate_basis_coefficients-Tuple{Expr,Tuple}",
    "page": "API",
    "title": "FEMBasis.calculate_basis_coefficients",
    "category": "Method",
    "text": "function calculate_basis_coefficients(polynomial::String, coordinates::Vararg{Tuple})\n\nCalculate \"interpolate coefficient matrix\" for some polynomial p.\n\nExamples\n\nThat is, if we have polynomial p = 1 + u + v and coordinates (0,0), (1,0), (0,1), we find A*p such that first row is the first coordinate, second row is second coordinate and so on:\n\njulia> p = \"1 + u + w\"\njulia> X = ((0.0,0.0), (1.0,0.0), (0.0,1.0))\njulia> calculate_basis_coefficient(p, X)\n[1.0 0.0 0.0 # <-- p(0.0,0.0) = 1.0     = [1.0 0.0 0.0] * [1.0, u, v]\n 1.0 1.0 0.0 # <-- p(1.0,0.0) = 1.0 + u = [1.0 1.0 0.0] * [1.0, u, v]\n 1.0 0.0 1.0] # <- p(0.0,1.0) = 1.0 + v = [1.0 0.0 1.0] * [1.0, u, v]\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.calculate_interpolation_polynomial_derivatives-Tuple{Array{T,1} where T,Int64}",
    "page": "API",
    "title": "FEMBasis.calculate_interpolation_polynomial_derivatives",
    "category": "Method",
    "text": "Calculate derivatives of basis functions with respect to parameters u, v, w.\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.calculate_interpolation_polynomials-Tuple{Expr,Array{T,2} where T}",
    "page": "API",
    "title": "FEMBasis.calculate_interpolation_polynomials",
    "category": "Method",
    "text": "\n\n"
},

{
    "location": "api.html#FEMBasis.jacobian-Tuple{Any,Any,Any}",
    "page": "API",
    "title": "FEMBasis.jacobian",
    "category": "Method",
    "text": "jacobian(B, X, xi)\n\nGiven basis B, calculate jacobian at xi.\n\nExample\n\nB = Quad4()\nX = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])\njacobian(B, X, (0.0, 0.0))\n\n# output\n\n2×2 Array{Float64,2}:\n 0.5  0.0\n 0.0  0.5\n\n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLHEX1}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 1 point rule on hexahedron. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLHEX243}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 243 point rule on quadrilateral. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLHEX27}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 27 point rule on hexahedron. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLHEX81}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 81 point rule on hexahedron. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLHEX8}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 8 point rule on hexahedron. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLPYR5B}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 5 point rule on pyramid. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLPYR5}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 5 point rule on pyramid. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLQUAD16}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 16 point rule on quadrilateral. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLQUAD1}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 1 point rule on quadrilateral. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLQUAD25}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 25 point rule on quadrilateral. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLQUAD4}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 4 point rule on quadrilateral. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLQUAD9}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 9 point rule on quadrilateral. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLSEG1}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 1 point rule on segment. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLSEG2}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 2 point rule on segment. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLSEG3}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 3 point rule on segment. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLSEG4}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 4 point rule on segment. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLSEG5}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 5 point rule on segment. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTET15}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 15 point rule on tetrahedron. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTET1}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 1 point rule on tetrahedron. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTET4}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 4 point rule on tetrahedron. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTET5}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 5 point rule on tetrahedron. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTRI12}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 12 point rule on triangle. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTRI1}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 1 point rule on triangle. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTRI3B}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 3 point rule on triangle. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTRI3}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 3 point rule on triangle. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTRI4B}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 4 point rule on triangle. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTRI4}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 4 point rule on triangle. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTRI6}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 6 point rule on triangle. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLTRI7}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 7 point rule on triangle. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLWED21}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 21 point rule on wedge. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLWED6B}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 6 point rule on wedge. \n\n\n\n"
},

{
    "location": "api.html#FEMQuad.get_quadrature_points-Tuple{Type{Val{:GLWED6}}}",
    "page": "API",
    "title": "FEMQuad.get_quadrature_points",
    "category": "Method",
    "text": "Gauss-Legendre quadrature, 6 point rule on wedge. \n\n\n\n"
},

{
    "location": "api.html#Index-1",
    "page": "API",
    "title": "Index",
    "category": "section",
    "text": "DocTestSetup = quote\n    using FEMBase\n    using FEMBasis\n    using FEMQuad\nendModules = [FEMBase, FEMBasis, FEMQuad]"
},

]}
