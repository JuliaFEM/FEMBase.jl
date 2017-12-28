var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#FEMBase.jl-1",
    "page": "Introduction",
    "title": "FEMBase.jl",
    "category": "section",
    "text": "(Image: Build Status)(Image: Coverage Status)(Image: )(Image: )(Image: Issues)FEMBase.jl is a JuliaFEM base package. It includes all basic data structures so that one can start programming own finite element models."
},

{
    "location": "fields.html#",
    "page": "Fields",
    "title": "Fields",
    "category": "page",
    "text": "using FEMBaseFrom the beginning a clear concept: \"everything is a field\". We think that constant is just a special case of field which does not vary in temporal and spatial direction. Fields can vary in spatial direction (constant or variable) and in temporal direction (time variant and in-variant).From here we can deduct for types of (discrete) fields:discrete, constant, time invariant (DCTI)\ndiscrete, variable, time invariant (DVTI)\ndiscrete, constant, time variant (DCTV)\ndiscrete, variable, time variant (DVTV)Field itself can be anything. However, usually either scalar, vector or tensor (matrix)."
},

{
    "location": "fields.html#Creating-fields-1",
    "page": "Fields",
    "title": "Creating fields",
    "category": "section",
    "text": "For discrete fields that are varying in spatial direction, value for each discrete point is defined using NTuple. The order of points is implicitly assumed to be same than node ordering in ABAQUS. Usual use case for variable field is that the field is interpolated to the element volume using interpolation polynomials, i.e. u(xi) = u_i N_i(xi), where N_i is the basis function for that node and u_i is the discrete value.For example, (1, 2, 3, 4) is a scalar field having length of 4 and ([1,2],[2,3],[3,4],[4,5]) is a vector field having length of 4.For fields that are varying in temporal direction, time => value syntax is used. The first item in pair is time and second item is value attached to that time. For example, 0.0 => 1.0 is a time-dependent scalar field having value 1.0 at time 0.0.The most simple field is a field that is constant in both time and spatial direction. Discrete, constant, time invariant (DCTI):a = DCTI(1.0)FEMBase.DCTI{Float64}(1.0)Then we have discrete, variable, time invariant fields (DVTI). Note the use of NTuple when defining field.b = DVTI( (1.0, 2.0) )FEMBase.DVTI{2,Float64}((1.0, 2.0))Discrete, constant, time variant field (DCTV) is constant in spatial direction but can vary in temporal direction. Here, => syntax is used. New values can be added to field using update!. If there already exists value for that particular time, it will be overridden. It is assumed that content of field in time direction is monotonically increasing, i.e. t_i-1  t_i  t_i+1. For the sake of clarity let's also mention that update! works for time invariant fields as well if content needs to be updated.c = DCTV(0.0 => 1.0, 1.0 => 2.0)\nupdate!(c, 2.0 => 3.0)\ncFEMBase.DCTV{Float64}(Pair{Float64,Float64}[0.0=>1.0, 1.0=>2.0, 2.0=>3.0])Discrete, variable, time variant (DVTV) field is the most general one, allowing values of field to vary in both spatial and time direction.d = DVTV(0.0 => (1.0, 2.0), 1.0 => (2.0, 3.0))\nupdate!(d, 2.0 => (3.0, 4.0))3-element Array{Pair{Float64,Tuple{Float64,Float64}},1}:\n 0.0=>(1.0, 2.0)\n 1.0=>(2.0, 3.0)\n 2.0=>(3.0, 4.0)In examples above, all fields was scalar fields. Defining vector or tensor fields goes in the same spirit. Only difference is that now we define vectors and tensors, not a single scalar value. They can vary in spatial and time direction in the same way than scalar fields. Here is example of defining the abovementioned vector fields:a = DCTI([1.0, 2.0])\nb = DVTI(([1.0, 2.0], [2.0, 3.0]))\nc = DCTV(0.0 => [1.0, 2.0], 1.0 => [2.0, 3.0])\nd = DVTV(0.0 => ([1.0, 2.0], [2.0, 3.0]), 1.0 => ([2.0, 3.0], [3.0, 4.0]))FEMBase.DVTV{2,Array{Float64,1}}(Pair{Float64,Tuple{Array{Float64,1},Array{Float64,1}}}[0.0=>([1.0, 2.0], [2.0, 3.0]), 1.0=>([2.0, 3.0], [3.0, 4.0])])"
},

{
    "location": "fields.html#Accessing-fields-1",
    "page": "Fields",
    "title": "Accessing fields",
    "category": "section",
    "text": "Accessing fields is done using a single command: interpolate. For time varying fields, one can interpolate in time direction. For example, if we have (constant) 12 at time t=00 and 34 at time t=10, linear interpolation yieldsc = DCTV(0.0 => [1.0,2.0], 1.0 => [3.0,4.0])\ninterpolate(c, 0.5)2-element Array{Float64,1}:\n 2.0\n 3.0If field is spatially varying, a Tuple will be returned, typically having one value for each \"node\". This can then be interpolated in spatial direction, typically using basis functions defined on element, e.g. u = N_i u_i:d = DVTV(0.0 => (1.0,2.0), 1.0 => (3.0,4.0))\ninterpolate(d, 0.5)(2.0, 3.0)Although the two fields above looks quite same, the key difference is that in DCTV field we have a constant vectorial value (defined using square brackets []) and in latter DVTV field we have a scalar value (defined using round brackets) changing in spatial direction from 1.0 to 2.0 at time t=00 and changing from 3.0 to 4.0 at time t=10.One should be always able to interpolate in time direction, even if field is time invariant, to get trivial solution:interpolate(a, 0.5), interpolate(b, 0.5),\ninterpolate(c, 0.5), interpolate(d, 0.5)([1.0, 2.0], ([1.0, 2.0], [2.0, 3.0]), [2.0, 3.0], (2.0, 3.0))For spatially varying fields, one can access to ith element using getindex:getindex(a, 1), getindex(b, 1), getindex(b, 2)([1.0, 2.0], [1.0, 2.0], [2.0, 3.0])First time interpolation, then spatial lookup, i.e.getindex(interpolate(d, 0.5), 1)2.0By passing already evaluated basis functions to interpolate, field is intepolated in both temporal and spatial direction.d = DVTV(0.0 => (1.0,2.0), 1.0 => (3.0,4.0))\nxi = 0.0\nN = [1/2*(1-xi) 1/2*(1+xi)] # evaluated basis functions, linear interpolation\nt = 1.0 # time\ninterpolate(d, t, N)3.5Internally, each field is a subtype of AbstractField. Internally each field has a data which be accessed directly for hacking.d.data2-element Array{Pair{Float64,Tuple{Float64,Float64}},1}:\n 0.0=>(1.0, 2.0)\n 1.0=>(3.0, 4.0)"
},

{
    "location": "fields.html#Continuous-fields-1",
    "page": "Fields",
    "title": "Continuous fields",
    "category": "section",
    "text": "Then we have continuous fields which may be useful when defining analytical boundary conditions. For that we have CVTV, where \"C\" stands for continuous.f(xi,t) = xi[1]*xi[2]*t\ng = CVTV(f)\ng([1.0,2.0],3.0)6.0"
},

{
    "location": "fields.html#Dictionary-fields-1",
    "page": "Fields",
    "title": "Dictionary fields",
    "category": "section",
    "text": "Usually it is assumed that discrete field has values number of nodes or degrees of freedom in element, typically something small like 1-30. There might be cases where it is more practical to define field in a sense that indexing is not continuous or starting from 1. For example, we might want to define field common for a set of elements. In that case it's natural to think that each index in field corresponds to the certain id-number of node. For example, if we have triangle element connecting nodes 1, 1000 and 100000, we still want to access to that field naturally using getindex, e.g. f[1], f[1000] and f[100000]. In that case, more appropriate internal structure for field is based on dictionary.It only makes sense to define dictionary fields for spatially varying fields. Two new fields are introduced: DVTId and DVTVd, where last \"d\" stands for \"dictionary\".Keep on mind, that this type of field has one restriction. If and when this field is typically defined on nodes of several elements, field must be continuous between elements. That is, if field value in node 1000 is for example 1.0, then it's 1.0 in all elements connecting to that node. To define jumps on field, one must define field element-wise. # Define eg. \"geometry\" for nodes 1,2,3,4.\nX = Dict(1=>[0.0,0.0], 1000=>[1.0,0.0], 100000=>[1.0,1.0])\nG = DVTId(X)\nG[1], G[1000], G[100000]([0.0, 0.0], [1.0, 0.0], [1.0, 1.0])Y = Dict(1=>[1.0,1.0], 1000=>[2.0,1.0], 100000=>[2.0,2.0])\nF = DVTVd(0.0 => X, 1.0 => Y)\ninterpolate(F,0.5)[100000]2-element Array{Float64,1}:\n 1.5\n 1.5"
},

{
    "location": "fields.html#Using-common-constructor-field-1",
    "page": "Fields",
    "title": "Using common constructor field",
    "category": "section",
    "text": "Now we have introduced total of 7 fields: DCTI, DCTV, DVTI, DVTV, CVTV, DVTId, DVTVd. A good question arises that how to remember all this stuff and is it even necessary? Luckily not, because one can use a single constructor called field to create all kind of fields. Type of field is inspected from data type. You really don't have to know about this technical stuff, just declare new field using intuition and field command.f1 = field(1.0)\ntypeof(f1)FEMBase.DCTI{Float64}f2 = field((1.0, 2.0))\ntypeof(f2)FEMBase.DVTI{2,Float64}f3 = field(0.0 => 1.0)\ntypeof(f3)FEMBase.DCTV{Float64}f4 = field(0.0 => (1.0, 2.0), 1.0 => (2.0, 3.0))\ntypeof(f4)FEMBase.DVTV{2,Float64}f5 = field((xi,t) -> xi[1]*t)\ntypeof(f5)FEMBase.CVTVf6 = field(1 => 1.0, 2 => 2.0)\ntypeof(f6)FEMBase.DVTId{Float64}f7 = field(0.0 => (1=>1.0,10=>2.0), 1.0 => (1=>2.0,10=>3.0))\ntypeof(f7)FEMBase.DVTVd{Float64}"
},

{
    "location": "fields.html#Developing-new-fields-1",
    "page": "Fields",
    "title": "Developing new fields",
    "category": "section",
    "text": "If the default ones are not enough, it's always possible to define new ones. Minimum requirements is that field is a subtype of AbstractField and interpolate has been defined to it."
},

{
    "location": "basis.html#",
    "page": "Basis functions",
    "title": "Basis functions",
    "category": "page",
    "text": ""
},

{
    "location": "basis.html#Shape-functions-1",
    "page": "Basis functions",
    "title": "Shape functions",
    "category": "section",
    "text": "Also known as basis functions, interpolation polynomials and so on. Typically unknown field variable is interpolated from element nodal values using continuous functions. That is, Standard Lagrange polynomials as supported.Linear shape functions:Seg2\nTri3\nQuad4\nTet4\nPyr5\nWedge6\nHex8Quadratic and biquadratic shape functions:Seg3\nTri6, Tri7\nQuad8, Quad9\nTet10\nWedge15\nHex20, Hex27NURBS shape functions:NSeg\nNSurf\nNSolidEvaluating basis and derivative of basis functions with respect to dimensionless coordinates:using FEMBase\nusing FEMBase.FEMBasis: eval_dbasis!, jacobian, grad, interpolate,\n                        get_reference_element_coordinates, create_basis\nB = Quad4()FEMBasis.Quad4()length(B)4size(B)(2, 4)For fast evaluations, one must allocate array outside of the hot loops to get speed.N = zeros(1, 4)\ndN = zeros(2, 4)\nxi = (0.0, 0.0)(0.0, 0.0)eval_basis!(B, N, xi)1×4 Array{Float64,2}:\n 0.25  0.25  0.25  0.25eval_dbasis!(B, dN, xi)2×4 Array{Float64,2}:\n -0.25   0.25  0.25  -0.25\n -0.25  -0.25  0.25   0.25For Langrange interpolation polynomials, by definition, on each node shape function corresponding to that node gets value of 1 and the rest is zero. Node ordering follows the same defined in e.g. in ABAQUS and in many other FEM softwares.get_reference_element_coordinates(Quad4)((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0))for xi in get_reference_element_coordinates(Quad4)\n    eval_basis!(B, N, xi)\n    println(\"$N at $xi\")\nend[1.0 0.0 0.0 0.0] at (-1.0, -1.0)\n[0.0 1.0 0.0 0.0] at (1.0, -1.0)\n[0.0 0.0 1.0 0.0] at (1.0, 1.0)\n[0.0 0.0 0.0 1.0] at (-1.0, 1.0)"
},

{
    "location": "basis.html#Mathematics-1",
    "page": "Basis functions",
    "title": "Mathematics",
    "category": "section",
    "text": "Without knowing anything about the real shape of domain, eval_dbasis! is calculating gradient with respect to dimensionless coordinates xi_i, i.e.\\begin{equation} \\frac{\\partial\\boldsymbol{N}}{\\partial\\boldsymbol{\\xi}}=\\left[\\begin{array}{cccc} \\frac{\\partial N_{1}}{\\partial\\xi_{1}} & \\frac{\\partial N_{2}}{\\partial\\xi_{1}} & \\cdots & \\frac{\\partial N_{n}}{\\partial\\xi_{1}}\\\n\\frac{\\partial N_{1}}{\\partial\\xi_{2}} & \\frac{\\partial N_{2}}{\\partial\\xi_{2}} & \\cdots & \\frac{\\partial N_{n}}{\\partial\\xi_{2}}\\\n\\frac{\\partial N_{1}}{\\partial\\xi_{3}} & \\frac{\\partial N_{2}}{\\partial\\xi_{3}} & \\cdots & \\frac{\\partial N_{n}}{\\partial\\xi_{3}} \\end{array}\\right] \\end{equation}Usually this is not wanted, but instead gradient of basis functions is calculated with respect to natural coordinates X_i,\\begin{equation} \\frac{\\partial\\boldsymbol{N}}{\\partial\\boldsymbol{X}}=\\left[\\begin{array}{cccc} \\frac{\\partial N_{1}}{\\partial X_{1}} & \\frac{\\partial N_{2}}{\\partial X_{1}} & \\cdots & \\frac{\\partial N_{n}}{\\partial X_{1}}\\\n\\frac{\\partial N_{1}}{\\partial X_{2}} & \\frac{\\partial N_{2}}{\\partial X_{2}} & \\cdots & \\frac{\\partial N_{n}}{\\partial X_{2}}\\\n\\frac{\\partial N_{1}}{\\partial X_{3}} & \\frac{\\partial N_{2}}{\\partial X_{3}} & \\cdots & \\frac{\\partial N_{n}}{\\partial X_{3}} \\end{array}\\right] \\end{equation}To get this, inverse of Jacobian matrix is needed.X = ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])\nxi = (0.0, 0.0)\nJ = jacobian(B, X, xi)2×2 Array{Float64,2}:\n 0.5  0.0\n 0.0  0.5inv(J) * dN2×4 Array{Float64,2}:\n -0.5   0.5  0.5  -0.5\n -0.5  -0.5  0.5   0.5Or directly:dNdX = grad(B, X, xi)2×4 Array{Float64,2}:\n -0.5   0.5  0.5  -0.5\n -0.5  -0.5  0.5   0.5If interpolation domain is manifold, Jacobian is not square and inverse cannot be taken.X2 = ([0.0,0.0,0.0], [1.0, 0.0,1.0], [1.0,1.0,1.0], [0.0,1.0,0.0])\nxi = (0.0, 0.0)\nJ = jacobian(B, X2, xi)2×3 Array{Float64,2}:\n 0.5  0.0  0.5\n 0.0  0.5  0.0One can use Jacobian to calculate surface integral:\\begin{equation} \\iint_{S}f\\,\\mathrm{d}\\Sigma=\\iint_{T}f\\left(\\boldsymbol{x}\\left(s,t\\right)\\right)\\left\\Vert \\frac{\\partial\\boldsymbol{x}}{\\partial s}\\times\\frac{\\partial\\boldsymbol{x}}{\\partial t}\\right\\Vert \\,\\mathrm{d}s\\mathrm{d}t \\end{equation}4*norm(cross(J[1,:], J[2,:])), sqrt(2) # area of manifold(1.4142135623730951, 1.4142135623730951)Gradient of e.g. displacement field or temperature field can be also evaluated:u = ([0.0, 0.0], [1.0, -1.0], [2.0, 3.0], [0.0, 0.0])\nT = (1.0, 2.0, 3.0, 4.0)\ngrad(B, u, X, xi)2×2 Array{Float64,2}:\n 1.5  0.5\n 1.0  2.0grad(B, T, X, xi)1×2 RowVector{Float64,Array{Float64,1}}:\n 0.0  2.0One can interpolate fields using basis:interpolate(B, u, xi)2-element Array{Float64,1}:\n 0.75\n 0.5interpolate(B, T, xi)2.5At last to avoid unnecessary memory allocations, a struct BasisInfo is introduced, containing memory space for calculations.bi = BasisInfo(Quad4)\neval_basis!(bi, X, xi)FEMBasis.BasisInfo{FEMBasis.Quad4,Float64}([0.25 0.25 0.25 0.25], [-0.25 0.25 0.25 -0.25; -0.25 -0.25 0.25 0.25], [-0.5 0.5 0.5 -0.5; -0.5 -0.5 0.5 0.5], [0.5 0.0; 0.0 0.5], [2.0 -0.0; -0.0 2.0], 0.25)bi.J2×2 Array{Float64,2}:\n 0.5  0.0\n 0.0  0.5bi.N1×4 Array{Float64,2}:\n 0.25  0.25  0.25  0.25bi.dN2×4 Array{Float64,2}:\n -0.25   0.25  0.25  -0.25\n -0.25  -0.25  0.25   0.25bi.detJ0.25bi.grad2×4 Array{Float64,2}:\n -0.5   0.5  0.5  -0.5\n -0.5  -0.5  0.5   0.5bi.invJ2×2 Array{Float64,2}:\n  2.0  -0.0\n -0.0   2.0gradu = zeros(2, 2)\ngrad!(bi, gradu, u)2×2 Array{Float64,2}:\n 1.5  0.5\n 1.0  2.0"
},

{
    "location": "basis.html#Defining-custom-shape-functions-1",
    "page": "Basis functions",
    "title": "Defining custom shape functions",
    "category": "section",
    "text": "Depending from the type of shape functions, they can be created more or less automatic way. For Lagrange type interpolation, ones needs only to define polynomial and corner points for domain. For example, if domain is 01^2, one can use create_basis:code = create_basis(\n    :MyQuad,\n    \"My special domain\",\n    (\n        (0.0, 0.0),\n        (1.0, 0.0),\n        (1.0, 1.0),\n        (0.0, 1.0),\n    ),\n    \"1 + u + v + u*v\"\n)\neval(code)B = MyQuad()\nxi = (0.5, 0.5)\neval_basis!(B, N, xi)1×4 Array{Float64,2}:\n 0.25  0.25  0.25  0.25In this case partial derivatives of shape functions are with respect to X, because interpolation polynomials are calculated against real domain and not \"reference domain\":eval_dbasis!(B, dN, xi)2×4 Array{Float64,2}:\n -0.5   0.5  0.5  -0.5\n -0.5  -0.5  0.5   0.5J = jacobian(B, X, xi)2×2 Array{Float64,2}:\n 1.0  0.0\n 0.0  1.0u = ([0.0, 0.0], [1.0, -1.0], [2.0, 3.0], [0.0, 0.0])\ngrad(B, u, X, xi)2×2 Array{Float64,2}:\n 1.5  0.5\n 1.0  2.0Shape functions can be defined manually and calculate partial derivatives automatically. C1-continuous Hermite shape functions can be defined as:code = create_basis(\n    :C1Hermite,\n    \"C1-continuous Hermite shape functions\",\n    (\n        (0.0,),\n        (0.0,),\n        (1.0,),\n        (1.0,)\n    ),\n    (\n        \"2*u^3 - 3*u^2 + 1\",\n        \"u^3 - 2*u^2 + u\",\n        \"-2*u^3 + 3*u^2\",\n        \"u^3 - u^2\"\n    )\n)\neval(code)B = C1Hermite()\nxi = (0.0,)\neval_basis!(B, N, xi)1×4 Array{Float64,2}:\n 1.0  0.0  0.0  0.0dN = zeros(1, 4)\neval_dbasis!(B, dN, xi)1×4 Array{Float64,2}:\n 0.0  1.0  0.0  0.0xi = (1.0,)\neval_basis!(B, N, xi)1×4 Array{Float64,2}:\n 0.0  0.0  1.0  0.0eval_dbasis!(B, dN, xi)1×4 Array{Float64,2}:\n 0.0  0.0  0.0  1.0The last option is to build everything from scratch. For that, one must import and define following functions:Base.size\nBase.length\nFEMBase.FEMBasis.get_reference_element_coordinates\nFEMBase.FEMBasis.eval_basis!\nFEMBase.FEMBasis.eval_dbasis!A simple implementation of P-hierarchical 1d-basis would then beimport Base: size, length\nimport FEMBase: get_reference_element_coordinates,\n                eval_basis!, eval_dbasis!,\n                AbstractBasis\n\ntype PSeg <: AbstractBasis\n    order :: Int\nend\n\nfunction PSeg()\n    return PSeg(1)\nend\n\nfunction length(basis::PSeg)\n    return basis.order+1\nend\n\nfunction size(basis::PSeg)\n    return (1, basis.order+1)\nend\n\nfunction get_reference_element_coordinates(basis::PSeg)\n    return ((-1.0,), (1.0,))\nend\n\n\"\"\"\n    get_legendre_polynomial(n)\n\nReturn Legendgre polynomial of order `n` to inverval ξ ∈ [1, 1].\n\nImplementation uses Bonnet's recursion formula. See\nhttps://en.wikipedia.org/wiki/Legendre_polynomials\n\"\"\"\nfunction get_legendre_polynomial(n)\n    n == 0 && return xi -> 1.0\n    n == 1 && return xi -> xi\n    Pm1 = get_legendre_polynomial(n-1)\n    Pm2 = get_legendre_polynomial(n-2)\n    A(xi) = (2.0*n-1.0)*xi*Pm1(xi)\n    B(xi) = (n-1.0)*Pm2(xi)\n    return xi -> (A(xi)-B(xi))/n\nend\n\n\"\"\"\n    get_legendre_polynomial_derivative(n)\n\nReturn derivative of Legendgre polynomial of order `n` to\ninverval ξ ∈  [-1, 1]\n\"\"\"\nfunction get_legendre_polynomial_derivative(n)\n    n == 0 && return xi -> 0.0\n    n == 1 && return xi -> 1.0\n    Pm1 = get_legendre_polynomial_derivative(n-1)\n    Pm2 = get_legendre_polynomial_derivative(n-2)\n    A(xi) = (2.0*(n-1.0)+1.0)*xi*Pm1(xi)\n    B(xi) = (n+1.0-1.0)*Pm2(xi)\n    return xi -> (A(xi)-B(xi))/(n-1.0)\nend\n\nfunction eval_basis!{T}(basis::PSeg, N::Matrix{T}, xi::Tuple{T})\n    n = length(basis)\n    t = xi[1]\n    N[1] = 0.5*(1-t)\n    N[2] = 0.5*(1+t)\n    n < 3 && return N\n    for i=3:n\n        j = i-1\n        P1 = get_legendre_polynomial(j)\n        P2 = get_legendre_polynomial(j-2)\n        N[i] = 1.0/sqrt(2.0*(2.0*j-1.0))*(P1(t)-P2(t))\n    end\n    return N\nend\n\nfunction eval_dbasis!{T}(basis::PSeg, dN::Matrix{T}, xi::Tuple{T})\n    n = length(basis)\n    t = xi[1]\n    dN[1] = -0.5\n    dN[2] = 0.5\n    n < 3 && return N\n    for i=3:n\n        j = i-1\n        P1 = get_legendre_polynomial_derivative(j)\n        P2 = get_legendre_polynomial_derivative(j-2)\n        dN[i] = 1.0/sqrt(2.0*(2.0*j-1.0))*(P1(t)-P2(t))\n    end\n    return dN\nendeval_dbasis! (generic function with 21 methods)B = PSeg()PSeg(1)N = zeros(1, 2)\neval_basis!(B, N, (0.0,))1×2 Array{Float64,2}:\n 0.5  0.5N = zeros(1, 3)\nB.order = 2\neval_basis!(B, N, (0.0,))1×3 Array{Float64,2}:\n 0.5  0.5  -0.612372using PyPlot\nB.order = 6\nN = zeros(1, length(B))\nn = 50\nxi = linspace(-1, 1, n)\nNN = zeros(n, length(B))\nfor i=1:n\n    eval_basis!(B, N, (xi[i],))\n    NN[i,:] = N[:]\nendplot(NN)\ntitle(\"Hierarchical shape functions to order 6\");(Image: png)"
},

{
    "location": "elements.html#",
    "page": "Elements",
    "title": "Elements",
    "category": "page",
    "text": "In JuliaFEM, elements are \"containers\", combining fields and basis functions described above. Among that, element has information about topology (connectivity) and integration rule. These fundamentals forms a finite element, the backbone of finite element method, as the basic idea after all is just to discretize continuous domain to smaller topological entities like tetrahedrons and perform same operations to each element.using FEMBaseel = Element(Quad4, [1, 2, 3, 4])FEMBase.Element{FEMBasis.Quad4}(-1, [1, 2, 3, 4], FEMBase.Point{FEMBase.IntegrationPoint}[], Dict{String,FEMBase.AbstractField}(), FEMBasis.Quad4())Setting fields to element is done using a command update!, which either creates new field if does not already exist, or updates the old one. Typically, at least field called geometry needs to be defined to element as it's used to calculate Jacobian of element. Fields can be discrete, continuous, time invariant or variant, variable or constant, like described earlier.X = Dict(1 => [0.0,0.0], 2=>[1.0,0.0], 3=>[1.0,1.0], 4=>[0.0,1.0])\nupdate!(el, \"geometry\", X)FEMBase.DVTId{Array{Float64,1}}(Dict(4=>[0.0, 1.0],2=>[1.0, 0.0],3=>[1.0, 1.0],1=>[0.0, 0.0]))el.fieldsDict{String,FEMBase.AbstractField} with 1 entry:\n  \"geometry\" => FEMBase.DVTId{Array{Float64,1}}(Dict(4=>[0.0, 1.0],2=>[1.0, 0.0…u0 = ([0.0,0.0], [0.0,0.0], [0.0,0.0], [0.0,0.0])\nu1 = ([0.0,0.0], [0.0,0.0], [0.5,0.0], [0.0,0.0])\nupdate!(el, \"displacement\", 0.0 => u0)\nupdate!(el, \"displacement\", 1.0 => u1)\nel.fieldsDict{String,FEMBase.AbstractField} with 2 entries:\n  \"geometry\"     => FEMBase.DVTId{Array{Float64,1}}(Dict(4=>[0.0, 1.0],2=>[1.0,…\n  \"displacement\" => FEMBase.DVTV{4,Array{Float64,1}}(Pair{Float64,NTuple{4,Arra…Interpolating of fields goes calling Element(field_name, xi, time). For example, position of material particle X in initial configuration and deformed configuration in the middle of the element at time t=1 can be found asxi = (0.0, 0.0)\ntime = 1.0\nX = el(\"geometry\", xi, time)\nu = el(\"displacement\", xi, time)\nx = X+u\nprintln(\"X = $X, x = $x\")X = [0.5, 0.5], x = [0.625, 0.5]Jacobian, determinant of Jacobian and gradient of field can be calculated adding extra argument Val{:Jacobian}, Val{:detJ}, Val{:Grad} to the above command and not passing field name, i.e.el(xi, time, Val{:Jacobian})2×2 Array{Float64,2}:\n 0.5  0.0\n 0.0  0.5el(xi, time, Val{:detJ})0.25el(xi, time, Val{:Grad})2×4 Array{Float64,2}:\n -0.5   0.5  0.5  -0.5\n -0.5  -0.5  0.5   0.5Usually what the user wants is still a gradient of some field. For example, displacement gradient:gradu = el(\"displacement\", xi, time, Val{:Grad})\ngradu2×2 Array{Float64,2}:\n 0.25  0.25\n 0.0   0.0Or temperature gradient:update!(el, \"temperature\", (1.0, 2.0, 3.0, 4.0))\ngradT = el(\"temperature\", xi, time, Val{:Grad})1×2 RowVector{Float64,Array{Float64,1}}:\n 0.0  2.0Accessing integration points of element is done using command get_integration_points. Combining interpolation and integration one can already calculate local matrices of a single element or, for example area and strain energy:update!(el, \"lambda\", 96.0)\nupdate!(el, \"mu\", 48.0)\n\nA = 0.0\nW = 0.0\nfor ip in get_integration_points(el)\n    detJ = el(ip, time, Val{:detJ})\n    A += ip.weight * detJ\n    ∇u = el(\"displacement\", ip, time, Val{:Grad})\n    E = 1/2*(∇u + ∇u')\n    λ = el(\"lambda\", ip, time)\n    μ = el(\"mu\", ip, time)\n    W += ip.weight * ( λ/2*trace(E*E') + μ*trace(E)^2) * detJ\nend\n\nprintln(\"Area: $A\")\nprintln(\"Strain energy: $W\")Area: 1.0\nStrain energy: 10.0Local stiffness matrix for Poisson problem:K = zeros(4,4)\nupdate!(el, \"coefficient\", 36.0)\nfor ip in get_integration_points(el)\n    dN = el(ip, time, Val{:Grad})\n    detJ = el(ip, time, Val{:detJ})\n    c = el(\"coefficient\", ip, time)\n    K += ip.weight * c*dN'*dN * detJ\nend\nK4×4 Array{Float64,2}:\n  24.0   -6.0  -12.0   -6.0\n  -6.0   24.0   -6.0  -12.0\n -12.0   -6.0   24.0   -6.0\n  -6.0  -12.0   -6.0   24.0"
},

{
    "location": "developing.html#",
    "page": "Developing",
    "title": "Developing",
    "category": "page",
    "text": ""
},

{
    "location": "developing.html#Developing-new-input/output-interfaces-1",
    "page": "Developing",
    "title": "Developing new input/output interfaces",
    "category": "section",
    "text": "E.g. mesh reader or results writer.type MeshReader <: AbstractMeshReader\nend\n\n\"\"\"\n    read_mesh(mesh::MeshReader, data)\n\nReads mesh from disk/memory/cloud/sql/etc. and returns Mesh object\n\"\"\"\nfunction read_mesh(mesh::MeshReader, filename::String)\n    # do something ...\nendtype ResultsWriter <: AbstractResultsWriter\nend\n\n\"\"\"\n    write_results!(results::ResultsWriter, data)\n\nGiven data, write calculation results back to disk/memory/cloud/sql/etc.\n\"\"\"\nfunction write_results!(results::ResultsWriter, data)\n    # write results ...\nend"
},

{
    "location": "developing.html#Developing-new-physical-models-1",
    "page": "Developing",
    "title": "Developing new physical models",
    "category": "section",
    "text": "Starting point is a weak formulation, which is then discretized to elements.type Heat <: AbstractProblem\nend\n\n\"\"\"\n    assemble!(settings::Heat, problem::Problem, assembly::Assembly, elements::Vector{Element}, time::Float64)\n\nGiven `problem` spesification, assemble `elements` to global stiffness matrix and force\nvector defined in `assembly` for some given `time`.\n\"\"\"\nfunction assemble!(settings::Heat, problem::Problem, assembly::Assembly,\n                   elements::Vector{Element}, time::Float64)\n    # integrate local matrices and add them to assembly\nend"
},

{
    "location": "developing.html#Developing-new-material-model-1",
    "page": "Developing",
    "title": "Developing new material model",
    "category": "section",
    "text": "Aim is to define material response given data.type LinearIsotropic <: AbstractMaterial\nend\n\n\"\"\"\n    calculate_material_response!(material::LinearIsotropic, data)\n\nGiven strain tensor and some other quantities, calculate material response\n\"\"\"\nfunction calculate_response!(material::LinearIsotropic, data)\n    # given strain tensor and some other quantities, calculate stress\nend"
},

{
    "location": "developing.html#Developing-new-interpolation-functions-1",
    "page": "Developing",
    "title": "Developing new interpolation functions",
    "category": "section",
    "text": "E.g. basis functionstype LinQuad4Basis <: AbstractBasis\nend\n\n\"\"\"\n    evaluate_basis!(basis::LinQuad4Basis, element::Element, xi, time, N::Matrix{Float64})\n\nEvaluate basis functions at some point $\\xi$ and store results to `N`.\n\"\"\"\nfunction evaluate_basis!(basis::LinQuad4Basis, element::Element, xi, time, N::Matrix{Float64})\n    # populate N with new basis\nend"
},

{
    "location": "developing.html#Developing-new-integration-rules-1",
    "page": "Developing",
    "title": "Developing new integration rules",
    "category": "section",
    "text": "type Quad4PointGaussLegendre <: AbstractIntegrationRule\nend\n\n\"\"\"\n    get_integration_points(q::Quad4PointQaussLegendre)\n\nReturn integration point locations and weights.\n\"\"\"\nfunction get_integration_points(q::Quad4PointGaussLegendre)\n    # return integration points\nend"
},

{
    "location": "developing.html#Developing-new-solver-1",
    "page": "Developing",
    "title": "Developing new solver",
    "category": "section",
    "text": "type ImplicitTimeSolver <: AbstractSolver\nend\n\n\"\"\"\n    solve!(settings::ImplicitTimeSolver, solver::Solver, time)\n\nAssemble problems, solve problem, update problems, write results and so on.\n\"\"\"\nfunction solve!(settings::ImplicitTimeSolver, solver::Solver, time)\n    # do solution\nend"
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
    "location": "api.html#FEMBase.CVTV",
    "page": "API",
    "title": "FEMBase.CVTV",
    "category": "Type",
    "text": "CVTV(f)\n\nContinuous, variable, time variant field.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.DCTI",
    "page": "API",
    "title": "FEMBase.DCTI",
    "category": "Type",
    "text": "DCTI(T)\n\nDiscrete, constant, time-invariant field. This field is constant in both spatial direction and time direction, i.e. df/dX = 0 and df/dt = 0.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.DVTI",
    "page": "API",
    "title": "FEMBase.DVTI",
    "category": "Type",
    "text": "Discrete, variable, time-invariant field. This is constant in time direction, but not in spatial direction, i.e. df/dt = 0 but df/dX != 0. The basic structure of data is Vector, and it is implicitly assumed that length of field matches to the number of shape functions, so that interpolation in spatial direction works.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.DVTId",
    "page": "API",
    "title": "FEMBase.DVTId",
    "category": "Type",
    "text": "DVTId(X::Dict)\n\nDiscrete, variable, time invariant dictionary field.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.DVTVd",
    "page": "API",
    "title": "FEMBase.DVTVd",
    "category": "Type",
    "text": "DVTVd(time => data::Dict)\n\nDiscrete, variable, time variant dictionary field.\n\n\n\n"
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
    "location": "api.html#FEMBase.assemble_mass_matrix!-Tuple{FEMBase.Problem,Array{FEMBase.Element{FEMBasis.Tet10},1},Any}",
    "page": "API",
    "title": "FEMBase.assemble_mass_matrix!",
    "category": "Method",
    "text": "assemble_mass_matrix!(problem, elements::Vector{Element{Tet10}}, time)\n\nAssemble Tet10 mass matrices using special method. If Tet10 has constant metric if can be integrated analytically to gain performance.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.field-Tuple{Any}",
    "page": "API",
    "title": "FEMBase.field",
    "category": "Method",
    "text": "field(x)\n\nCreate new discrete, constant, time invariant field from value x.\n\nExample\n\nf = field(1.0)\n\n\n\n"
},

{
    "location": "api.html#FEMBase.field-Tuple{Function}",
    "page": "API",
    "title": "FEMBase.field",
    "category": "Method",
    "text": "field(x::Function)\n\nCreate new, continuous, variable, time variant field from function x.\n\nExample\n\nf = field( (xi,t) -> xi[1]*t )\n\n\n\n"
},

{
    "location": "api.html#FEMBase.field-Union{Tuple{Dict{Int64,T}}, Tuple{T}} where T",
    "page": "API",
    "title": "FEMBase.field",
    "category": "Method",
    "text": "field(x::Dict)\n\nCreate new discrete, variable, time invariant dictionary field from dictionary x.\n\nExample\n\nX = Dict(1 => [0.0,0.0], 2 => [1.0,0.0])\nf = field(X)\n\n\n\n"
},

{
    "location": "api.html#FEMBase.field-Union{Tuple{N}, Tuple{Tuple{Vararg{T,N}}}, Tuple{T}} where T where N",
    "page": "API",
    "title": "FEMBase.field",
    "category": "Method",
    "text": "field(x::NTuple{N,T})\n\nCreate new discrete, variable, time invariant field from tuple x.\n\nExample\n\nf = field( (1.0, 2.0) )\n\n\n\n"
},

{
    "location": "api.html#FEMBase.field-Union{Tuple{N}, Tuple{T}, Tuple{Vararg{Pair{Float64,Tuple{Vararg{Pair{Int64,T},N}}},N} where N}} where T where N",
    "page": "API",
    "title": "FEMBase.field",
    "category": "Method",
    "text": "field(x::Pair{Float64, NTuple{N, Pair{Int64, T}}})\n\nCreate new discrete, variable, time variant dictionary field from x.\n\nExample\n\nX1 = (1 => 1.0, 2 => 2.0)\nX2 = (1 => 2.0, 2 => 3.0)\nf = field(0.0 => X1, 1.0 => X2)\n\n\n\n"
},

{
    "location": "api.html#FEMBase.field-Union{Tuple{N}, Tuple{T}, Tuple{Vararg{Pair{Float64,Tuple{Vararg{T,N}}},N} where N}} where T where N",
    "page": "API",
    "title": "FEMBase.field",
    "category": "Method",
    "text": "field(x::Pair{Float64,NTuple{N,T})\n\nCreate new discrete, variable, time variant field from pair x.\n\nExample\n\nf = field(1.0=>(1.0,2.0), 2.0=>(2.0,3.0))\n\n\n\n"
},

{
    "location": "api.html#FEMBase.field-Union{Tuple{T}, Tuple{Vararg{Pair{Float64,Dict{Int64,T}},N} where N}} where T",
    "page": "API",
    "title": "FEMBase.field",
    "category": "Method",
    "text": "field(t::Float64 => x::Dict, ...)\n\nCreate new discrete, variable, time variant dictionary field from pair of time and dictionary.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.field-Union{Tuple{T}, Tuple{Vararg{Pair{Float64,T},N} where N}} where T",
    "page": "API",
    "title": "FEMBase.field",
    "category": "Method",
    "text": "field(x::Pair{Float64,T})\n\nCreate new discrete, constant, time variant field from pair x.\n\nExample\n\nf = field(1.0=>1.0, 2.0=>2.0)\n\n\n\n"
},

{
    "location": "api.html#FEMBase.field-Union{Tuple{T}, Tuple{Vararg{Pair{Int64,T},N} where N}} where T",
    "page": "API",
    "title": "FEMBase.field",
    "category": "Method",
    "text": "field(x::Pair{Int64,T})\n\nCreate new discrete, variable, time invariant dictionary field from x.\n\nExample\n\nf = field(1 => 1.0, 2 => 2.0)\n\n\n\n"
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
    "location": "api.html#FEMBase.get_unknown_field_dimension-Tuple{FEMBase.Problem}",
    "page": "API",
    "title": "FEMBase.get_unknown_field_dimension",
    "category": "Method",
    "text": "get_unknown_field_dimension(problem)\n\nReturn the dimension of the unknown field of this problem.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_unknown_field_name-Union{Tuple{FEMBase.Problem{P}}, Tuple{P}} where P",
    "page": "API",
    "title": "FEMBase.get_unknown_field_name",
    "category": "Method",
    "text": "Return the name of the unknown field of this problem. \n\n\n\n"
},

{
    "location": "api.html#FEMBase.get_unknown_field_name-Union{Tuple{P}, Tuple{P}} where P<:FEMBase.AbstractProblem",
    "page": "API",
    "title": "FEMBase.get_unknown_field_name",
    "category": "Method",
    "text": "get_unknown_field_name(problem)\n\nDefault function if unknown field name is not defined for some problem.\n\n\n\n"
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
    "location": "api.html#FEMBase.update!-Tuple{FEMBase.DCTI,Any}",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "update!(f::DCTI, data)\n\nUpdate new value to field.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Tuple{FEMBase.DVTI,Any}",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "update!(f::DVTI, data)\n\nUpdate new value to field.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Tuple{FEMBase.Problem,FEMBase.Assembly,Array{T,1} where T,Array{T,1} where T}",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "update!(problem, assembly, u, la)\n\nUpdate the problem solution vector for assembly.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Union{Tuple{FEMBase.DCTV,Pair{Float64,T}}, Tuple{T}} where T",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "update!(f::DCTV, time => data)\n\nUpdate new value to field.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Union{Tuple{FEMBase.DVTId{T},Dict{Int64,T}}, Tuple{T}} where T",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "update!(field::DVTId, data::Dict)\n\nUpdate data to field.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Union{Tuple{FEMBase.DVTV,Pair{Float64,Tuple{Vararg{T,N}}}}, Tuple{N}, Tuple{T}} where T where N",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "update!(f::DVTV, time => data)\n\nUpdate new value to field.\n\n\n\n"
},

{
    "location": "api.html#FEMBase.update!-Union{Tuple{FEMBase.DVTVd,Pair{Float64,Dict{Int64,T}}}, Tuple{T}} where T",
    "page": "API",
    "title": "FEMBase.update!",
    "category": "Method",
    "text": "update!(f::DCTVd, time => data)\n\nUpdate new value to dictionary field.\n\n\n\n"
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
    "location": "api.html#FEMBasis.interpolate-Tuple{Any,Any}",
    "page": "API",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(a, b)\n\nA helper function for interpolate routines. Given iterables a and b, calculate c = aᵢbᵢ. Length of a can be less than b, but not vice versa.\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.interpolate-Tuple{FEMBase.DCTI,Vararg{Any,N} where N}",
    "page": "API",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(f::DCTI, args...)\n\nInterpolate constant, time-invariant DCTI field in time direction. That is trivially only the data itself.\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.interpolate-Tuple{FEMBase.DCTV,Any}",
    "page": "API",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(f::DCTV, time)\n\nInterpolate constant time variant DCTV field in time direction.\n\nNotes\n\nFirst check that is outside of range -> extrapolate Secondly check is \"exact match\" in time At last, find the correct bin and use linear interpolation\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.interpolate-Tuple{FEMBase.DVTId,Any}",
    "page": "API",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(f::DVTId, time)\n\nInterpolate DVTId, returns trivially the content as this is time invariant field.\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.interpolate-Union{Tuple{FEMBase.DVTI{N,T},Any,Any}, Tuple{N}, Tuple{T}} where T where N",
    "page": "API",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(f::DVTI, t, B)\n\nInterpolate variable, time-invariant DVTI field in time and spatial direction.\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.interpolate-Union{Tuple{FEMBase.DVTI{N,T},Any}, Tuple{N}, Tuple{T}} where T where N",
    "page": "API",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(f::DVTI, t)\n\nInterpolate variable, time-invariant DVTI field in time direction.\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.interpolate-Union{Tuple{FEMBase.DVTVd{T},Any}, Tuple{T}} where T",
    "page": "API",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(f::DVTVd, time)\n\nInterpolate variable, time variant DVTVd dictionary field in time direction.\n\nNotes\n\nFirst check that is outside of range -> extrapolate Secondly check is \"exact match\" in time At last, find the correct bin and use linear interpolation\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.interpolate-Union{Tuple{FEMBase.DVTV{N,T},Any,Any}, Tuple{N}, Tuple{T}} where T where N",
    "page": "API",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(f::DVTV, time, basis)\n\nInterpolate variable, time variant DVTV field in both time and spatial direction.\n\nNotes\n\nFirst check that is outside of range -> extrapolate Secondly check is \"exact match\" in time At last, find the correct bin and use linear interpolation\n\n\n\n"
},

{
    "location": "api.html#FEMBasis.interpolate-Union{Tuple{FEMBase.DVTV{N,T},Any}, Tuple{N}, Tuple{T}} where T where N",
    "page": "API",
    "title": "FEMBasis.interpolate",
    "category": "Method",
    "text": "interpolate(f::DVTV, time)\n\nInterpolate variable, time variant DVTV field in time direction.\n\nNotes\n\nFirst check that is outside of range -> extrapolate Secondly check is \"exact match\" in time At last, find the correct bin and use linear interpolation\n\n\n\n"
},

{
    "location": "api.html#Base.SparseArrays.sparse-Tuple{FEMBase.SparseMatrixCOO}",
    "page": "API",
    "title": "Base.SparseArrays.sparse",
    "category": "Method",
    "text": "Convert from COO format to CSC.\n\nParameters\n\ntol     used to drop near zero values less than tol.\n\n\n\n"
},

{
    "location": "api.html#Base.haskey-Tuple{FEMBase.Element,Any}",
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
    "location": "api.html#FEMBase.get_global_solution-Tuple{FEMBase.Problem,FEMBase.Assembly}",
    "page": "API",
    "title": "FEMBase.get_global_solution",
    "category": "Method",
    "text": "get_global_solution(problem, assembly)\n\nReturn a global solution (u, la) for a problem.\n\nNotes\n\nIf the length of solution vector != number of nodes, i.e. the field dimension is something else than 1, reshape vectors so that their length matches to the number of nodes. This helps to get nodal results easily.\n\n\n\n"
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
