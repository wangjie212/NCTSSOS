var documenterSearchIndex = {"docs":
[{"location":"trace/#Trace-Polynomial-Optimization","page":"Trace Polynomial Optimization","title":"Trace Polynomial Optimization","text":"","category":"section"},{"location":"trace/","page":"Trace Polynomial Optimization","title":"Trace Polynomial Optimization","text":"Trace polynomial optimization concerns minimizing a trace polynomial subject to a tuple of noncommutative trace polynomial inequality constraints and equality constraints, which in general takes the form:","category":"page"},{"location":"trace/","page":"Trace Polynomial Optimization","title":"Trace Polynomial Optimization","text":"mathrminf_mathbfxinmathcalB(mathcalH)^n f(mathbfx) text st  g_1(mathbfx)ge0ldotsg_m(mathbfx)ge0h_1(mathbfx)=0ldotsh_ell(mathbfx)=0","category":"page"},{"location":"trace/","page":"Trace Polynomial Optimization","title":"Trace Polynomial Optimization","text":"where f is a trace polynomial in noncommuting variables mathbfx, g_1ldotsg_mh_1ldotsh_ell are noncommutative trace polynomials in noncommuting variables mathbfx, and mathcalH is an (infinite dimensional) seperable Hilbert space.","category":"page"},{"location":"trace/","page":"Trace Polynomial Optimization","title":"Trace Polynomial Optimization","text":"To illustrate how to solve a trace polynomial optimization problem with NCTSSOS, let us consider the following simple example.","category":"page"},{"location":"trace/","page":"Trace Polynomial Optimization","title":"Trace Polynomial Optimization","text":"using NCTSSOS\nn = 4 # number of variables\n# define the objective function: tr(x1y1) + tr(x1y2) + tr(x2y1) − tr(x2y2)\nsupp = [[[1;3]], [[1;4]], [[2;3]], [[2;4]]] # define the support data\ncoe = [-1; -1; -1; 1] # define the coefficient data\nd = 1 # set the relaxation order\nopt,data = ptraceopt_first(supp, coe, n, d, constraint=\"unipotent\") # compute the first TS step of the NCTSSOS hierarchy\nopt,data = ptraceopt_higher!(data) # compute higher TS steps of the NCTSSOS hierarchy","category":"page"},{"location":"trace/#Keyword-arguments","page":"Trace Polynomial Optimization","title":"Keyword arguments","text":"","category":"section"},{"location":"trace/","page":"Trace Polynomial Optimization","title":"Trace Polynomial Optimization","text":"Argument Description Default value\nnumeq Specify the last numeq constraints to be equality constraints 0\nTS Types of chordal extensions used in term sparsity iterations: \"block\"(maximal chordal extension), \"MD\" (approximately smallest chordal extension), false (invalidating term sparsity exploitation) \"block\"\nsolver Specify an SDP solver: \"Mosek\" or \"COSMO\" \"Mosek\"\ncosmo_setting Parameters for the COSMO solver: cosmopara(epsabs, epsrel, maxiter) cosmo_para(1e-5, 1e-5, 1e4)\nQUIET Silence the output false\nsolve Solve the SDP relaxation true\nGram Output Gram matrices false\nconstraint nothing or \"projection\" (assume x_i^2=x_i for all i) or \"unipotent\" (assume x_i^2=1 for all i) nothing","category":"page"},{"location":"trace/#References","page":"Trace Polynomial Optimization","title":"References","text":"","category":"section"},{"location":"trace/","page":"Trace Polynomial Optimization","title":"Trace Polynomial Optimization","text":"[1] Optimization over trace polynomials, 2021.  ","category":"page"},{"location":"state/#State-Polynomial-Optimization","page":"State Polynomial Optimization","title":"State Polynomial Optimization","text":"","category":"section"},{"location":"state/","page":"State Polynomial Optimization","title":"State Polynomial Optimization","text":"State polynomial optimization concerns minimizing a state polynomial subject to a tuple of noncommutative state polynomial inequality constraints and equality constraints, which in general takes the form:","category":"page"},{"location":"state/","page":"State Polynomial Optimization","title":"State Polynomial Optimization","text":"mathrminf_mathbfxinmathcalB(mathcalH)^n varsigmainmathcalS(mathcalH)^n f(mathbfx varsigma) text st  g_1(mathbfx varsigma)ge0ldotsg_m(mathbfx varsigma)ge0h_1(mathbfx varsigma)=0ldotsh_ell(mathbfx varsigma)=0","category":"page"},{"location":"state/","page":"State Polynomial Optimization","title":"State Polynomial Optimization","text":"where f is a state polynomial in noncommuting variables mathbfx, g_1ldotsg_mh_1ldotsh_ell are noncommutative state polynomials in noncommuting variables mathbfx, and mathcalH is an (infinite dimensional) seperable Hilbert space.","category":"page"},{"location":"state/","page":"State Polynomial Optimization","title":"State Polynomial Optimization","text":"To illustrate how to solve a state polynomial optimization problem with NCTSSOS, let us consider the following simple example.","category":"page"},{"location":"state/","page":"State Polynomial Optimization","title":"State Polynomial Optimization","text":"using NCTSSOS\nn = 4 # number of variables\n# define the objective function: ς(x1y1) + ς(x1y2) + ς(x2y1) − ς(x2y2)\nsupp = [[[1;3]], [[1;4]], [[2;3]], [[2;4]]] # define the coefficient data\ncoe = [1; 1; 1; -1] # define the coefficient data\nd = 1 # set the relaxation order\nopt,data = pstateopt_first(supp, coe, n, d, vargroup=[2;2], constraint=\"unipotent\") # vargroup=[2;2] constrains [xi, yj] = 0\nopt,data = pstateopt_higher!(data) # compute higher TS steps of the NCTSSOS hierarchy","category":"page"},{"location":"state/#Keyword-arguments","page":"State Polynomial Optimization","title":"Keyword arguments","text":"","category":"section"},{"location":"state/","page":"State Polynomial Optimization","title":"State Polynomial Optimization","text":"Argument Description Default value\nnumeq Specify the last numeq constraints to be equality constraints 0\nTS Types of chordal extensions used in term sparsity iterations: \"block\"(maximal chordal extension), \"MD\" (approximately smallest chordal extension), false (invalidating term sparsity exploitation) \"block\"\nvargroup Partition the variables into several groups such that variables from different groups commute [n]\nsolver Specify an SDP solver: \"Mosek\" or \"COSMO\" \"Mosek\"\ncosmo_setting Parameters for the COSMO solver: cosmopara(epsabs, epsrel, maxiter) cosmo_para(1e-5, 1e-5, 1e4)\nQUIET Silence the output false\nsolve Solve the SDP relaxation true\nGram Output Gram matrices false\nconstraint nothing or \"projection\" (assume x_i^2=x_i for all i) or \"unipotent\" (assume x_i^2=1 for all i) nothing","category":"page"},{"location":"state/#References","page":"State Polynomial Optimization","title":"References","text":"","category":"section"},{"location":"state/","page":"State Polynomial Optimization","title":"State Polynomial Optimization","text":"[1] State polynomials: positivity, optimization and nonlinear Bell inequalities, 2023. ","category":"page"},{"location":"ncpop/#Noncommutative-polynomial-optimization","page":"Noncommutative Polynomial Optimization","title":"Noncommutative polynomial optimization","text":"","category":"section"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"Noncommutative polynomial optimization concerns minimizing the smallest eigenvalue or the trace of a noncommutative polynomial subject to a tuple of noncommutative polynomial inequality constraints and equality constraints, which in general takes the form:","category":"page"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"mathrminf_mathbfxinmathcalB(mathcalH)^n lambda_min(f(mathbfx)) (textor  mathrmtr(f(mathbfx)))  text st  g_1(mathbfx)ge0ldotsg_m(mathbfx)ge0h_1(mathbfx)=0ldotsh_ell(mathbfx)=0","category":"page"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"where fg_1ldotsg_mh_1ldotsh_ellinmathbbRlanglemathbfxrangle are noncommutative polynomials in noncommuting variables mathbfx, and mathcalH is an (infinite dimensional) seperable Hilbert space.","category":"page"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"To illustrate how to solve a noncommutative polynomial optimization problem with NCTSSOS, let us consider the following simple example.","category":"page"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"using DynamicPolynomials\nusing NCTSSOS\n@ncpolyvar x[1:2]\nf = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2\nineq = [4 - x[1]^2 - x[2]^2]\neq = [x[1]*x[2] + x[2]*x[1] - 2]\npop = [f; ineq; eq]\nd = 2 # set the relaxation order\nopt,data = nctssos_first(pop, x, d, numeq=1, obj=\"eigen\") # compute the first TS step of the NCTSSOS hierarchy\nopt,data = nctssos_higher!(data) # compute higher TS steps of the NCTSSOS hierarchy","category":"page"},{"location":"ncpop/#Keyword-arguments","page":"Noncommutative Polynomial Optimization","title":"Keyword arguments","text":"","category":"section"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"Argument Description Default value\nnumeq Specify the last numeq constraints to be equality constraints 0\nreducebasis Reduce the monomial bases false\nobj \"eigen\" (perform eigenvalue minimization) or \"trace\" (perform trace minimization) \"eigen\"\nTS Types of chordal extensions used in term sparsity iterations: \"block\"(maximal chordal extension), \"MD\" (approximately smallest chordal extension), false (invalidating term sparsity exploitation) \"block\"\nnormality Impose normality condtions false\nsoc Impose state optimality condtions false\nmerge Merge overlapping PSD blocks false\nmd Parameter for tunning the merging strength 3\nsolver Specify an SDP solver: \"Mosek\" or \"COSMO\" \"Mosek\"\ncosmo_setting Parameters for the COSMO solver: cosmopara(epsabs, epsrel, maxiter) cosmo_para(1e-5, 1e-5, 1e4)\nQUIET Silence the output false\nsolve Solve the SDP relaxation true\nGram Output Gram matrices false\npartition Assume that the first partition variables commute with the remaining varibles 0\nconstraint nothing or \"projection\" (assume x_i^2=x_i for all i) or \"unipotent\" (assume x_i^2=1 for all i) nothing","category":"page"},{"location":"ncpop/#Correlative-sparsity","page":"Noncommutative Polynomial Optimization","title":"Correlative sparsity","text":"","category":"section"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"The following is an example where one exploits correlative sparsity and term sparsity simultaneously.","category":"page"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"using DynamicPolynomials\nusing NCTSSOS\nn = 10\n@ncpolyvar x[1:n]\nf = 0.0\nfor i = 1:n\n    jset = max(1, i-5) : min(n, i+1)\n    jset = setdiff(jset, i)\n    f += (2x[i] + 5*x[i]^3 + 1)^2\n    f -= sum([4x[i]*x[j] + 10x[i]^3*x[j] + 2x[j] + 4x[i]*x[j]^2 + 10x[i]^3*x[j]^2 + 2x[j]^2 for j in jset])\n    f += sum([x[j]*x[k] + 2x[j]^2*x[k] + x[j]^2*x[k]^2 for j in jset for k in jset])\nend\npop = [f]\nfor i = 1:n\n    push!(pop, 1 - x[i]^2)\n    push!(pop, x[i] - 1/3)\nend\nd = 3 # set the relaxation order\nopt,data = cs_nctssos_first(pop, x, d) # compute the first TS step of the CS-NCTSSOS hierarchy\nopt,data = cs_nctssos_higher!(data) # compute higher TS steps of the CS-NCTSSOS hierarchy","category":"page"},{"location":"ncpop/#Keyword-arguments-2","page":"Noncommutative Polynomial Optimization","title":"Keyword arguments","text":"","category":"section"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"Argument Description Default value\nnumeq Specify the last numeq constraints to be equality constraints 0\nobj \"eigen\" (perform eigenvalue minimization) or \"trace\" (perform trace minimization) \"eigen\"\nCS Types of chordal extensions in exploiting correlative sparsity: \"MF\" (approximately smallest chordal extension), \"NC\" (not performing chordal extension), false (invalidating correlative sparsity exploitation) \"MF\"\nTS Types of chordal extensions used in term sparsity iterations: \"block\"(maximal chordal extension), \"MD\" (approximately smallest chordal extension), false (invalidating term sparsity exploitation) \"block\"\nsolver Specify an SDP solver: \"Mosek\" or \"COSMO\" \"Mosek\"\ncosmo_setting Parameters for the COSMO solver: cosmopara(epsabs, epsrel, maxiter) cosmo_para(1e-5, 1e-5, 1e4)\nQUIET Silence the output false\nsolve Solve the SDP relaxation true\nGram Output Gram matrices false\npartition Assume that the first partition variables commute with the remaining varibles 0\nconstraint nothing or \"projection\" (assume x_i^2=x_i for all i) or \"unipotent\" (assume x_i^2=1 for all i) nothing","category":"page"},{"location":"ncpop/#Methods","page":"Noncommutative Polynomial Optimization","title":"Methods","text":"","category":"section"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"nctssos_first\nnctssos_higher!\ncs_nctssos_first\ncs_nctssos_higher!","category":"page"},{"location":"ncpop/#NCTSSOS.nctssos_first","page":"Noncommutative Polynomial Optimization","title":"NCTSSOS.nctssos_first","text":"opt,data = nctssos_first(f::Polynomial{false, T} where T<:Number, x::Vector{PolyVar{false}};\n    newton=true, reducebasis=true, TS=\"block\", obj=\"eigen\", merge=false, md=3, solve=true, Gram=false, QUIET=false)\n\nCompute the first step of the NCTSSOS hierarchy for unconstrained noncommutative polynomial optimization. If newton=true, then compute a monomial basis by the Newton chip method. If reducebasis=true, then remove monomials from the monomial basis by diagonal inconsistency. If TS=\"block\", use maximal chordal extensions; if TS=\"MD\", use approximately smallest chordal extensions. If obj=\"eigen\", minimize the eigenvalue; if obj=\"trace\", then minimize the trace. If merge=true, perform the PSD block merging. Return the optimum and other auxiliary data.\n\nArguments\n\nf: the objective function for unconstrained noncommutative polynomial optimization.\nx: the set of noncommuting variables.\nmd: the tunable parameter for merging blocks.\n\n\n\n\n\n","category":"function"},{"location":"ncpop/#NCTSSOS.nctssos_higher!","page":"Noncommutative Polynomial Optimization","title":"NCTSSOS.nctssos_higher!","text":"opt,data = nctssos_higher!(data, TS=\"block\", merge=false, md=3, solve=true, Gram=false, QUIET=false)\n\nCompute higher steps of the NCTSSOS hierarchy. Return the optimum and other auxiliary data.\n\n\n\n\n\n","category":"function"},{"location":"ncpop/#NCTSSOS.cs_nctssos_first","page":"Noncommutative Polynomial Optimization","title":"NCTSSOS.cs_nctssos_first","text":"opt,data = cs_nctssos_first(pop, x, d; numeq=0, CS=\"MF\", TS=\"block\", QUIET=false, obj=\"eigen\", solve=true, Gram=false)\n\nCompute the first step of the CS-NCTSSOS hierarchy for constrained noncommutative polynomial optimization with relaxation order d. Return the optimum and other auxiliary data.\n\nArguments\n\npop: the vector of the objective function, inequality constraints, and equality constraints.\nx: the set of noncommuting variables.\nd: the relaxation order of the moment-SOHS hierarchy.\nnumeq: the number of equality constraints.\n\n\n\n\n\nopt,data = cs_nctssos_first(supp::Vector{Vector{Vector{UInt16}}}, coe, n::Int, d::Int; numeq=0,\nCS=\"MF\", TS=\"block\", QUIET=false, obj=\"eigen\", solve=true, Gram=false)\n\nCompute the first step of the CS-NCTSSOS hierarchy for constrained noncommutative polynomial optimization with relaxation order d. Here the polynomial optimization problem is defined by supp and coe, corresponding to the supports and coeffients of pop respectively. Return the optimum and other auxiliary data.\n\nArguments\n\nsupp: the supports of the polynomial optimization problem.\ncoe: the coeffients of the polynomial optimization problem.\nd: the relaxation order of the moment-SOHS hierarchy.\nnumeq: the number of equality constraints.\n\n\n\n\n\n","category":"function"},{"location":"ncpop/#NCTSSOS.cs_nctssos_higher!","page":"Noncommutative Polynomial Optimization","title":"NCTSSOS.cs_nctssos_higher!","text":"opt,data = cs_nctssos_higher!(data; TS=\"block\", QUIET=false, solve=true, Gram=false)\n\nCompute higher steps of the CS-NCTSSOS hierarchy. Return the optimum and other auxiliary data.\n\n\n\n\n\n","category":"function"},{"location":"ncpop/#References","page":"Noncommutative Polynomial Optimization","title":"References","text":"","category":"section"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"[1] Exploiting Term Sparsity in Noncommutative Polynomial Optimization, 2021.     [2] Sparse polynomial optimization: theory and practice, 2023.    [3] Optimization of polynomials in non-commuting variables, 2016. ","category":"page"},{"location":"sohs/#Sum-Of-Hermitian-Squares-Optimization","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"","category":"section"},{"location":"sohs/","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"A general Sum-of-Hermitian-squares optimization (including noncommutative polynomial optimization as a special case) problem takes the form:","category":"page"},{"location":"sohs/","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"mathrminf_mathbfyinmathbbR^n mathbfc^intercalmathbfy","category":"page"},{"location":"sohs/","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"mathrmst a_k0+y_1a_k1+cdots+y_na_kninmathrmSOHS k=1ldotsm","category":"page"},{"location":"sohs/","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"where mathbfcinmathbbR^n and a_kiinmathbbRlanglemathbfxrangle are noncommutative polynomials. In NCTSSOS, SOHS constraints could be handled with the routine add_psatz!:","category":"page"},{"location":"sohs/","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"model,info = add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order; obj=\"eigen\", CS=false, cliques=[], TS=\"block\", SO=1, partition=0, constraint=nothing, QUIET=false, constrs=nothing)","category":"page"},{"location":"sohs/","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"where nonneg is a nonnegative noncommutative polynomial constrained to admit a Putinar's style SOHS representation on the noncommutative semialgebraic set defined by ineq_cons and eq_cons, and SO is the sparse order.","category":"page"},{"location":"sohs/","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"The following is a simple exmaple.","category":"page"},{"location":"sohs/","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"using JuMP\nusing MosekTools\nusing DynamicPolynomials\nusing MultivariatePolynomials\nusing NCTSSOS\n@ncpolyvar x[1:2]\nf = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2\ng = 4 - x[1]^2 - x[2]^2\nh = x[1]*x[2] + x[2]*x[1] - 2\nd = 2 # set the relaxation order\nmodel = Model(optimizer_with_attributes(Mosek.Optimizer))\nset_optimizer_attribute(model, MOI.Silent(), false)\nλ = @variable(model)\nmodel,info1 = add_psatz!(model, f - λ, x, [g], [h], d, QUIET=true, TS=false, constrs=\"con1\")\n@objective(model, Max, λ)\noptimize!(model)\nobjv = objective_value(model)\n@show objv","category":"page"},{"location":"sohs/#Keyword-arguments","page":"Sum-Of-Hermitian-Squares Optimization","title":"Keyword arguments","text":"","category":"section"},{"location":"sohs/","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"Argument Description Default value\nobj \"eigen\" (perform eigenvalue minimization) or \"trace\" (perform trace minimization) \"eigen\"\nCS Types of chordal extensions in exploiting correlative sparsity: \"MF\" (approximately smallest chordal extension), \"NC\" (not performing chordal extension), false (invalidating correlative sparsity exploitation) \"MF\"\ncliques Use customized variable cliques []\nTS Types of chordal extensions used in term sparsity iterations: \"block\"(maximal chordal extension), \"signsymmetry\" (sign symmetries), \"MD\" (approximately smallest chordal extension), false (invalidating term sparsity iterations) \"block\"\nSO Specify the sparse order 1\nQUIET Silence the output false\npartition Assume that the first partition variables commute with the remaining varibles 0\nconstraint nothing or \"projection\" (assume x_i^2=x_i for all i) or \"unipotent\" (assume x_i^2=1 for all i) nothing","category":"page"},{"location":"sohs/#Methods","page":"Sum-Of-Hermitian-Squares Optimization","title":"Methods","text":"","category":"section"},{"location":"sohs/","page":"Sum-Of-Hermitian-Squares Optimization","title":"Sum-Of-Hermitian-Squares Optimization","text":"add_psatz!","category":"page"},{"location":"#NCTSSOS","page":"Home","title":"NCTSSOS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"NCTSSOS aims to provide a user-friendly and efficient tool for solving optimization problems with non-commutative/trace/state polynomials, which is based on the structured moment-SOHS hierarchy.","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Jie Wang, Academy of Mathematics and Systems Science, Chinese Academy of Sciences.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"NCTSSOS could be installed by running","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/wangjie212/NCTSSOS","category":"page"},{"location":"#Related-packages","page":"Home","title":"Related packages","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DynamicPolynomials: Polynomial definition\nMultivariatePolynomials: Polynomials manipulations\nTSSOS: Commutative polynomial optimization\nChordalGraph: Chordal graphs and chordal extentions","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"[1] Exploiting Term Sparsity in Noncommutative Polynomial Optimization, 2021.     [2] Sparse polynomial optimization: theory and practice, 2023.      [3] Optimization of polynomials in non-commuting variables, 2016.   [4] State polynomials: positivity, optimization and nonlinear Bell inequalities, 2023. ","category":"page"},{"location":"bell/#Bell-inequalities","page":"Bell inequalities","title":"Bell inequalities","text":"","category":"section"},{"location":"bell/#Linear-bell-inequalities","page":"Bell inequalities","title":"Linear bell inequalities","text":"","category":"section"},{"location":"bell/","page":"Bell inequalities","title":"Bell inequalities","text":"The CHSH inequality:","category":"page"},{"location":"bell/","page":"Bell inequalities","title":"Bell inequalities","text":"using DynamicPolynomials\nusing NCTSSOS\n@ncpolyvar x[1:2]\n@ncpolyvar y[1:2]\nf = x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2]\nopt,data = nctssos_first([-f], [x;y], 1, TS=false, partition=2, constraint=\"unipotent\")","category":"page"},{"location":"bell/","page":"Bell inequalities","title":"Bell inequalities","text":"The I_3322 inequality:","category":"page"},{"location":"bell/","page":"Bell inequalities","title":"Bell inequalities","text":"using DynamicPolynomials\nusing NCTSSOS\n@ncpolyvar x[1:3]\n@ncpolyvar y[1:3]\nf = x[1]*(y[1] + y[2] + y[3]) + x[2]*(y[1] + y[2] - y[3]) + x[3]*(y[1] - y[2]) - x[1] - 2*y[1] - y[2]\nopt,data = nctssos_first([-f], [x;y], 3, TS=false, normality=true, partition=3, constraint=\"projection\")","category":"page"},{"location":"bell/#Nonlinear-bell-inequalities","page":"Bell inequalities","title":"Nonlinear bell inequalities","text":"","category":"section"}]
}
