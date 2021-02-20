var documenterSearchIndex = {"docs":
[{"location":"ncpop/#Noncommutative-polynomial-optimization","page":"Noncommutative Polynomial Optimization","title":"Noncommutative polynomial optimization","text":"","category":"section"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"nctssos_first\r\nnctssos_higher!\r\ncs_nctssos_first\r\ncs_nctssos_higher!","category":"page"},{"location":"ncpop/#NCTSSOS.nctssos_first","page":"Noncommutative Polynomial Optimization","title":"NCTSSOS.nctssos_first","text":"opt,data = nctssos_first(f::Polynomial{false, T} where T<:Number, x::Vector{PolyVar{false}};\n    newton=true, reducebasis=true, TS=\"block\", obj=\"eigen\", merge=false, md=3, QUIET=false)\n\nCompute the first step of the NCTSSOS hierarchy for unconstrained noncommutative polynomial optimization. If newton=true, then compute a monomial basis by the Newton chip method. If reducebasis=true, then remove monomials from the monomial basis by diagonal inconsistency. If TS=\"block\", use maximal chordal extensions; if TS=\"MD\", use approximately smallest chordal extensions. If obj=\"eigen\", minimize the eigenvalue; if obj=\"trace\", then minimize the trace. If merge=true, perform the PSD block merging. Return the optimum and other auxiliary data.\n\nArguments\n\nf: the objective function for unconstrained noncommutative polynomial optimization.\nx: the set of noncommuting variables.\nmd: the tunable parameter for merging blocks.\n\n\n\n\n\n","category":"function"},{"location":"ncpop/#NCTSSOS.nctssos_higher!","page":"Noncommutative Polynomial Optimization","title":"NCTSSOS.nctssos_higher!","text":"opt,data = nctssos_higher!(data, TS=\"block\", merge=false, md=3, QUIET=false)\n\nCompute higher steps of the NCTSSOS hierarchy. Return the optimum and other auxiliary data.\n\n\n\n\n\n","category":"function"},{"location":"ncpop/#NCTSSOS.cs_nctssos_first","page":"Noncommutative Polynomial Optimization","title":"NCTSSOS.cs_nctssos_first","text":"opt,data = cs_nctssos_first(pop, x, d; numeq=0, CS=\"MF\", TS=\"block\", merge=false, md=3,\nQUIET=false, obj=\"eigen\", solve=true)\n\nCompute the first step of the CS-NCTSSOS hierarchy for constrained noncommutative polynomial optimization with relaxation order d. Return the optimum and other auxiliary data.\n\nArguments\n\npop: the vector of the objective function, inequality constraints, and equality constraints.\nx: the set of noncommuting variables.\nd: the relaxation order of the moment-SOHS hierarchy.\nnumeq: the number of equality constraints.\n\n\n\n\n\nopt,data = cs_nctssos_first(supp::Vector{Vector{Vector{UInt16}}}, coe, n::Int, d::Int; numeq=0,\nCS=\"MF\", TS=\"block\", merge=false, md=3, QUIET=false, obj=\"eigen\", solve=true)\n\nCompute the first step of the CS-NCTSSOS hierarchy for constrained noncommutative polynomial optimization with relaxation order d. Here the polynomial optimization problem is defined by supp and coe, corresponding to the supports and coeffients of pop respectively. Return the optimum and other auxiliary data.\n\nArguments\n\nsupp: the supports of the polynomial optimization problem.\ncoe: the coeffients of the polynomial optimization problem.\nd: the relaxation order of the moment-SOHS hierarchy.\nnumeq: the number of equality constraints.\n\n\n\n\n\n","category":"function"},{"location":"ncpop/#NCTSSOS.cs_nctssos_higher!","page":"Noncommutative Polynomial Optimization","title":"NCTSSOS.cs_nctssos_higher!","text":"opt,data = cs_nctssos_higher!(data; TS=\"block\", QUIET=false, merge=false, md=3, solve=true)\n\nCompute higher steps of the CS-NCTSSOS hierarchy. Return the optimum and other auxiliary data.\n\n\n\n\n\n","category":"function"},{"location":"#NCTSSOS","page":"Home","title":"NCTSSOS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"NCTSSOS is a noncommutative polynomial optimization package based on the sparsity adapted moment-SOHS hierarchies.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Jie Wang,  LAAS-CNRS.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"NCTSSOS is simply installed by running","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/wangjie212/NCTSSOS","category":"page"},{"location":"#Related-packages","page":"Home","title":"Related packages","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DynamicPolynomials: Polynomial definition\nMultivariatePolynomials: Polynomials manipulations\nTSSOS: Commutative polynomial optimization\nChordalGraph: Chordal graphs and chordal extentions\nSparseJSR: Computing joint spetral radius","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Exploiting term sparsity in Noncommutative Polynomial Optimization, Jie Wang and Victor Magron, 2020.\nTSSOS: A Moment-SOS hierarchy that exploits term sparsity, Jie Wang, Victor Magron and Jean B. Lasserre, 2020.\nChordal-TSSOS: a moment-SOS hierarchy that exploits term sparsity with chordal extension, Jie Wang, Victor Magron and Jean B. Lasserre, 2020.\nCS-TSSOS: Correlative and term sparsity for large-scale polynomial optimization, Jie Wang, Victor Magron, Jean B. Lasserre and Ngoc H. A. Mai, 2020.","category":"page"},{"location":"spop/#Noncommutative-polynomial-optimization","page":"Noncommutative polynomial optimization","title":"Noncommutative polynomial optimization","text":"","category":"section"},{"location":"spop/","page":"Noncommutative polynomial optimization","title":"Noncommutative polynomial optimization","text":"nctssos_first\r\nnctssos_higher!\r\ncs_nctssos_first\r\ncs_nctssos_higher!","category":"page"}]
}
