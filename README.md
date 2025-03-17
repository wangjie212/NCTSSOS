# NCTSSOS
NCTSSOS aims to provide a user-friendly and efficient tool for solving optimization problems with non-commutative/trace/state polynomials, which is based on the structured moment-SOHS hierarchy. To use NCTSSOS in Julia, run
```Julia
pkg> add https://github.com/wangjie212/NCTSSOS
 ```

 | **Documentation** |
 |:-----------------:|
 | [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://wangjie212.github.io/NCTSSOS/dev) |

## Dependencies
- [Julia](https://julialang.org/)
- [JuMP](https://github.com/jump-dev/JuMP.jl)
- [Mosek](https://www.mosek.com/) or [COSMO](https://github.com/oxfordcontrol/COSMO.jl)
- [ChordalGraph](https://github.com/wangjie212/ChordalGraph)

NCTSSOS has been tested on Ubuntu and Windows.

## Usage
### Unconstrained non-commutative polynomial optimization
Taking $f=1+x_1^4+x_2^4+x_3^4+x_1x_2+x_2x_1+x_2x_3+x_3x_2$ as an example, to compute the first step of the NCTSSOS hierarchy, run

```Julia
using NCTSSOS
using DynamicPolynomials
@ncpolyvar x[1:3]
f = 1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2] + x[2]*x[1] + x[2]*x[3] + x[3]*x[2]

pop = PolynomialOptimizationProblem(f, typeof(f)[])

corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, nothing)

cliques_term_sparsities = [
    [TermSparsity(Vector{Monomial{false}}(),[basis]) for basis in idx_basis]
        for idx_basis in corr_sparsity.cliques_idcs_bases
]

moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, cliques_term_sparsities)

set_optimizer(moment_problem.model, Clarabel.Optimizer)
optimize!(moment_problem.model)
```

To use term sparsity 

```Julia
using CliqueTrees
ts_algo = MMD()

cliques_objective = [reduce(+, [issubset(effective_variables(mono), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

# prepare the support for each term sparse localizing moment
initial_activated_supp = [
    # why does order matter here?
    sorted_union(symmetric_canonicalize.(monomials(obj_part)), [neat_dot(b, b) for b in idcs_bases[1]])
    for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)
]

cliques_term_sparsities_s = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
    [iterate_term_sparse_supp(activated_supp, poly, basis, ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
end

moment_problem_s = moment_relax(pop, corr_sparsity.cliques_cons, cliques_term_sparsities_s)
set_optimizer(moment_problem_s.model, Clarabel.Optimizer)
optimize!(moment_problem_s.model)
```

### Constrained non-commutative polynomial optimization
Taking the objective $f=2-x_1^2+x_1x_2^2x_1-x_2^2$ and constraints $g=4-x_1^2-x_2^2\ge0$, $h=x_1x_2+x_2x_1-2=0$ as an example, to solve the first step of the NCTSSOS hierarchy, run

```Julia
@ncpolyvar x[1:2]
f = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
g = 4 - x[1]^2 - x[2]^2
h = x[1]*x[2] + x[2]*x[1] - 2
h2 = -h
order = 2

pop = PolynomialOptimizationProblem(f, [g, h1, h2])
```

To use Correlative Sparsity

```Julia
cs_algo = MD()
corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, cs_algo)

cliques_term_sparsities = [
    [TermSparsity(Vector{Monomial{false}}(),[basis]) for basis in idx_basis]
    for idx_basis in corr_sparsity.cliques_idcs_bases
]

moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, cliques_term_sparsities)

set_optimizer(moment_problem.model, Clarabel.Optimizer)
optimize!(moment_problem.model)
```

To exploit correlative sparsity and term sparsity simultaneously, do

```Julia
ts_algo = MMD()

cliques_objective = [reduce(+, [issubset(effective_variables(mono), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

# prepare the support for each term sparse localizing moment
initial_activated_supp = [
    # why does order matter here?
    sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, pop.constraints[cons_idx]), [neat_dot(b, b) for b in idcs_bases[1]])
    for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)
]

cliques_term_sparsities_s = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
    [iterate_term_sparse_supp(activated_supp, poly, basis, ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
end

moment_problem_s = moment_relax(pop, corr_sparsity.cliques_cons, cliques_term_sparsities_s)
set_optimizer(moment_problem_s.model, Clarabel.Optimizer)
optimize!(moment_problem_s.model)
```

## References
[1] [Exploiting Term Sparsity in Noncommutative Polynomial Optimization](https://arxiv.org/abs/2010.06956), 2021.  
[2] [Sparse polynomial optimization: theory and practice](https://arxiv.org/abs/2208.11158), 2023.  
[3] [State polynomials: positivity, optimization and nonlinear Bell inequalities](https://arxiv.org/abs/2301.12513), 2023. 

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn
