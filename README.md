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
- [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl)
- [CliqueTrees](https://github.com/AlgebraicJulia/CliqueTrees.jl)
- [ChordalGraph](https://github.com/wangjie212/ChordalGraph)

NCTSSOS has been tested on Ubuntu, MacOS, and Windows.

## Usage
### Unconstrained non-commutative polynomial optimization
Taking $f=1+x_1^4+x_2^4+x_3^4+x_1x_2+x_2x_1+x_2x_3+x_3x_2$ as an example, to compute the first step of the NCTSSOS hierarchy, run

```Julia
using NCTSSOS
using DynamicPolynomials
using Clarabel

@ncpolyvar x[1:3]
f = 1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2] + x[2]*x[1] + x[2]*x[3] + x[3]*x[2]

pop = PolynomialOptimizationProblem(f)

mom_order = 2

problem = cs_nctssos(pop, mom_order)
myans = solve_problem(problem,Clarabel.Optimizer)
```

To use term sparsity

```Julia
using CliqueTrees
ts_order = 1
ts_algo = MMD()

problem = cs_nctssos(pop, mom_order, ts_order=ts_order, ts_algo=ts_algo)
```

To use correlative sparsity

```Julia
using CorrelativeSparsity
cs_algo = MMD()

problem = cs_nctssos(pop, mom_order, cs_algo=cs_algo)
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

problem = cs_nctssos(pop, mom_order)
myans = solve_problem(problem,Clarabel.Optimizer)
```

To use Correlative Sparsity

```Julia
cs_algo = MD()

problem = cs_nctssos(pop, order, cs_algo=cs_algo)
myans = solve_problem(problem,Clarabel.Optimizer)
```

To exploit correlative sparsity and term sparsity simultaneously, do

```Julia
ts_algo = MMD()
ts_algo = 1

problem = cs_nctssos(pop, order, cs_algo=cs_algo,ts_order=1, ts_algo=ts_algo)
myans = solve_problem(problem,Clarabel.Optimizer)
```

## References
[1] [Exploiting Term Sparsity in Noncommutative Polynomial Optimization](https://arxiv.org/abs/2010.06956), 2021.
[2] [Sparse polynomial optimization: theory and practice](https://arxiv.org/abs/2208.11158), 2023.
[3] [State polynomials: positivity, optimization and nonlinear Bell inequalities](https://arxiv.org/abs/2301.12513), 2023.

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn
