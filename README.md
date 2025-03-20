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
f = 1.0 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2] + x[2]*x[1] + x[2]*x[3] + x[3]*x[2]

pop =  PolyOpt(f)

solver_config_dense = SolverConfig(optimizer=Clarabel.Optimizer)

result_dense = cs_nctssos(pop, solver_config_dense)
```

To use correlative sparsity

```Julia
result_cs = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF()))
```

To use term sparsity

```Julia
result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))
```


### Constrained non-commutative polynomial optimization
Taking the objective $f=2-x_1^2+x_1x_2^2x_1-x_2^2$ and constraints $g=4-x_1^2-x_2^2\ge0$, $h=x_1x_2+x_2x_1-2=0$ as an example, to solve the first step of the NCTSSOS hierarchy, run

```Julia
@ncpolyvar x[1:2]
f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
g = 4.0 - x[1]^2 - x[2]^2
h1 = x[1]*x[2] + x[2]*x[1] - 2.0
h2 = -h1

pop =  PolyOpt(f, [g, h1, h2])

result_dense = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer))
```

To use Correlative Sparsity

```Julia
result_cs = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF()))
```

To exploit correlative sparsity and term sparsity simultaneously, do

```Julia
result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))
```

To exploit higher iteration of term sparsity, do

```Julia
result_cs_ts_higher = cs_nctssos_higher(pop, result_cs_ts, SolverConfig(optimizer=Clarabel.Optimizer, mom_order=2, cs_algo=MF(), ts_algo=MMD()))
```

## References
[1] [Exploiting Term Sparsity in Noncommutative Polynomial Optimization](https://arxiv.org/abs/2010.06956), 2021.
[2] [Sparse polynomial optimization: theory and practice](https://arxiv.org/abs/2208.11158), 2023.
[3] [State polynomials: positivity, optimization and nonlinear Bell inequalities](https://arxiv.org/abs/2301.12513), 2023.

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn
