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
opt,data = nctssos_first(f, x, obj="eigen")
```

To sovle higher steps of the NCTSSOS hierarchy, repeatedly run

```Julia
opt,data = nctssos_higher!(data)
```

Options:  
**obj**: "eigen" by default (perform eigenvalue minimization), "trace" (perform trace minimization)  
**TS**: "block" by default (maximal chordal extension), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations) 
**partition**: specify that the first *partition* variables commute with the remaining variables  
**comm_var**: specify that the first *comm_var* variables commute each other  
**constraint**: nothing by default or "projection" (satisfying $x_i^2=x_i$) or "unipotent" (satisfying $x_i^2=1$)  

### Constrained non-commutative polynomial optimization
Taking the objective $f=2-x_1^2+x_1x_2^2x_1-x_2^2$ and constraints $g=4-x_1^2-x_2^2\ge0$, $h=x_1x_2+x_2x_1-2=0$ as an example, to solve the first step of the NCTSSOS hierarchy, run

```Julia
@ncpolyvar x[1:2]
f = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
ineq = [4 - x[1]^2 - x[2]^2]
eq = [x[1]*x[2] + x[2]*x[1] - 2]
pop = [f; ineq; eq]
d = 2 # set the relaxation order
opt,data = nctssos_first(pop, x, d, numeq=1, obj="eigen")
```

To solve higher steps of the NCTSSOS hierarchy, repeatedly run

```Julia
opt,data = nctssos_higher!(data)
```

Options:  
**obj**: "eigen" by default (perform eigenvalue minimization), "trace" (perform trace minimization)  
**TS**: "block" by default (maximal chordal extension), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations)  
**partition**: specify that the first *partition* variables commute with the remaining variables  
**comm_var**: specify that the first *comm_var* variables commute each other  
**constraint**: nothing by default or "projection" (satisfying $x_i^2=x_i$) or "unipotent" (satisfying $x_i^2=1$)  

To exploit correlative sparsity and term sparsity simultaneously, run

```Julia
opt,data = cs_nctssos_first(pop, x, d, obj="eigen")
opt,data = cs_nctssos_higher!(data)
```

### Trace polynomial optimization
Check out /examples/traceopt.jl.

### State polynomial optimization
Check out /examples/stateopt.jl for state polynomial optimization over real numbers and /examples/state_complex.jl for state polynomial optimization over complex numbers.

## References
[1] [Exploiting Term Sparsity in Noncommutative Polynomial Optimization](https://arxiv.org/abs/2010.06956), 2021.  
[2] [Sparse polynomial optimization: theory and practice](https://arxiv.org/abs/2208.11158), 2023.  
[3] [State polynomials: positivity, optimization and nonlinear Bell inequalities](https://arxiv.org/abs/2301.12513), 2023. 

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn
