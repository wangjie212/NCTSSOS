# Sum-Of-Hermitian-Squares Optimization

A general Sum-of-Hermitian-squares optimization (including noncommutative polynomial optimization as a special case) problem takes the form:

$$\mathrm{inf}_{\mathbf{y}\in\mathbb{R}^n}\ \mathbf{c}^{\intercal}\mathbf{y}$$

$$\mathrm{s.t.}\ a_{k0}+y_1a_{k1}+\cdots+y_na_{kn}\in\mathrm{SOHS},\ k=1,\ldots,m.$$

where $\mathbf{c}\in\mathbb{R}^n$ and $a_{ki}\in\mathbb{R}\langle\mathbf{x}\rangle$ are noncommutative polynomials. In NCTSSOS, SOHS constraints could be handled with the routine **add_psatz!**:

```Julia
model,info = add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order; obj="eigen", CS=false, cliques=[], TS="block", SO=1, partition=0, constraint=nothing, QUIET=false, constrs=nothing)
```
where **nonneg** is a nonnegative noncommutative polynomial constrained to admit a Putinar's style SOHS representation on the noncommutative semialgebraic set defined by **ineq_cons** and **eq_cons**, and **SO** is the sparse order.

The following is a simple exmaple.

```Julia
using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using NCTSSOS
@ncpolyvar x[1:2]
f = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
g = 4 - x[1]^2 - x[2]^2
h = x[1]*x[2] + x[2]*x[1] - 2
d = 2 # set the relaxation order
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), false)
λ = @variable(model)
model,info1 = add_psatz!(model, f - λ, x, [g], [h], d, QUIET=true, TS=false, constrs="con1")
@objective(model, Max, λ)
optimize!(model)
objv = objective_value(model)
@show objv
```

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
obj | "eigen" (perform eigenvalue minimization) or "trace" (perform trace minimization) | "eigen"
CS | Types of chordal extensions in exploiting correlative sparsity: "MF" (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation) | "MF"
cliques | Use customized variable cliques | []
TS | Types of chordal extensions used in term sparsity iterations: "block"(maximal chordal extension), "signsymmetry" (sign symmetries), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations) | "block"
SO | Specify the sparse order | 1
QUIET | Silence the output| false
partition | Assume that the first **partition** variables commute with the remaining varibles | 0
constraint | nothing or "projection" (assume $x_i^2=x_i$ for all $i$) or "unipotent" (assume $x_i^2=1$ for all $i$) | nothing

## Methods
```@docs
add_psatz!
```