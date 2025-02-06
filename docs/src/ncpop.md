# Noncommutative polynomial optimization

Noncommutative polynomial optimization concerns minimizing the smallest eigenvalue or the trace of a noncommutative polynomial subject to a tuple of noncommutative polynomial inequality constraints and equality constraints, which in general takes the form:

$$\mathrm{inf}_{\mathbf{x}\in\mathcal{B}(\mathcal{H})^n}\ \lambda_{\min}(f(\mathbf{x}))\ (\text{or } \mathrm{tr}(f(\mathbf{x}))) \ \text{ s.t. }\ g_1(\mathbf{x})\ge0,\ldots,g_m(\mathbf{x})\ge0,h_1(\mathbf{x})=0,\ldots,h_{\ell}(\mathbf{x})=0,$$

where $f,g_1,\ldots,g_m,h_1,\ldots,h_{\ell}\in\mathbb{R}\langle\mathbf{x}\rangle$ are noncommutative polynomials in noncommuting variables $\mathbf{x}$, and $\mathcal{H}$ is an (infinite dimensional) seperable Hilbert space.

To illustrate how to solve a noncommutative polynomial optimization problem with NCTSSOS, let us consider the following simple example.

```Julia
using DynamicPolynomials
using NCTSSOS
@ncpolyvar x[1:2]
f = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
ineq = [4 - x[1]^2 - x[2]^2]
eq = [x[1]*x[2] + x[2]*x[1] - 2]
pop = [f; ineq; eq]
d = 2 # set the relaxation order
opt,data = nctssos_first(pop, x, d, numeq=1, obj="eigen") # compute the first TS step of the NCTSSOS hierarchy
opt,data = nctssos_higher!(data) # compute higher TS steps of the NCTSSOS hierarchy
```

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
numeq | Specify the last **numeq** constraints to be equality constraints | 0
reducebasis | Reduce the monomial bases | false
obj | "eigen" (perform eigenvalue minimization) or "trace" (perform trace minimization) | "eigen"
TS | Types of chordal extensions used in term sparsity iterations: "block"(maximal chordal extension), "MD" (approximately smallest chordal extension), false (invalidating term sparsity exploitation) | "block"
normality | Impose normality condtions | false
soc | Impose state optimality condtions | false
merge | Merge overlapping PSD blocks | false
md | Parameter for tunning the merging strength | 3
solver | Specify an SDP solver: "Mosek" or "COSMO" | "Mosek"
cosmo_setting | Parameters for the COSMO solver: cosmo_para(eps_abs, eps_rel, max_iter) | cosmo_para(1e-5, 1e-5, 1e4)
QUIET | Silence the output| false
solve | Solve the SDP relaxation | true
Gram | Output Gram matrices | false
partition | Assume that the first **partition** variables commute with the remaining varibles | 0
constraint | nothing or "projection" (assume $x_i^2=x_i$ for all $i$) or "unipotent" (assume $x_i^2=1$ for all $i$) | nothing

## Correlative sparsity
The following is an example where one exploits correlative sparsity and term sparsity simultaneously.

```Julia
using DynamicPolynomials
using NCTSSOS
n = 10
@ncpolyvar x[1:n]
f = 0.0
for i = 1:n
    jset = max(1, i-5) : min(n, i+1)
    jset = setdiff(jset, i)
    f += (2x[i] + 5*x[i]^3 + 1)^2
    f -= sum([4x[i]*x[j] + 10x[i]^3*x[j] + 2x[j] + 4x[i]*x[j]^2 + 10x[i]^3*x[j]^2 + 2x[j]^2 for j in jset])
    f += sum([x[j]*x[k] + 2x[j]^2*x[k] + x[j]^2*x[k]^2 for j in jset for k in jset])
end
pop = [f]
for i = 1:n
    push!(pop, 1 - x[i]^2)
    push!(pop, x[i] - 1/3)
end
d = 3 # set the relaxation order
opt,data = cs_nctssos_first(pop, x, d) # compute the first TS step of the CS-NCTSSOS hierarchy
opt,data = cs_nctssos_higher!(data) # compute higher TS steps of the CS-NCTSSOS hierarchy
```

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
numeq | Specify the last **numeq** constraints to be equality constraints | 0
obj | "eigen" (perform eigenvalue minimization) or "trace" (perform trace minimization) | "eigen"
CS | Types of chordal extensions in exploiting correlative sparsity: "MF" (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation) | "MF"
TS | Types of chordal extensions used in term sparsity iterations: "block"(maximal chordal extension), "MD" (approximately smallest chordal extension), false (invalidating term sparsity exploitation) | "block"
solver | Specify an SDP solver: "Mosek" or "COSMO" | "Mosek"
cosmo_setting | Parameters for the COSMO solver: cosmo_para(eps_abs, eps_rel, max_iter) | cosmo_para(1e-5, 1e-5, 1e4)
QUIET | Silence the output| false
solve | Solve the SDP relaxation | true
Gram | Output Gram matrices | false
partition | Assume that the first **partition** variables commute with the remaining varibles | 0
constraint | nothing or "projection" (assume $x_i^2=x_i$ for all $i$) or "unipotent" (assume $x_i^2=1$ for all $i$) | nothing

## Methods
```@docs
nctssos_first
nctssos_higher!
cs_nctssos_first
cs_nctssos_higher!
```

### References
[1] [Exploiting Term Sparsity in Noncommutative Polynomial Optimization](https://arxiv.org/abs/2010.06956), 2021.  
[2] [Sparse polynomial optimization: theory and practice](https://arxiv.org/abs/2208.11158), 2023.  
[3] [Optimization of polynomials in non-commuting variables](https://link.springer.com/content/pdf/10.1007/978-3-319-33338-0.pdf), 2016. 