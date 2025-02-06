# Trace Polynomial Optimization

Trace polynomial optimization concerns minimizing a trace polynomial subject to a tuple of noncommutative trace polynomial inequality constraints and equality constraints, which in general takes the form:

$$\mathrm{inf}_{\mathbf{x}\in\mathcal{B}(\mathcal{H})^n}\ f(\mathbf{x})\ \text{ s.t. }\ g_1(\mathbf{x})\ge0,\ldots,g_m(\mathbf{x})\ge0,h_1(\mathbf{x})=0,\ldots,h_{\ell}(\mathbf{x})=0,$$

where $f$ is a trace polynomial in noncommuting variables $\mathbf{x}$, $g_1,\ldots,g_m,h_1,\ldots,h_{\ell}$ are noncommutative trace polynomials in noncommuting variables $\mathbf{x}$, and $\mathcal{H}$ is an (infinite dimensional) seperable Hilbert space.

To illustrate how to solve a trace polynomial optimization problem with NCTSSOS, let us consider the following simple example.

```Julia
using NCTSSOS
n = 4 # number of variables
# define the objective function: tr(x1y1) + tr(x1y2) + tr(x2y1) − tr(x2y2)
supp = [[[1;3]], [[1;4]], [[2;3]], [[2;4]]] # define the support data
coe = [-1; -1; -1; 1] # define the coefficient data
d = 1 # set the relaxation order
opt,data = ptraceopt_first(supp, coe, n, d, constraint="unipotent") # compute the first TS step of the NCTSSOS hierarchy
opt,data = ptraceopt_higher!(data) # compute higher TS steps of the NCTSSOS hierarchy
```

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
numeq | Specify the last **numeq** constraints to be equality constraints | 0
TS | Types of chordal extensions used in term sparsity iterations: "block"(maximal chordal extension), "MD" (approximately smallest chordal extension), false (invalidating term sparsity exploitation) | "block"
solver | Specify an SDP solver: "Mosek" or "COSMO" | "Mosek"
cosmo_setting | Parameters for the COSMO solver: cosmo_para(eps_abs, eps_rel, max_iter) | cosmo_para(1e-5, 1e-5, 1e4)
QUIET | Silence the output| false
solve | Solve the SDP relaxation | true
Gram | Output Gram matrices | false
constraint | nothing or "projection" (assume $x_i^2=x_i$ for all $i$) or "unipotent" (assume $x_i^2=1$ for all $i$) | nothing

### References
[1] [Optimization over trace polynomials](https://link.springer.com/article/10.1007/s00023-021-01095-4), 2021.  
