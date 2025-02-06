# State Polynomial Optimization

State polynomial optimization concerns minimizing a state polynomial subject to a tuple of noncommutative state polynomial inequality constraints and equality constraints, which in general takes the form:

$$\mathrm{inf}_{\mathbf{x}\in\mathcal{B}(\mathcal{H})^n, \varsigma\in\mathcal{S}(\mathcal{H})^n}\ f(\mathbf{x}; \varsigma)\ \text{ s.t. }\ g_1(\mathbf{x}; \varsigma)\ge0,\ldots,g_m(\mathbf{x}; \varsigma)\ge0,h_1(\mathbf{x}; \varsigma)=0,\ldots,h_{\ell}(\mathbf{x}; \varsigma)=0,$$

where $f$ is a state polynomial in noncommuting variables $\mathbf{x}$, $g_1,\ldots,g_m,h_1,\ldots,h_{\ell}$ are noncommutative state polynomials in noncommuting variables $\mathbf{x}$, and $\mathcal{H}$ is an (infinite dimensional) seperable Hilbert space.

To illustrate how to solve a state polynomial optimization problem with NCTSSOS, let us consider the following simple example.

```Julia
using NCTSSOS
n = 4 # number of variables
# define the objective function: ς(x1y1) + ς(x1y2) + ς(x2y1) − ς(x2y2)
supp = [[[1;3]], [[1;4]], [[2;3]], [[2;4]]] # define the coefficient data
coe = [1; 1; 1; -1] # define the coefficient data
d = 1 # set the relaxation order
opt,data = pstateopt_first(supp, coe, n, d, vargroup=[2;2], constraint="unipotent") # vargroup=[2;2] constrains [xi, yj] = 0
opt,data = pstateopt_higher!(data) # compute higher TS steps of the NCTSSOS hierarchy
```

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
numeq | Specify the last **numeq** constraints to be equality constraints | 0
TS | Types of chordal extensions used in term sparsity iterations: "block"(maximal chordal extension), "MD" (approximately smallest chordal extension), false (invalidating term sparsity exploitation) | "block"
vargroup | Partition the variables into several groups such that variables from different groups commute | [n]
solver | Specify an SDP solver: "Mosek" or "COSMO" | "Mosek"
cosmo_setting | Parameters for the COSMO solver: cosmo_para(eps_abs, eps_rel, max_iter) | cosmo_para(1e-5, 1e-5, 1e4)
QUIET | Silence the output| false
solve | Solve the SDP relaxation | true
Gram | Output Gram matrices | false
constraint | nothing or "projection" (assume $x_i^2=x_i$ for all $i$) or "unipotent" (assume $x_i^2=1$ for all $i$) | nothing

### References
[1] [State polynomials: positivity, optimization and nonlinear Bell inequalities](https://arxiv.org/abs/2301.12513), 2023. 
