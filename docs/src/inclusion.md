# Inclusion constants for free spectrahedra

## One dichotomic and one arbitrary measurement

The following script allows to compute the bound for k=3 displayed in Table 2 from the paper [arXiv:2512.](https://arxiv.org/abs/2512.)

```Julia
using NCTSSOS
using DynamicPolynomials
using COSMO
settings = cosmo_para()
settings.eps_abs = 1e-7 # absolute residual tolerance
settings.eps_rel = 1e-7 # relative residual tolerance
settings.max_iter = 1e7 # maximum number of iterations

@ncpolyvar x[1:3] b[1:3]
obj = x[1]*b[1]+x[2]*b[2]+x[3]*b[3];
ineq = [1+1.5*b[2], 1+1.5*b[3], 1-1.5*(b[2]+b[3]), 1-(x[1]+4/3*x[2]-2/3*x[3]),1-(-x[1]+4/3*x[2]-2/3*x[3]),1-(x[1]-2/3*x[2]+4/3*x[3]),1-(-x[1]-2/3*x[2]+4/3*x[3]),1-(x[1]-2/3*x[2]-2/3*x[3]),1-(-x[1]-2/3*x[2]-2/3*x[3])];
eq = [b[1]^2-1];
eq2=[(b[2]+2/3)*(b[2]-4/3), (b[3]+2/3)*(b[3]-4/3), b[2]*b[3]-b[3]*b[2]];
pop=[-obj;ineq;eq;eq2];
opt,data = nctssos_first(pop, [x;b], 4, numeq=4, partition=3, TS=false, obj="eigen");
1/opt # 0.683012
```
