# Bell inequalities

## Linear bell inequalities

The CHSH inequality:

```Julia
using DynamicPolynomials
using NCTSSOS
@ncpolyvar x[1:2]
@ncpolyvar y[1:2]
f = x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2]
opt,data = nctssos_first([-f], [x;y], 1, TS=false, partition=2, constraint="unipotent")
```

The I_3322 inequality:

```Julia
using DynamicPolynomials
using NCTSSOS
@ncpolyvar x[1:3]
@ncpolyvar y[1:3]
f = x[1]*(y[1] + y[2] + y[3]) + x[2]*(y[1] + y[2] - y[3]) + x[3]*(y[1] - y[2]) - x[1] - 2*y[1] - y[2]
opt,data = nctssos_first([-f], [x;y], 3, TS=false, normality=true, partition=3, constraint="projection")
```

The I_3322 inequality solved with the sequential NPA hierarchy:

```Julia
using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using NCTSSOS
@ncpolyvar y[1:3]
y1=y[1]; y2=y[2]; y3=y[3];

### f= -1 + x1/4 + x2/4 - y1/4 - y2/4 + x1*y1/4 + x1*y2/4 - x1*y3/4 + x2*y1/4 + x2*y2/4 + x2*y3/4 - x3*y1/4 + x3*y2/4
### Objective function of I_3322 after the change of variable x -> 2x-1, y-> 2y-1 to optimize over unitaries

### The 3 observables x are not needed anymore.
### Instead, one defines 6 linear functionals L_ij satisfying
### L_0j + L_1j = L, for all j=1,2,3
### L_j := L_0j - L_1j, for all j=1,2,3
### L_j(yk) = <xj yk> = L_0j(yk) - L_1j(yk)
### L(1) = <1> = L_01(1) + L_11(1)
### L(yk) = <yk> = L_0j(yk) + L_1j(yk)
### L_j(1) = <xj> = L_0j(1) - L_1j(1)
### Then one can write <f> = sum_ij L_ij (fij), with fij given below:

f01 = -1 + y1/4 + y2/4 - y3/4 + 1/4 - y1/4 ;
f11 = -1 + -y1/4 - y2/4 + y3/4 - 1/4 - y1/4;
f02 = y1/4 + y2/4 + y3/4 + 1/4 - y2/4;
f12 = -y1/4 - y2/4 - y3/4 - 1/4 - y2/4 ;
f03 = -y1/4 + y2/4;
f13 = y1/4 - y2/4;

d = 5;
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), false)
λ = @variable(model)
p = add_poly!(model,y,2*d,constraint="unipotent");
q = add_poly!(model,y,2*d,constraint="unipotent");

info = add_psatz!(model, λ-f01+p,y,[], [],d, QUIET=true, TS=false, constraint="unipotent");
info = add_psatz!(model,λ-f11+p,y,[],[],d, QUIET=true, TS=false, constraint="unipotent");
info = add_psatz!(model,-f02-p+q,y,[],[],d, QUIET=true, TS=false, constraint="unipotent");
info = add_psatz!(model,-f12-p+q,y,[],[],d, QUIET=true, TS=false, constraint="unipotent");
info = add_psatz!(model,-f03-q,y,[],[],d, QUIET=true, TS=false, constraint="unipotent");
info = add_psatz!(model,-f13-q,y,[],[],d, QUIET=true, TS=false, constraint="unipotent");
@objective(model, Min, λ)
optimize!(model)
objv = objective_value(model)
@show objv
```

## Nonlinear bell inequalities
