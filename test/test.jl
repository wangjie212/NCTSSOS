using DynamicPolynomials
using NCTSSOS

n = 3
@ncpolyvar x[1:3]
f = x[1]^2-x[1]*x[2]-x[2]*x[1]+3x[2]^2-2x[1]*x[2]*x[1]+2x[1]*x[2]^2*x[1]-x[2]*x[3]-x[3]*x[2]+
6x[3]^2+9x[2]^2*x[3]+9x[3]*x[2]^2-54x[3]*x[2]*x[3]+142x[3]*x[2]^2*x[3]

opt,data = nctssos_first(f, x, newton=true, reducebasis=true, TS="MD", obj="eigen", QUIET=true)
# optimum = -0.0035512

opt,data = nctssos_first(f, x, newton=true, TS="MD", obj="trace", QUIET=true)
# optimum = -0.0035512

n = 2
@ncpolyvar x[1:2]
f = 2-x[1]^2+x[1]*x[2]^2*x[1]-x[2]^2
g = 4-x[1]^2-x[2]^2
h = x[1]*x[2]+x[2]*x[1]-2
pop = [f, g, h]

opt,data = nctssos_first(pop, x, 2, numeq=1, TS="MD", obj="eigen", QUIET=true)
# optimum = -0.9999999

opt,data = nctssos_higher!(data, TS="MD", QUIET=true)
# optimum = -0.9999999

opt,data = nctssos_first(pop, x, 2, numeq=1, TS="MD", obj="trace", QUIET=true)
# optimum = -0.9999999

opt,data = nctssos_higher!(data, TS="MD", QUIET=true)
# optimum = -0.9999999

n = 10
@ncpolyvar x[1:n]
f = 0.0
for i = 1:n
    jset = max(1,i-5):min(n,i+1)
    jset = setdiff(jset,i)
    f += (2x[i]+5*x[i]^3+1)^2
    f -= sum([4x[i]*x[j]+10x[i]^3*x[j]+2x[j]+4x[i]*x[j]^2+10x[i]^3*x[j]^2+2x[j]^2 for j in jset])
    f += sum([x[j]*x[k]+2x[j]^2*x[k]+x[j]^2*x[k]^2 for j in jset for k in jset])
end

opt,data = cs_nctssos_first([f], x, 3, TS="MD", obj="trace", QUIET=true)
# optimum = 6.5e-8

opt,data = cs_nctssos_higher!(data, TS="MD", QUIET=true)
# optimum = 6.5e-7

pop = [f]
for i = 1:n
    push!(pop, 1-x[i]^2)
    push!(pop, x[i]-1/3)
end

opt,data = cs_nctssos_first(pop, x, 3, TS="MD", obj="eigen", QUIET=true)
# optimum = 3.011288

opt,data = cs_nctssos_first(pop, x, 3, TS="MD", obj="trace", QUIET=true)
# optimum = 3.011288

println("Run successfully!")
