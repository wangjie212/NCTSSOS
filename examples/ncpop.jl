using DynamicPolynomials
using NCTSSOS

n = 3
@ncpolyvar x[1:3]
f = x[1]^2 - x[1]*x[2] - x[2]*x[1] + 3.0x[2]^2 - 2x[1]*x[2]*x[1] + 2x[1]*x[2]^2*x[1] - x[2]*x[3] - x[3]*x[2] +
6x[3]^2 + 9x[2]^2*x[3] + 9x[3]*x[2]^2 - 54x[3]*x[2]*x[3] + 142x[3]*x[2]^2*x[3]
pop = [f]
opt,data = ncpop(pop, x, 2, CS=false, TS="block", newton=true, Gram=true)

for i = 1:n
    push!(pop, 1 - x[i]^2)
    push!(pop, x[i] - 1/3)
end
opt,data = ncpop(pop, x, 2, TS="block", Gram=true)
opt,data = ncpop(data, TS="block", Gram=true)

@ncpolyvar x[1:2]
f = 3 + x[1]^2 + 2x[1]^3 + 2x[1]^4 + x[1]^6 - 4x[1]^4*x[2] + x[1]^4*x[2]^2 + 4x[1]^3*x[2] + 2x[1]^3*x[2]^2 -
2x[1]^3*x[2]^3 + 2x[1]^2*x[2] - x[1]^2*x[2]^2 + 8x[1]*x[2]*x[1]*x[2] + 2x[1]^2*x[2]^3 - 4x[1]*x[2] +
4x[1]*x[2]^2 + 6x[1]*x[2]^4 - 2x[2] + x[2]^2 - 4x[2]^3 + 2x[2]^4 + 2x[2]^6
opt,data = ncpop([f], x, 3, CS=false, TS=false, newton=false, Gram=true)

n = 2
@ncpolyvar x[1:2]
f = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2 + x[1]*x[2]*x[1]*x[2] + x[2]*x[1]*x[2]*x[1] +
x[1]^3*x[2] + x[2]*x[1]^3 + x[1]*x[2]^3 + x[2]^3*x[1]
g1 = 1 - x[1]^2
g2 = 1 - x[2]^2
pop = [f, g1, g2]
ncpop(pop, x, 2, TS="block")

f = x[2]*x[1] + x[2]*x[1]^2*x[2] - 2x[2]*x[1]^2 + x[2]*x[1]*x[2]*x[1] - 2x[2]*x[1]*x[2] +
x[1]*x[2] + x[1]*x[2]*x[1]*x[2] - 2x[1]^2*x[2] + x[2]*x[1]^2*x[2] - 2x[2]*x[1]*x[2]
pop = [f, 4-x[1]^2, 4-x[2]^2]
opt,data = ncpop(pop, x, 2, TS="block", obj="trace")

@ncpolyvar x[1:6]
f = x[1]*(x[4] + x[5] + x[6]) + x[2]*(x[4] + x[5] - x[6]) + x[3]*(x[4] - x[5]) - x[1] - 2x[4] - x[5]
ncpop([-f], x, 2, partition=3, constraint="projection", TS=false, QUIET=true, obj="eigen")
ncpop([-f], x, 5, soc=true, partition=3, constraint="projection", TS=false, QUIET=true, obj="eigen")

# Broyden banded polynomial
n = 100
@ncpolyvar x[1:n]
f = 0
for i = 1:n
    jset = max(1, i-5) : min(n, i+1)
    jset = setdiff(jset, i)
    g = sum(x[j] + x[j]^2 for j in jset)
    f += (2*x[i] + 5*x[i]^3 + 1 - g)^2
end

# Chained singular polynomial
n = 4
f = 0
@ncpolyvar x[1:n]
for i = 1:2:n-3
    f += (x[i] + 10*x[i+1])^2 + 5*(x[i+2] - x[i+3])^2 + (x[i+1] - 2*x[i+2])^4 + 10*(x[i] - 10*x[i+3])^4
end
