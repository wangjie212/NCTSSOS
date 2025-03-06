# This script shows how to solve complex state polynomial optimization problems using NCTSSOS.
# Examples below are taken from the paper: "State polynomials: positivity, optimization and nonlinear Bell inequalities"
using NCTSSOS
using LinearAlgebra

## Example 7.2.0
n = 4 # number of variables
# define the cost: ς^re(x1y1) + ς^re(x1y2) + ς^re(x2y1) − ς^re(x2y2)
supp = [
    [[[1; 3]], Vector{Int}[]],
    [[[1; 4]], Vector{Int}[]],
    [[[2; 3]], Vector{Int}[]],
    [[[2; 4]], Vector{Int}[]],
]
coe = [-1; -1; -1; 1]
d = 1 # relaxation order
opt, data = cpstateopt_first(
    supp, coe, n, d; vargroup=[2; 2], TS="block", constraint="unipotent"
)
# vargroup = [2;2] defines [xi, yj] = 0
# constraint = "unipotent" defines xi^2 = 1, yj^2 = 1

## Example 7.2.1
n = 4
supp = [
    [[[1; 4], [1; 4]], Vector{Int}[]],
    [[[2; 3], [2; 3]], Vector{Int}[]],
    [[[1; 4], [2; 3]], Vector{Int}[]],
    [[[1; 3], [1; 3]], Vector{Int}[]],
    [[[2; 4], [2; 4]], Vector{Int}[]],
    [[[1; 3], [2; 4]], Vector{Int}[]],
]
coe = [-1; -1; -2; -1; -1; 2]
d = 3
opt, data = cpstateopt_first(
    supp, coe, n, d; vargroup=[2; 2], TS="block", constraint="unipotent"
)

## Example 7.2.2
n = 6
supp = [
    [[[1; 4]], Vector{Int}[]],
    [[[1], [4]], Vector{Int}[]],
    [[[1; 5]], Vector{Int}[]],
    [[[1], [5]], Vector{Int}[]],
    [[[1; 6]], Vector{Int}[]],
    [[[1], [6]], Vector{Int}[]],
    [[[2; 4]], Vector{Int}[]],
    [[[2], [4]], Vector{Int}[]],
    [[[2; 5]], Vector{Int}[]],
    [[[2], [5]], Vector{Int}[]],
    [[[2; 6]], Vector{Int}[]],
    [[[2], [6]], Vector{Int}[]],
    [[[3; 4]], Vector{Int}[]],
    [[[3], [4]], Vector{Int}[]],
    [[[3; 5]], Vector{Int}[]],
    [[[3], [5]], Vector{Int}[]],
]
coe = [-1; 1; -1; 1; -1; 1; -1; 1; -1; 1; 1; -1; -1; 1; 1; -1]
d = 2
opt, data = cpstateopt_first(
    supp, coe, n, d; vargroup=[3; 3], TS="block", constraint="unipotent"
)

## Example 7.2.3
n = 4
supp = [
    [[[2]], Vector{Int}[]],
    [[[3]], Vector{Int}[]],
    [[[4]], Vector{Int}[]],
    [[[1; 3]], Vector{Int}[]],
    [[[2; 3]], Vector{Int}[]],
    [[[1; 4]], Vector{Int}[]],
    [[[2; 4]], Vector{Int}[]],
    [[[1], [3]], Vector{Int}[]],
    [[[2], [3]], Vector{Int}[]],
    [[[2], [4]], Vector{Int}[]],
    [[[1], [1]], Vector{Int}[]],
    [[[4], [4]], Vector{Int}[]],
]
coe = -[1; 1; 1; -1; 1; 1; 1; -1; -1; -1; -1; -1]
d = 2
opt, data = cpstateopt_first(
    supp, coe, n, d; vargroup=[2; 2], TS="block", constraint="unipotent"
)

# quantum bilocal networks
## Example 8.1.1
n = 6
supp = [
    [[[1; 3; 5], [1; 3; 5]], Vector{Int}[]],
    [[[1; 3; 6], [1; 3; 6]], Vector{Int}[]],
    [[[2; 3; 5], [2; 3; 5]], Vector{Int}[]],
    [[[2; 3; 6], [2; 3; 6]], Vector{Int}[]],
    [[[1; 4; 5], [1; 4; 5]], Vector{Int}[]],
    [[[1; 4; 6], [1; 4; 6]], Vector{Int}[]],
    [[[2; 4; 5], [2; 4; 5]], Vector{Int}[]],
    [[[2; 4; 6], [2; 4; 6]], Vector{Int}[]],
    [[[1; 3; 5], [1; 3; 6]], Vector{Int}[]],
    [[[1; 3; 5], [2; 3; 5]], Vector{Int}[]],
    [[[1; 3; 5], [2; 3; 6]], Vector{Int}[]],
    [[[1; 3; 5], [1; 4; 5]], Vector{Int}[]],
    [[[1; 3; 5], [1; 4; 6]], Vector{Int}[]],
    [[[1; 3; 5], [2; 4; 5]], Vector{Int}[]],
    [[[1; 3; 5], [2; 4; 6]], Vector{Int}[]],
    [[[1; 3; 6], [2; 3; 5]], Vector{Int}[]],
    [[[1; 3; 6], [2; 3; 6]], Vector{Int}[]],
    [[[1; 3; 6], [1; 4; 5]], Vector{Int}[]],
    [[[1; 3; 6], [1; 4; 6]], Vector{Int}[]],
    [[[1; 3; 6], [2; 4; 5]], Vector{Int}[]],
    [[[1; 3; 6], [2; 4; 6]], Vector{Int}[]],
    [[[2; 3; 5], [2; 3; 6]], Vector{Int}[]],
    [[[2; 3; 5], [1; 4; 5]], Vector{Int}[]],
    [[[2; 3; 5], [1; 4; 6]], Vector{Int}[]],
    [[[2; 3; 5], [2; 4; 5]], Vector{Int}[]],
    [[[2; 3; 5], [2; 4; 6]], Vector{Int}[]],
    [[[2; 3; 6], [1; 4; 5]], Vector{Int}[]],
    [[[2; 3; 6], [1; 4; 6]], Vector{Int}[]],
    [[[2; 3; 6], [2; 4; 5]], Vector{Int}[]],
    [[[2; 3; 6], [2; 4; 6]], Vector{Int}[]],
    [[[1; 4; 5], [1; 4; 6]], Vector{Int}[]],
    [[[1; 4; 5], [2; 4; 5]], Vector{Int}[]],
    [[[1; 4; 5], [2; 4; 6]], Vector{Int}[]],
    [[[1; 4; 6], [2; 4; 5]], Vector{Int}[]],
    [[[1; 4; 6], [2; 4; 6]], Vector{Int}[]],
    [[[2; 4; 5], [2; 4; 6]], Vector{Int}[]],
    [[[1; 3; 5]], Vector{Int}[]],
    [[[1; 3; 6]], Vector{Int}[]],
    [[[2; 3; 5]], Vector{Int}[]],
    [[[2; 3; 6]], Vector{Int}[]],
    [[[1; 4; 5]], Vector{Int}[]],
    [[[1; 4; 6]], Vector{Int}[]],
    [[[2; 4; 5]], Vector{Int}[]],
    [[[2; 4; 6]], Vector{Int}[]],
]
coe = [
    1 / 8 * [
        1
        1
        1
        1
        1
        1
        1
        1
        2
        2
        2
        -2
        2
        2
        -2
        2
        2
        -2
        2
        2
        -2
        2
        -2
        2
        2
        -2
        -2
        2
        2
        -2
        -2
        -2
        2
        2
        -2
        -2
    ]
    -[1; 1; 1; 1; 1; -1; -1; 1]
]
d = 3
@time begin
    opt, data = cpstateopt_first(
        supp,
        coe,
        n,
        d;
        constraint="unipotent",
        vargroup=[2; 2; 2],
        TS="block",
        solve=true,
        bilocal=[2; 5],
    )
end

## Example 8.1.2
n = 8
a = [[1; 4; 7], [1; 4; 8], [2; 4; 7], [2; 4; 8], [3; 4; 7], [3; 4; 8], [4; 7], [4; 8]]
ca = [1 / 2; 1 / 2; 1 / 2; 1 / 2; 1 / 2; 1 / 2; 1 / 2; 1 / 2]
b = [
    [1; 5; 7],
    [1; 5; 8],
    [2; 5; 7],
    [2; 5; 8],
    [3; 5; 7],
    [3; 5; 8],
    [5; 7],
    [5; 8],
    [1; 6; 7],
    [1; 6; 8],
    [2; 6; 7],
    [2; 6; 8],
]
cb = [
    1 / 2
    -1 / 2
    1 / 2
    -1 / 2
    -1 / 2
    1 / 2
    1 / 2
    -1 / 2
    1 / 2
    -1 / 2
    -1 / 2
    1 / 2
]
c = [[1], [2]]
cc = [1; 1]
supp = Vector{Vector{Vector{Vector{Int}}}}(undef, 275)
coe = zeros(275)
k = 1
for i in 1:8, j in 1:12
    supp[k] = [[a[i], b[j]], Vector{Int}[]]
    coe[k] = -2 * ca[i] * cb[j]
    k += 1
end
for i in 1:8, j in 1:2
    supp[k] = [[a[i], c[j]], Vector{Int}[]]
    coe[k] = -2 * ca[i] * cc[j]
    k += 1
end
for i in 1:12, j in 1:2
    supp[k] = [[b[i], c[j]], Vector{Int}[]]
    coe[k] = -2 * cb[i] * cc[j]
    k += 1
end
for i in 1:8
    supp[k:(k + 1)] = [[[a[i]], Vector{Int}[]], [[a[i], a[i]], Vector{Int}[]]]
    coe[k:(k + 1)] = [-8 * ca[i], ca[i]^2]
    k += 2
    for j in (i + 1):8
        supp[k] = [[a[i], a[j]], Vector{Int}[]]
        coe[k] = 2 * ca[i] * ca[j]
        k += 1
    end
end
for i in 1:12
    supp[k:(k + 1)] = [[[b[i]], Vector{Int}[]], [[b[i], b[i]], Vector{Int}[]]]
    coe[k:(k + 1)] = [-8 * cb[i], cb[i]^2]
    k += 2
    for j in (i + 1):12
        supp[k] = [[b[i], b[j]], Vector{Int}[]]
        coe[k] = 2 * cb[i] * cb[j]
        k += 1
    end
end
for i in 1:2
    supp[k:(k + 1)] = [[[c[i]], Vector{Int}[]], [[c[i], c[i]], Vector{Int}[]]]
    coe[k:(k + 1)] = [8 * cc[i], cc[i]^2]
    k += 2
    for j in (i + 1):2
        supp[k] = [[c[i], c[j]], Vector{Int}[]]
        coe[k] = 2 * cc[i] * cc[j]
        k += 1
    end
end
d = 3
@time begin
    opt, data = cpstateopt_first(
        supp,
        coe,
        n,
        d;
        constraint="unipotent",
        vargroup=[3; 3; 2],
        TS="block",
        solve=false,
        bilocal=[3; 7],
    )
    opt, data = cpstateopt_higher!(data; TS="block", bilocal=[3; 7])
end
println(-opt - 16)

## Example 8.1.3
n = 9
supp = [
    [[[4; 7]], Vector{Int}[]],
    [[[5; 8]], Vector{Int}[]],
    [[[6; 9]], Vector{Int}[]],
    [[[1; 4]], Vector{Int}[]],
    [[[2; 5]], Vector{Int}[]],
    [[[3; 6]], Vector{Int}[]],
    [[[1; 5; 9]], Vector{Int}[]],
    [[[1; 6; 8]], Vector{Int}[]],
    [[[2; 4; 9]], Vector{Int}[]],
    [[[2; 6; 7]], Vector{Int}[]],
    [[[3; 4; 8]], Vector{Int}[]],
    [[[3; 5; 7]], Vector{Int}[]],
]
coe = -[1 / 3; 1 / 3; 1 / 3; -1 / 3; -1 / 3; -1 / 3; -1; -1; -1; -1; -1; -1]
d = 3
@time begin
    opt, data = cpstateopt_first(
        supp,
        coe,
        n,
        d;
        constraint="unipotent",
        vargroup=[3; 3; 3],
        TS="block",
        solve=false,
        bilocal=[3; 7],
        zero_moments=true,
    )
    opt, data = cpstateopt_higher!(data; TS="block", bilocal=[3; 7], zero_moments=true)
end
# check flatness
k = 4
ind = [
    sum(length.(data.ptsupp[data.tbasis[1][data.wbasis[1][i][1]]])) +
    length(data.basis[1][data.wbasis[1][i][2]]) <= 1 for
    i in 1:length(data.wbasis[1][data.blocks[1][k]])
]
r0 = count(eigvals(data.moment[k][ind, ind]) .> 1e-3)
r1 = count(eigvals(data.moment[k]) .> 1e-3)
println([r0, r1])
