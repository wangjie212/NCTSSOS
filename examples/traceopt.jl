using DelimitedFiles
include("E:\\Programs\\NCTSSOS\\src\\NCTSSOS.jl")
using .NCTSSOS

# The toy example
n = 3
tr_supp = [[[1;2;3]], [[1;2], [3]]]
coe = [1; 1]
d = 2
opt,data = ptraceopt_first(tr_supp, coe, n, d, TS=false, constraint="projection")

# Example 6.2.0
n = 4
tr_supp = [[[1;3]], [[1;4]], [[2;3]], [[2;4]]]
coe = [-1; -1; -1; 1]
d = 1
opt,data = ptraceopt_first(tr_supp, coe, n, d, TS="block", constraint="unipotent")
# opt,data = ptraceopt_higher!(data, TS="block")

# Example 6.2.1
n = 4
tr_supp = [[[1;4], [1;4]], [[2;3], [2;3]], [[1;4], [2;3]], [[1;3], [1;3]], [[2;4], [2;4]], [[1;3], [2;4]]]
coe = [-1; -1; -2; -1; -1; 2]
d = 2
opt,data = ptraceopt_first(tr_supp, coe, n, d, TS="block", constraint="unipotent")
# opt,data = ptraceopt_higher!(data, TS="block")

# Example 6.2.2
n = 6
tr_supp = [[[1;4]], [[1], [4]], [[1;5]], [[1], [5]], [[1;6]], [[1], [6]], [[2;4]], [[2], [4]],
[[2;5]], [[2], [5]], [[2;6]], [[2], [6]], [[3;4]], [[3], [4]], [[3;5]], [[3], [5]]]
coe = [-1; 1; -1; 1; -1; 1; -1; 1; -1; 1; 1; -1; -1; 1; 1; -1]
d = 2
opt,data = ptraceopt_first(tr_supp, coe, n, d, TS=false, constraint="unipotent")
# opt,data = ptraceopt_higher!(data, TS="block")

# Werner witness example
n = 4
a = [0 -1 1 1 1 -1 1 1 0 0 -1 1 0 1 -1 -1]
Y = 1/12*a'*a
io = open("E:\\Programs\\NCTSSOS\\examples\\Werner_matrix.csv", "r")
A = readdlm(io, ',')
close(io)
io = open("E:\\Programs\\NCTSSOS\\examples\\Werner_sigma.csv", "r")
c = readdlm(io, ',')
close(io)
dY = Vector{Vector{Array{Float64, 2}}}(undef, 52)
for i = 1:52
    dY[i] = Vector{Array{Float64, 2}}(undef, 4)
    for j = 1:4
        dY[i][5-j] = A[2j-1:2j, 2i-1:2i]
    end
end

d = 2
opt,data = Werner_witness_first(dY, c, n, d, TS=false)
opt,data = Werner_witness_higher!(data, TS="block")

sY = sum(c[i]*kron(kron(kron(dY[i][1], dY[i][2]), dY[i][3]), dY[i][4]) for i=1:52)
opnorm(sY-Y, Inf)
