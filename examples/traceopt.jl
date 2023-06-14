using DelimitedFiles
using NCTSSOS

# The toy example
n = 3
supp = [[[1;2;3]], [[1;2], [3]]]
coe = [1; 1]
d = 3
opt,data = ptraceopt_first(supp, coe, n, d, TS=false, constraint="projection")

# Example 6.2.0
n = 4
supp = [[[1;3]], [[1;4]], [[2;3]], [[2;4]]]
coe = [-1; -1; -1; 1]
d = 1
opt,data = ptraceopt_first(supp, coe, n, d, TS="block", constraint="nilpotent")
# opt,data = ptraceopt_higher!(data, TS="block")

# Example 6.2.1
n = 4
supp = [[[1;4], [1;4]], [[2;3], [2;3]], [[1;4], [2;3]], [[1;3], [1;3]], [[2;4], [2;4]], [[1;3], [2;4]]]
coe = [-1; -1; -2; -1; -1; 2]
d = 2
opt,data = ptraceopt_first(supp, coe, n, d, TS=false, constraint="nilpotent")
# opt,data = ptraceopt_higher!(data, TS="block")

# Example 6.2.2
n = 6
supp = [[[1;4]], [[1], [4]], [[1;5]], [[1], [5]], [[1;6]], [[1], [6]], [[2;4]], [[2], [4]],
[[2;5]], [[2], [5]], [[2;6]], [[2], [6]], [[3;4]], [[3], [4]], [[3;5]], [[3], [5]]]
coe = [-1; 1; -1; 1; -1; 1; -1; 1; -1; 1; 1; -1; -1; 1; 1; -1]
d = 2
opt,data = ptraceopt_first(supp, coe, n, d, TS=false, constraint="nilpotent")
# opt,data = ptraceopt_higher!(data, TS="block")

# Example 6.2.3
n = 8
supp = [[[1;3], [1;3], [5;7], [5;7]], [[2;3], [2;3], [5;7], [5;7]], [[1;3], [1;3], [5;8], [5;8]], [[2;3], [2;3], [5;8], [5;8]],
[[1;4], [1;4], [6;7], [6;7]], [[2;4], [2;4], [6;7], [6;7]], [[1;4], [1;4], [6;8], [6;8]], [[2;4], [2;4], [6;8], [6;8]],
[[1;3], [2;3], [5;7], [5;7]], [[1;3], [1;3], [5;7], [5;8]], [[1;3], [2;3], [5;7], [5;8]], [[1;3], [1;4], [5;7], [6;7]],
[[1;3], [2;4], [5;7], [6;7]], [[1;3], [1;4], [5;7], [6;8]], [[1;3], [2;4], [5;7], [6;8]],
[[2;3], [2;3], [5;7], [5;8]], [[2;3], [1;4], [5;7], [6;7]], [[2;3], [2;4], [5;7], [6;7]], [[2;3], [1;4], [5;7], [6;8]],
[[2;3], [2;4], [5;7], [6;8]], [[1;3], [2;3], [5;8], [5;8]], [[1;3], [1;4], [5;8], [6;7]], [[1;3], [2;4], [5;8], [6;7]],
[[1;3], [1;4], [5;8], [6;8]], [[1;3], [2;4], [5;8], [6;8]], [[2;3], [1;4], [5;8], [6;7]], [[2;3], [2;4], [5;8], [6;7]],
[[2;3], [1;4], [5;8], [6;8]], [[2;3], [2;4], [5;8], [6;8]], [[1;4], [2;4], [6;7], [6;7]], [[1;4], [1;4], [6;7], [6;8]],
[[1;4], [2;4], [6;7], [6;8]], [[2;4], [2;4], [6;7], [6;8]], [[1;4], [2;4], [6;8], [6;8]],
[[1;3], [5;7]], [[2;3], [5;7]], [[1;3], [5;8]], [[2;3], [5;8]], [[1;4], [6;7]], [[2;4], [6;7]], [[1;4], [6;8]], [[2;4], [6;8]]]
coe = -[-1/8*[1; 1; 1; 1; 1; 1; 1; 1; 2; 2; 4; 2; -2; -2; 2; 2; 2; -2; -2; 2; 2; 2; -2; -2; 2; 2; -2; -2; 2; -2; -2; 4; -2; -2];
1; 1; 1; 1; 1; -1; -1; 1]
d = 4
opt,data = ptraceopt_first(supp, coe, n, d, TS="block", constraint="projection")

# Werner witness example
n = 4
a = [0 -1 1 1 1 -1 1 1 0 0 -1 1 0 1 -1 -1]
Y = 1/12*a'*a
io = open("D:\\Programs\\NCTSSOS\\examples\\Werner_matrix.csv", "r")
A = readdlm(io, ',')
close(io)
io = open("D:\\Programs\\NCTSSOS\\examples\\Werner_sigma.csv", "r")
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
opt,data = Werner_witness_first(dY, c, n, d, monosquare=true, TS=false)
opt,data = Werner_witness_higher!(data, TS="block")

n = 4
supp = Vector{Vector{Union{Vector{Vector{Int}}, mixword}}}(undef, 3)
supp[1] = [[[1;4], [1;4]], [[2;3], [2;3]], [[1;4], [2;3]], [[1;3], [1;3]], [[2;4], [2;4]], [[1;3], [2;4]]]
supp[2] = [mixword([1;2], []), mixword([3;4], []), mixword([], [])]
supp[3] = [mixword([], [[1], [2]]), mixword([], [[3], [4]]), mixword([], [])]
coe = [[-1; -1; -2; -1; -1; 2], [1; 1; -1], [1; 1; -1]]
d = 2
opt,data = traceopt_first(supp, coe, n, d, numeq=1, TS=false, solve=true, constraint="nilpotent")
