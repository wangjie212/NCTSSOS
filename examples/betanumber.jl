# using PyCall
# using Conda
using NCTSSOS
# Conda.add("networkx")
# nx = pyimport("networkx")

# graph6 code to igraph
# ig6 = "GUzrvW"
# edges = [[i+1,j+1] for (i,j) in nx.from_graph6_bytes(codeunits(ig6)).edges()]

# edges = [[1;2], [2;3], [3;4], [4;5], [1;5], [1;6], [2;6], [3;6], [4;6], [5;6]]
# edges = [[1,2],[1,4],[1,5],[1,7],[2,3],[2,4],[2,5],[3,4],[3,6],[4,5],[4,7],[5,6],[5,7],[6,7]]
edges = [[1, 2], [2, 3], [3, 4], [4, 5], [1, 5], [1, 6], [2, 6], [3, 6], [4, 6], [5, 6]]
sort!(edges)
n = 6
supp = Vector{Vector{Union{Vector{Vector{Int}},mixword}}}(undef, Int(n * (n - 1) / 2) + 1)
coe = Vector{Vector{Float64}}(undef, Int(n * (n - 1) / 2) + 1)
supp[1] = [[[i], [i]] for i in 1:n]
coe[1] = -ones(n)
k = 2
for i in 1:n, j in (i + 1):n
    supp[k] = [mixword([i; j], []), mixword([j; i], [])]
    if NCTSSOS.bfind(edges, length(edges), [i; j]) === nothing
        coe[k] = [1; -1]
    else
        coe[k] = [1; 1]
    end
    k += 1
end

d = 2 # relaxation order
@time opt, data = stateopt_first(
    supp,
    coe,
    n,
    d;
    numeq=Int(n * (n - 1) / 2),
    QUIET=true,
    TS="block",
    constraint="unipotent",
)
# @time opt,data = stateopt_higher!(data, QUIET=true, TS="block")
