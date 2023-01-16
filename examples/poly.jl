# the Broyden banded polynomial
n = 100
@ncpolyvar x[1:n]
f = 0.
for i = 1:n
    jset = max(1, i-5) : min(n, i+1)
    jset = setdiff(jset, i)
    g = sum(x[j] + x[j]^2 for j in jset)
    global f += (2*x[i] + 5*x[i]^3 + 1 - g)^2
end

# the Broyden banded polynomial
n = 4
supp = [UInt16[]]
coe = Float64[n]
for i = 1:n
    jset = max(1,i-5):min(n,i+1)
    jset = setdiff(jset,i)
    push!(supp, [i], [i;i], [i;i;i], [i;i;i;i], [i;i;i;i;i;i])
    push!(coe, 4, 4, 10, 20, 25)
    for j in jset
        push!(supp, [j], [j;j], [i;j], [i;j;j], [i;i;i;j], [i;i;i;j;j])
        push!(coe, -2, -2, -4, -4, -10, -10)
        for k in jset
            push!(supp, [j;k], [j;j;k], [j;j;k;k])
            push!(coe, 1, 2, 1)
        end
    end
end

# the generalized Rosenbrock polynomial
n = 4000
supp[1] = [UInt16[]]
coe[1] = Float64[n]
for i = 2:n
    push!(supp[1], [i-1;i-1;i-1;i-1], [i-1;i-1;i], [i], [i;i])
    push!(coe[1], 100, -200, -2, 101)
end

# the chained Wood polynomial
n = 600
supp[1] = [UInt16[]]
coe[1] = Float64[21*n-41]
for i = 1:2:n-3
    push!(supp[1], [i], [i;i], [i;i;i;i], [i;i;i+1], [i+1], [i+1;i+1], [i+1;i+3],
    [i+2], [i+2;i+2], [i+2;i+2;i+2;i+2], [i+2;i+2;i+3], [i+3], [i+3;i+3])
    push!(coe[1], -2, 1, 100, -200, -40, 110.1, 19.8, -2, 1, 90, -180, -40, 100.1)
end

# the Broyden tridiagonal polynomial
n = 20
supp[1] = [UInt16[]]
coe[1] = Float64[n]
push!(supp[1], [1;1], [1;1;1;1], [2;2], [1;1;1], [1;2], [1], [1;1;2], [2])
push!(coe[1], 5, 4, 4, -12, -12, 6, 8, -4)
for i = 2:n-1
    push!(supp[1], [i;i], [i;i;i;i], [i-1;i-1], [i+1;i+1], [i;i;i], [i-1;i],
    [i;i+1], [i], [i-1;i;i], [i;i;i+1], [i-1;i+1], [i-1], [i+1])
    push!(coe[1], 5, 4, 1, 4, -12, -6, -12, 6, 4, 8, 4, -2, -4)
end
push!(supp[1], [n;n], [n;n;n;n], [n-1;n-1], [n;n;n], [n-1;n], [n], [n-1;n;n], [n-1])
push!(coe[1], 5, 4, 1, -12, -6, 6, 4, -2)

# the Chained singular polynomial
n = 4
f = 0
@ncpolyvar x[1:n]
for i = 1:2:n-3
    global f += (x[i] + 10*x[i+1])^2 + 5*(x[i+2] - x[i+3])^2 + (x[i+1] - 2*x[i+2])^4 + 10*(x[i] - 10*x[i+3])^4
end

# the Chained singular polynomial
n = 80
supp = Vector{UInt16}[]
coe = Float64[]
for i = 1:2:n-3
    push!(supp, [i;i;i;i], [i;i;i;i+3], [i;i;i+3;i+3], [i;i;i+3;i], [i;i+3;i+3;i+3], [i;i+3;i+3;i], [i;i+3;i;i],
    [i;i+3;i;i+3], [i+1;i+1;i+1;i+1], [i+1;i+1;i+1;i+2], [i+1;i+1;i+2;i+2], [i+1;i+1;i+2;i+1], [i+1;i+2;i+2;i+2],
    [i+1;i+2;i+2;i+1], [i+1;i+2;i+1;i+1], [i+1;i+2;i+1;i+2], [i+2;i+2;i+2;i+2], [i+2;i+2;i+2;i+1], [i+2;i+2;i+1;i+1],
    [i+2;i+2;i+1;i+2], [i+2;i+1;i+1;i+1], [i+2;i+1;i+1;i+2], [i+2;i+1;i+2;i+2], [i+2;i+1;i+2;i+1],
    [i+3;i+3;i+3;i+3], [i+3;i+3;i+3;i], [i+3;i+3;i;i], [i+3;i+3;i;i+3], [i+3;i;i;i], [i+3;i;i;i+3], [i+3;i;i+3;i+3],
    [i+3;i;i+3;i], [i;i], [i;i+1], [i+1;i+1], [i+1;i], [i+2;i+2], [i+2;i+3], [i+3;i+3], [i+3;i+2])
    append!(coe, [10;-100;1000;-100;-10000;1000;-100;1000;1;-2;4;-2;-8;4;-2;4;16;-8;4;-8;-2;4;-8;4;100000;-10000;
    1000;-10000;-100;1000;-10000;1000;1;10;100;10;5;-5;5;-5])
end

supp[1] = [UInt16[]]
coe[1] = Float64[0]
for i = 1:2:n-3
    push!(supp[1], [i;i], [i;i+1], [i+1;i+1], [i+2;i+2], [i+3;i+3], [i+2;i+3], [i+1;i+1;i+1;i+1], [i+2;i+1;i+1;i+2], [i+2;i+2;i+2;i+2], [i+1;i+1;i+1;i+2], [i+1;i+1;i+2;i+2], [i+2;i+1;i+2;i+2], [i;i;i;i], [i+3;i;i;i+3], [i+3;i+3;i+3;i+3], [i;i;i;i+3], [i;i;i+3;i+3], [i+3;i;i+3;i+3])
    push!(coe[1], 1, 20 ,100, 5, 5, -10, 1, 16, 16, -8, 8, -32, 10, 40, 10, -40, 20, -40)
end

# randomly generated examples
l = 2
b = 15
n = (b-5)*l+5
supp = Vector{Vector{Vector{UInt16}}}(undef, l+1)
coe = Vector{Vector{Float64}}(undef, l+1)
supp[1] = Vector{Vector{UInt16}}[]
for i = 1:l, j = 1:b
    push!(supp[1], ceil.(UInt16, ((b-5)*i-(b-5)).+rand(4).*b))
end
coe[1] = 2*rand(b*l).-1
dg = Vector{Int}(undef, l)
for i = 1:l
    dg[i] = 2
    supp[i+1] = [[UInt16[]]; [UInt16[j;j] for j=(b-5)*i-(b-6):(b-5)*i+5]]
    coe[i+1] = Float64[1; -ones(b)]
end

# io = open("E:\\Programs\\NCTSSOS\\examples\\supp_$n.txt", "r")
# temp = readdlm(io, UInt16)
# close(io)
# supp[1] = [temp[i,:] for i in 1:size(temp,1)]
# io = open("E:\\Programs\\NCTSSOS\\examples\\coe_$n.txt", "r")
# temp = readdlm(io, Float64)
# close(io)
# coe[1] = temp[:,1]

dg = Vector{Int}(undef, 2*n)
for i = 1:n
    dg[i] = 2
    dg[i+n] = 1
    supp[i+1] = [UInt16[], UInt16[i;i]]
    coe[i+1] = Float64[1, -1]
    supp[i+n+1] = [UInt16[], UInt16[i]]
    coe[i+n+1] = Float64[-1/3, 1]
end

# Bell inequality
a = [-2 0 0 0]
b = [0 0 -1 -1 -1]
C = [1 1 1 1 1; -1 0 -1 1 0; -1 -1 0 0 1; 0 -1 1 0 0]
nx,ny = length(a),length(b)
supp = Vector{UInt16}[]
coe = Float64[]
for i = 1:nx
    if a[i] != 0
        push!(supp, UInt16[i])
        push!(coe, -a[i])
    end
end
for j = 1:ny
    if b[j] != 0
        push!(supp, UInt16[nx+j])
        push!(coe, -b[j])
    end
end
for i = 1:nx, j = 1:ny
    if C[i,j] != 0
        push!(supp, UInt16[i; nx+j])
        push!(coe, -C[i,j])
    end
end

@time begin
bell(supp,coe,nx+ny,nx,2,CS=false,TS=false)
end

yvar=UInt16[nx+i for i=1:ny]
cliques=Vector{Vector{UInt16}}(undef, nx)
for i=1:nx
    cliques[i]=[i;yvar]
end
cliquesize=(ny+1)*ones(Int, nx)
cql=nx

function bell(supp, coe, n, nx, d; CS="MD", TS=false, QUIET=false)
    cliques,cql,cliquesize=clique_decomp(n,supp,alg=CS)
    blocks,cl,blocksize,_,_,basis,_=get_blocks_mix(d,supp,cliques,cql,cliquesize,TS=TS,nx=nx)
    opt,_=blockupop_mix(n,d,supp,coe,basis,cliques,cql,cliquesize,blocks,cl,blocksize,nx=nx)
    return -opt
end

supp = Vector{UInt16}[[1;4], [1;5], [1;6], [2;4], [2;5], [2;6], [3;4], [3;5], [1], [4], [5]]
coe = -[1, 1, 1, 1, 1, -1, 1, -1, -1, -2, -1]
n = 6
supp = Vector{UInt16}[[1;4], [1;5], [1;6], [2;4], [2;5], [2;6], [3;4], [3;5], [1], [2], [4], [5], []]
coe = -[1/4, 1/4, 1/4, 1/4, 1/4, -1/4, 1/4, -1/4, 1/4, 1/4, -1/4, -1/4, -1]
n = 6
opt,data = nctssos_first(supp, coe, n, d=3, TS=false, monosquare=false, partition=3, constraint="projection")
opt,data = nctssos_higher!(data, TS="block")
