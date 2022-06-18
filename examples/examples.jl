using DynamicPolynomials
using DelimitedFiles
using Plots

include("D:\\Programs\\NCTSSOS\\src\\NCTSSOS.jl")
using .NCTSSOS

n = 3
@ncpolyvar x[1:3]
f = x[1]^2-x[1]*x[2]-x[2]*x[1]+3x[2]^2-2x[1]*x[2]*x[1]+2x[1]*x[2]^2*x[1]-x[2]*x[3]-x[3]*x[2]+
6x[3]^2+9x[2]^2*x[3]+9x[3]*x[2]^2-54x[3]*x[2]*x[3]+142x[3]*x[2]^2*x[3]

@time begin
opt,data = nctssos_first(f, x, newton=true, reducebasis=false, TS=false, QUIET=true, obj="eigen")
end

f = 3+x[1]^2+2x[1]^3+2x[1]^4+x[1]^6-4x[1]^4*x[2]+x[1]^4*x[2]^2+4x[1]^3*x[2]+2x[1]^3*x[2]^2-
2x[1]^3*x[2]^3+2x[1]^2*x[2]-x[1]^2*x[2]^2+8x[1]*x[2]*x[1]*x[2]+2x[1]^2*x[2]^3-4x[1]*x[2]+
4x[1]*x[2]^2+6x[1]*x[2]^4-2x[2]+x[2]^2-4x[2]^3+2x[2]^4+2x[2]^6

n = 10
supp = [UInt16[]]
coe = Float64[n]
for i = 2:n
    push!(supp, UInt16[i-1;i-1;i-1;i-1], UInt16[i-1;i-1;i], UInt16[i], UInt16[i;i])
    push!(coe, 100, -200, -2, 101)
end

@time begin
opt,data = nctssos_first(supp, coe, n, newton=true, QUIET=true, reducebasis=false, TS="MD", obj="eigen")
end

# @time begin
# opt,data=nctssos_first(supp,coe,newton=true,reducebasis=false,TS="block",obj="eigen")
# end
@time begin
opt,data = nctssos_first(supp,coe,n,2,dg,reducebasis=false,TS="MD",obj="eigen")
end

@time begin
opt,data = nctssos_higher!(data,TS="MD")
end

n = 2
@ncpolyvar x[1:2]
f = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2 + x[1]*x[2]*x[1]*x[2] + x[2]*x[1]*x[2]*x[1] +
x[1]^3*x[2] + x[2]*x[1]^3 + x[1]*x[2]^3 + x[2]^3*x[1]
g1 = 1 - x[1]^2
g2 = 1 - x[2]^2
# f = x[1]*x[2]*x[1]
# f = x[1]*x[2]^4*x[1]+x[2]*x[1]^4*x[2]-3x[1]*x[2]^2*x[1]+1
# g1 = 1-x[1]^2-x[2]^2
pop = [f, g1, g2]

# f = (1-x[1]^2)*(1-x[2]^2)+(1-x[2]^2)*(1-x[1]^2)
f = x[2]*x[1]+x[2]*x[1]^2*x[2]-2x[2]*x[1]^2+x[2]*x[1]*x[2]*x[1]-2x[2]*x[1]*x[2]+
x[1]*x[2]+x[1]*x[2]*x[1]*x[2]-2x[1]^2*x[2]+x[2]*x[1]^2*x[2]-2x[2]*x[1]*x[2]
pop = [f, 4-x[1]^2, 4-x[2]^2]

opt,data = nctssos_first(pop, x, 2, TS="MD", obj="trace")
opt,data = nctssos_higher!(data, TS="MD")

l = n
# b = 15
# n = (b-5)*l+5
# n = b
supp = Vector{Vector{Vector{UInt16}}}(undef, 2l+1)
coe = Vector{Vector{Float64}}(undef, 2l+1)
# io = open("D:\\Programs\\NCTSSOS\\examples\\random\\supp_$n.txt", "r")
# temp = readdlm(io, UInt16)
# close(io)
# supp[1] = [temp[i,:] for i in 1:size(temp,1)]
# io = open("D:\\Programs\\NCTSSOS\\examples\\random\\coe_$n.txt", "r")
# temp = readdlm(io, Float64)
# close(io)
# coe[1] = temp[:,1]
# supp[1] = Vector{Vector{UInt16}}[]
# for i = 1:l, j = 1:4b
#     push!(supp[1], ceil.(UInt16, ((b-5)*i-(b-5)).+rand(4).*b))
# end
# coe[1] = 2*rand(4b*l).-1
for i = 1:l
    supp[i+1] = [UInt16[], UInt16[i;i]]
    coe[i+1] = Float64[1;-1]
    supp[i+n+1] = [UInt16[], UInt16[i]]
    coe[i+n+1] = Float64[-1/3;1]
    # supp[i+1] = [[UInt16[]]; [UInt16[j;j] for j=(b-5)*i-(b-6):(b-5)*i+5]]
    # coe[i+1] = Float64[1; -ones(b)]
end

@time begin
opt,data = cs_nctssos_first(supp, coe, n, 3, TS="MD", QUIET=true, obj="trace")
end

@time begin
opt,data = nctssos_first(supp, coe, n, 2, TS="MD", obj="eigen")
end

@time begin
opt,data = cs_nctssos_higher!(data, QUIET=false, TS="MD")
end

io = open("D:\\Programs\\NCTSSOS\\examples\\supp_$n.txt", "w")
writedlm(io, convert(Vector{Vector{Int}}, supp[1]), ' ')
close(io)
io = open("D:\\Programs\\NCTSSOS\\examples\\coe_$n.txt", "w")
writedlm(io, coe[1], ' ')
close(io)

@time begin
opt,data = cs_nctssos_first(supp, coe, n, 2, TS="MD", obj="eigen")
end

@time begin
opt,data = cs_nctssos_higher!(data, TS="MD")
end

# @time begin
# opt,data = nctssos_first(pop,[x;y],2,reducebasis=false,numeq=length(pop)-10,TS=false,obj="eigen")
# end
@time begin
opt,data = nctssos_first(supp,coe,n,2,reducebasis=false,TS=false,obj="trace")
end
@time begin
opt,data = nctssos_higher!(data,TS="MD")
end
@time begin
opt,data = nctssos_first(supp,coe,n,2,reducebasis=false,TS="MD",obj="eigen")
end
