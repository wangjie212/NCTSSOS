using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using NCTSSOS

n = 2
@ncpolyvar x[1:2]
f = 2-x[1]^2+x[1]*x[2]^2*x[1]-x[2]^2
g = 4-x[1]^2-x[2]^2
h = x[1]*x[2]+x[2]*x[1]-2
d = 2

# modelling with add_psatz!
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), false)
λ = @variable(model)
model,info1 = add_psatz!(model, f - λ, x, [g], [h], d, QUIET=true, TS=false, constrs="con1")
@objective(model, Max, λ)
optimize!(model)
status = termination_status(model)
if status != MOI.OPTIMAL
    println("termination status: $status")
    status = primal_status(model)
    println("solution status: $status")
end
objv = objective_value(model)
@show objv

# retrieve Gram matrices
GramMat = Vector{Vector{Vector{Union{Float64,Matrix{Float64}}}}}(undef, info1.cql)
for i = 1:info1.cql
    GramMat[i] = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, 1+length(info1.I[i])+length(info1.J[i]))
    for j = 1:1+length(info1.I[i])+length(info1.J[i])
        GramMat[i][j] = [value.(info1.gram[i][j][l]) for l = 1:info1.cl[i][j]]
    end
end

# retrieve moment matrices
moment = [-dual(constraint_by_name(model, "con1[$i]")) for i=1:length(info1.tsupp)]
MomMat = get_moment_matrix(moment, info1.tsupp, info1.cql, info1.basis)


# modelling with add_SOHS!
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), false)
λ = @variable(model)
s0 = add_SOHS!(model, x, 2)
s1 = add_SOHS!(model, x, 1)
p = add_poly!(model, x, 2)
@constraint(model, arrange(f - λ - s0 - s1*g - p*h, x)[2] .== 0)
@objective(model, Max, λ)
optimize!(model)
status = termination_status(model)
if status != MOI.OPTIMAL
    println("termination status: $status")
    status = primal_status(model)
    println("solution status: $status")
end
objv = objective_value(model)
@show objv
