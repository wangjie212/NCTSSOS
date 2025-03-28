using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using NCTSSOS

n = 2
@ncpolyvar x[1:2]
f = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
g = 4 - x[1]^2 - x[2]^2
h = x[1]*x[2] + x[2]*x[1] - 2
d = 2

# modelling with nctssos
opt,data = nctssos_first([f;g;h], x, d, numeq=1, TS=false, Gram=true)

# modelling with add_psatz!
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), false)
λ = @variable(model)
info = add_psatz!(model, f - λ, x, [g], [h], d, QUIET=true, TS=false, constrs="con")
@objective(model, Max, λ)
optimize!(model)
objv = objective_value(model)
@show objv

# retrieve Gram matrices
GramMat = Vector{Vector{Vector{Union{Float64,Matrix{Float64}}}}}(undef, info.cql)
for i = 1:info.cql
    GramMat[i] = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, 1+length(info.I[i]))
    for j = 1:1+length(info.I[i])
        GramMat[i][j] = [value.(info.gram[i][j][l]) for l = 1:info.cl[i][j]]
    end
end

# retrieve moment matrices
moment = -dual(constraint_by_name(model, "con"))
MomMat = get_moment_matrix(moment, info)


# modelling with add_SOHS!
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), false)
λ = @variable(model)
s0 = add_SOHS!(model, x, 2) # generate an unknown sohs of degree 2d
s1 = add_SOHS!(model, x, 1)
p = add_poly!(model, x, 2) # generate an unknown nc polynomial of degree d
@constraint(model, arrange(f - λ - s0 - s1*g - p*h - h*star(p), x)[2] .== 0)
@objective(model, Max, λ)
optimize!(model)
objv = objective_value(model)
@show objv
