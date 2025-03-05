using Test, NCTSSOS, DynamicPolynomials, Clarabel
using JuMP

using NCTSSOS: get_C_α_j
using SparseArrays


# FIXME: This is not what it does 
# @testset "C_α_j" begin
# 	model = Model()
# 	@variable(model, x[1:4])

#     cons = @constraint(model,[x[1] - x[2] x[3] x[4] + x[1]; x[1] - x[2] x[3] x[4] + x[1]; x[1] - x[2] x[3] x[4] + x[1]] in PSDCone())

#     C_α_js = get_C_α_j(x, constraint_object(cons))

#     @test C_α_js == [sparse([1, 2, 3, 1, 2, 3], [1, 1, 1, 3, 3, 3], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 3, 3),
# 	 sparse([1, 2, 3], [1, 1, 1], [-1.0, -1.0, -1.0], 3, 3),
# 	 sparse([1, 2, 3], [2, 2, 2], [1.0, 1.0, 1.0], 3, 3),
# 	 sparse([1, 2, 3], [3, 3, 3], [1.0, 1.0, 1.0], 3, 3)]
# end

@testset "Dualization Trivial Example 2" begin
    n = 2
    true_min = 3
    @ncpolyvar x[1:n]

    f = x[1]^2 + x[1] * x[2] + x[2] * x[1] + x[2]^2 + true_min
    r = -10
    g1 = r - x[1]
    g2 = r - x[2]
    g3 = x[1] - r
    g4 = x[2] - r

    mom_method = MomentMethod(ceil(Int, maxdegree(f) / 2), identity, x)

    pop = PolynomialOptimizationProblem(f, [g1, g2, g3, g4], x)

    model = make_sdp(mom_method, pop)

    dual_model = dualize(model, get_total_basis2var_dict(mom_method))

	set_optimizer(model, Clarabel.Optimizer)
	optimize!(model)

	set_optimizer(dual_model, Clarabel.Optimizer)
	optimize!(dual_model)

	@test is_solved_and_feasible(model)
	@test is_solved_and_feasible(dual_model)

	@test isapprox(objective_value(model), objective_value(dual_model), atol=1e-5)
end

@testset "Dualization Example 2" begin
    n = 2
    @ncpolyvar x[1:n]
    f = 2-x[1]^2+x[1]*x[2]^2*x[1]-x[2]^2
    g = 4-x[1]^2-x[2]^2
    h1 = x[1]*x[2]+x[2]*x[1]-2
	h2 = -h1
    pop = PolynomialOptimizationProblem(f, [g, h1, h2], x)

	mom_method = MomentMethod(2, identity, x)

	model = make_sdp(mom_method, pop)
	dual_model = dualize(model, get_total_basis2var_dict(mom_method))

	set_optimizer(model, Clarabel.Optimizer)
	optimize!(model)

	set_optimizer(dual_model, Clarabel.Optimizer)
	optimize!(dual_model)

	@test is_solved_and_feasible(model)
	@test is_solved_and_feasible(dual_model)
	@test isapprox(objective_value(model), -1, atol=1e-6)
	@test isapprox(objective_value(dual_model), -1, atol=1e-6)
end

@testset "Dualization Trivial Example" begin
    n = 2
	true_min = 3
    @ncpolyvar x[1:n]

    f = x[1]^2 + x[1] * x[2] + x[2] * x[1] + x[2]^2 + true_min

    mom_method = MomentMethod(ceil(Int, maxdegree(f) / 2), identity, x)

    pop = PolynomialOptimizationProblem(f, x)

    model = make_sdp(mom_method, pop)

    dual_model = dualize(model, get_total_basis2var_dict(mom_method))

    set_optimizer(model, Clarabel.Optimizer)
    optimize!(model)
	objective_value(model)

	set_optimizer(dual_model, Clarabel.Optimizer)
	optimize!(dual_model)
	objective_value(dual_model)

	dual_model[:dual_variables]

	all_variables(dual_model)
	value.(all_variables(dual_model))

	all_constraints(dual_model, include_variable_in_set_constraints=true)

    @test is_solved_and_feasible(model)
    @test is_solved_and_feasible(dual_model)
    @test isapprox(objective_value(model), true_min, atol=1e-6)
    @test isapprox(objective_value(dual_model), true_min, atol=1e-6)

end

@testset "Dualization Example 1" begin 

    n = 3
    @ncpolyvar x[1:n]

    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    mom_method = MomentMethod(ceil(Int, maxdegree(f) / 2), identity, x)

    pop = PolynomialOptimizationProblem(f, x)

    model = make_sdp(mom_method, pop)

    dual_model = dualize(model,get_total_basis2var_dict(mom_method))

    set_optimizer(model, Clarabel.Optimizer)
    optimize!(model)

	set_optimizer(dual_model, Clarabel.Optimizer)
	optimize!(dual_model)

    @test is_solved_and_feasible(model)
    @test is_solved_and_feasible(dual_model)
    @test isapprox(objective_value(model), 4.372259295498716e-10, atol=1e-8)
    @test isapprox(objective_value(dual_model), 4.372259295498716e-10, atol=1e-8)
end


