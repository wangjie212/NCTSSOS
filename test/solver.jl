using Test, NCTSSOS
using MosekTools 
using DynamicPolynomials
using JuMP
using COSMO


@testset "Moment Method" begin
    n = 3
    @ncpolyvar x[1:n]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    mom_method = MomentMethod(ceil(Int, maxdegree(f) / 2), identity)

    pop = PolynomialOptimizationProblem(f, x)

    model = make_sdp(mom_method, pop)

	set_optimizer(model, COSMO.Optimizer)
	optimize!(model)

    value.(first(all_constraints(model, include_variable_in_set_constraints=true)))
	display(value.(all_constraints(model, include_variable_in_set_constraints=true)[2]))
	is_solved_and_feasible(model)

	all_variables(model)

	objective_value(model)
end

@testset "Moment Method Example 2" begin
    n = 2
    @ncpolyvar x[1:2]
    f = 2-x[1]^2+x[1]*x[2]^2*x[1]-x[2]^2
    g = 4-x[1]^2-x[2]^2
    h1 = x[1]*x[2]+x[2]*x[1]-2
	h2 = -h1
    pop = PolynomialOptimizationProblem(f, [g, h1, h2], x)

	mom_method = MomentMethod(2, identity)

	model = make_sdp(mom_method, pop)

	set_optimizer(model, Mosek.Optimizer)
	optimize!(model)
	is_solved_and_feasible(model)

	objective_value(model)

end


@testset "Init Moment Vector" begin
	@ncpolyvar x[1:3]

	model = Model()

	total_basis = get_basis(x,1)

	jump_vars = init_moment_vector!(model, total_basis)

	y = all_variables(model)

	@test jump_vars == Dict([one(x[1])=> y[1],x[1]=> y[2],x[2]=> y[3],x[3]=> y[4]])

	@ncpolyvar x[1:2]

	model = Model()

	total_basis = get_basis(x,2)

    jump_vars = init_moment_vector!(model, [row_idx * col_idx for row_idx in star.(total_basis) for col_idx in total_basis])

    for var in keys(jump_vars)
		display(var)
	end
end
