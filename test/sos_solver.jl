using Test, NCTSSOS, DynamicPolynomials, Clarabel
using JuMP

@test "Dualization" begin

    n = 3
    @ncpolyvar x[1:n]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

	g = sum(x.^2) - 10

    mom_method = MomentMethod(ceil(Int, maxdegree(f) / 2), identity, x)

    pop = PolynomialOptimizationProblem(f, [g,] , x)

    model = make_sdp(mom_method, pop)

    dual_model = dualize(mom_method, model)

    set_optimizer(model, Clarabel.Optimizer)

    optimize!(model)

	for (basis, jump_var) in get_total_basis2var_dict(mom_method)
		@show basis
		@show jump_var
	end

	objective_function(model)



	cmm = constraint_object(cons_mom_mtx)


	cmm.func[1]

	mom_mtx_cons_shape = JuMP.shape(constraint_object(cons_mom_mtx))

    cons_loc_mtx = constraint_by_name(model, "localizing_matrix_$(hash(g))")

	clm = constraint_object(cons_loc_mtx)

	aff_expr = clm.func[1]

	aff_expr.terms

	typeof(clm.func[1])

	loc_mtx_cons_shape = JuMP.shape(constraint_object(cons_loc_mtx))


	loc_mtx_cons_shape.side_dimension


	set_optimizer(dual_model, Clarabel.Optimizer)
	optimize!(dual_model)

    @test is_solved_and_feasible(model)
    @test is_solved_and_feasible(dual_model)
    @test isapprox(objective_value(model), 4.372259295498716e-10, atol=1e-8)
    @test isapprox(objective_value(dual_model), 4.372259295498716e-10, atol=1e-8)

end


