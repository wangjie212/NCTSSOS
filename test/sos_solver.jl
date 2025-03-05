using Test, NCTSSOS, DynamicPolynomials, Clarabel
using JuMP

@test "Dualization" begin

    n = 3
    @ncpolyvar x[1:n]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

	g = sum(x.^2) - 10

    mom_method = MomentMethod(ceil(Int, maxdegree(f) / 2), identity, x)

    pop = PolynomialOptimizationProblem(f, x)

    model = make_sdp(mom_method, pop)

    dual_model = dualize(model)

    set_optimizer(model, Clarabel.Optimizer)
    optimize!(model)

	set_optimizer(dual_model, Clarabel.Optimizer)
	optimize!(dual_model)

	objective_value(dual_model)

    constraint_object(all_constraints(dual_model, include_variable_in_set_constraints=true)[2]).func[1]


    @test is_solved_and_feasible(model)
    @test is_solved_and_feasible(dual_model)
    @test isapprox(objective_value(model), 4.372259295498716e-10, atol=1e-8)
    @test isapprox(objective_value(dual_model), 4.372259295498716e-10, atol=1e-8)

end


