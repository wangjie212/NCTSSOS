using Test, NCTSSOS
using DynamicPolynomials
using JuMP
using Clarabel
using Graphs

@testset "Moment Method Example 1" begin
    order = 2
    n = 3
    @ncpolyvar x[1:n]
    f =
        x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] +
        2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 +
        9x[2]^2 * x[3] +
        9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    pop = PolynomialOptimizationProblem(f, x)

    model = moment_relax(pop, order)

    set_optimizer(model, Clarabel.Optimizer)
    optimize!(model)

    display(value.(all_constraints(model; include_variable_in_set_constraints=true)[2]))

    is_solved_and_feasible(model)

    # NOTE: differs from original test case value since that one is a relaxed in terms of sparsity
    # This value here is obtained by running the master branch with no sparsity relaxation
    objective_value(model)
    @test isapprox(objective_value(model), 4.372259295498716e-10, atol=1e-8)
end

@testset "Moment Method Example 2" begin
    n = 2
    @ncpolyvar x[1:2]
    f = 2 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2
    h2 = -h1
    pop = PolynomialOptimizationProblem(f, [g, h1, h2], x)

    mom_method = MomentMethod(2, identity, x)

    model = make_sdp(mom_method, pop)

    set_optimizer(model, Clarabel.Optimizer)
    optimize!(model)
    is_solved_and_feasible(model)

    @test isapprox(objective_value(model), -1.0, atol=1e-6)
end

@testset "Moment Method Example 3" begin
    # NOTE: this is not doable in non-sparse case
    # n = 10
    # @ncpolyvar x[1:n]
    # f = 0.0
    # for i = 1:n
    #     jset = max(1,i-5):min(n,i+1)
    #     jset = setdiff(jset,i)
    #     f += (2x[i]+5*x[i]^3+1)^2
    #     f -= sum([4x[i]*x[j]+10x[i]^3*x[j]+2x[j]+4x[i]*x[j]^2+10x[i]^3*x[j]^2+2x[j]^2 for j in jset])
    #     f += sum([x[j]*x[k]+2x[j]^2*x[k]+x[j]^2*x[k]^2 for j in jset for k in jset])
    # end

    # pop = PolynomialOptimizationProblem(f, x)

    # mom_method = MomentMethod(3, identity, x)

    # model = make_sdp(mom_method, pop)

    # set_optimizer(model, Clarabel.Optimizer)
    # optimize!(model)

    # @test isapprox(objective_value(model), 0.0, atol=1e-4)
end

@testset "Init Moment Vector" begin
    @ncpolyvar x[1:3]

    model = Model()

    total_basis = get_basis(x, 1)

    jump_vars = init_moment_vector!(model, total_basis)

    y = all_variables(model)

    @test jump_vars == Dict([one(x[1]) => y[1], x[1] => y[2], x[2] => y[3], x[3] => y[4]])

    @ncpolyvar x[1:2]

    model = Model()

    total_basis = get_basis(x, 2)

    jump_vars = init_moment_vector!(
        model,
        [row_idx * col_idx for row_idx in star.(total_basis) for col_idx in total_basis],
    )

    for var in keys(jump_vars)
        display(var)
    end
end

@testset "Moment Method Heisenberg Model on Star Graph" begin
    num_sites = 4
    g = star_graph(num_sites)

    true_ans = -1.0

    vec_idx2ij = [(i, j) for i in 1:num_sites for j in (i + 1):num_sites]

    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

    @ncpolyvar pij[1:length(vec_idx2ij)]

    objective = sum(pij[[findvaridx(ee.src, ee.dst) for ee in edges(g)]])

    gs = [
        [(pij[findvaridx(i, j)]^2 - 1) for i in 1:num_sites for j in (i + 1):num_sites]
        [-(pij[findvaridx(i, j)]^2 - 1) for i in 1:num_sites for j in (i + 1):num_sites]
        [
            (
                pij[findvaridx(sort([i, j])...)] * pij[findvaridx(sort([j, k])...)] +
                pij[findvaridx(sort([j, k])...)] * pij[findvaridx(sort([i, j])...)] -
                pij[findvaridx(sort([i, j])...)] - pij[findvaridx(sort([j, k])...)] -
                pij[findvaridx(sort([i, k])...)] + 1
            ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
            (i != j && j != k && i != k)
        ]
        [
            -(
                pij[findvaridx(sort([i, j])...)] * pij[findvaridx(sort([j, k])...)] +
                pij[findvaridx(sort([j, k])...)] * pij[findvaridx(sort([i, j])...)] -
                pij[findvaridx(sort([i, j])...)] - pij[findvaridx(sort([j, k])...)] -
                pij[findvaridx(sort([i, k])...)] + 1
            ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
            (i != j && j != k && i != k)
        ]
    ]

    pop = PolynomialOptimizationProblem(objective, gs, pij)

    method = MomentMethod(1, identity, pij)

    model = make_sdp(method, pop)

    set_optimizer(model, Clarabel.Optimizer)
    optimize!(model)

    # FIXME: objective and dual seems to be converging why they say it's only 
    # nearly feasible? But it works for Mosek
    @test is_solved_and_feasible(model)
    @test isapprox(objective_value(model), true_ans, atol=1e-6)
end
