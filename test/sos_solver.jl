using Test, NCTSSOS
using DynamicPolynomials, Clarabel
using SparseArrays
using JuMP
using Graphs

using NCTSSOS: get_Cαj

@testset "Cαj" begin
    model = Model()
    @variable(model, x[1:4])

    cons = @constraint(model, [x[1]-x[2] x[3] x[4]+x[1]; x[1]-x[2] x[3] x[4]+x[1]; x[1]-x[2] x[3] x[4]+x[1]] in PSDCone())
    typeof(cons)

    C_α_js = get_Cαj(x, constraint_object(cons))

    @test C_α_js == [sparse([1, 2, 3, 1, 2, 3], [1, 1, 1, 3, 3, 3], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 3, 3),
        sparse([1, 2, 3], [1, 1, 1], [-1.0, -1.0, -1.0], 3, 3),
        sparse([1, 2, 3], [2, 2, 2], [1.0, 1.0, 1.0], 3, 3),
        sparse([1, 2, 3], [3, 3, 3], [1.0, 1.0, 1.0], 3, 3)]
end

@testset "Dualization Trivial Example 2" begin
    n = 2
    true_min = 3.0
    @ncpolyvar x[1:n]

    f = x[1]^2 + x[1] * x[2] + x[2] * x[1] + x[2]^2 + true_min
    r = -10.0
    g1 = r - x[1]
    g2 = r - x[2]
    g3 = x[1] - r
    g4 = x[2] - r

    pop = PolynomialOptimizationProblem(f, [g1, g2, g3, g4], x)
    order = 2

    moment_problem = moment_relax(pop, 2)

    sos_problem = sos_dualize(moment_problem)

    set_optimizer(sos_problem.model, Clarabel.Optimizer)
    optimize!(sos_problem.model)

    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)

    @test is_solved_and_feasible(moment_problem.model)
    @test is_solved_and_feasible(sos_problem.model)

    @test isapprox(objective_value(moment_problem.model), objective_value(sos_problem.model), atol=1e-5)
end

@testset "Dualization Example 2" begin
    n = 2
    @ncpolyvar x[1:n]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    h2 = -h1
    pop = PolynomialOptimizationProblem(f, [g, h1, h2], x)

    order = 2
    moment_method = moment_relax(pop, order)

    sos_method = sos_dualize(moment_method)

    set_optimizer(moment_method.model, Clarabel.Optimizer)
    optimize!(moment_method.model)

    set_optimizer(sos_method.model, Clarabel.Optimizer)
    optimize!(sos_method.model)

    @test is_solved_and_feasible(moment_method.model)
    @test is_solved_and_feasible(sos_method.model)
    @test isapprox(objective_value(moment_method.model), -1, atol=1e-6)
    @test isapprox(objective_value(sos_method.model), -1, atol=1e-6)
end

@testset "Dualization Trivial Example" begin
    n = 2
    true_min = 3.0
    @ncpolyvar x[1:n]

    f = x[1]^2 + x[1] * x[2] + x[2] * x[1] + x[2]^2 + true_min

    pop = PolynomialOptimizationProblem(f, x)
    order = 2

    moment_method = moment_relax(pop, order)

    sos_method = sos_dualize(moment_method)

    set_optimizer(moment_method.model, Clarabel.Optimizer)
    optimize!(moment_method.model)

    set_optimizer(sos_method.model, Clarabel.Optimizer)
    optimize!(sos_method.model)

    @test is_solved_and_feasible(moment_method.model)
    @test is_solved_and_feasible(sos_method.model)
    @test isapprox(objective_value(moment_method.model), true_min, atol=1e-6)
    @test isapprox(objective_value(sos_method.model), true_min, atol=1e-6)
end

@testset "Dualization Example 1" begin
    n = 3
    @ncpolyvar x[1:n]

    f =
        x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] +
        2.0 * x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 +
        9x[2]^2 * x[3] +
        9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    pop = PolynomialOptimizationProblem(f, x)
    order = 2

    moment_problem = moment_relax(pop, order)

    sos_problem = sos_dualize(moment_problem)

    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)

    set_optimizer(sos_problem.model, Clarabel.Optimizer)
    optimize!(sos_problem.model)

    @test is_solved_and_feasible(moment_problem.model)
    @test is_solved_and_feasible(sos_problem.model)
    @test isapprox(objective_value(moment_problem.model), 4.372259295498716e-10, atol=1e-6)
    @test isapprox(objective_value(sos_problem.model), 4.372259295498716e-10, atol=1e-6)
end

@testset "Dualization Heisenberg Model on Star Graph" begin
    num_sites = 8
    star = star_graph(num_sites)

    true_ans = -1.0

    vec_idx2ij = [(i, j) for i in 1:num_sites for j in (i+1):num_sites]

    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

    @ncpolyvar pij[1:length(vec_idx2ij)]

    objective = sum(1.0 * pij[[findvaridx(ee.src, ee.dst) for ee in edges(star)]])

    gs = [
        [(pij[findvaridx(i, j)]^2 - 1.0) for i in 1:num_sites for j in (i+1):num_sites]
        [-(pij[findvaridx(i, j)]^2 - 1.0) for i in 1:num_sites for j in (i+1):num_sites]
        [
            (
                pij[findvaridx(sort([i, j])...)] * pij[findvaridx(sort([j, k])...)] +
                pij[findvaridx(sort([j, k])...)] * pij[findvaridx(sort([i, j])...)] -
                pij[findvaridx(sort([i, j])...)] - pij[findvaridx(sort([j, k])...)] -
                pij[findvaridx(sort([i, k])...)] + 1.0
            ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
            (i != j && j != k && i != k)
        ]
        [
            -(
                pij[findvaridx(sort([i, j])...)] * pij[findvaridx(sort([j, k])...)] +
                pij[findvaridx(sort([j, k])...)] * pij[findvaridx(sort([i, j])...)] -
                pij[findvaridx(sort([i, j])...)] - pij[findvaridx(sort([j, k])...)] -
                pij[findvaridx(sort([i, k])...)] + 1.0
            ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
            (i != j && j != k && i != k)
        ]
    ]

    pop = PolynomialOptimizationProblem(objective, gs, pij)

    order = 1

    moment_problem = moment_relax(pop, order)

    sos_problem = sos_dualize(moment_problem)

    set_optimizer(sos_problem.model, Clarabel.Optimizer)
    optimize!(sos_problem.model)

    @test is_solved_and_feasible(sos_problem.model)
    @test isapprox(objective_value(sos_problem.model), true_ans, atol=1e-6)
end
