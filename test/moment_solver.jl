using Test, NCTSSOS
using DynamicPolynomials
using JuMP
using Clarabel
using Graphs
using NCTSSOS: get_basis, substitute_variables, constrain_moment_matrix!, remove_zero_degree, star

@testset "Replace DynamicPolynomials variables with JuMP variables" begin
    @ncpolyvar x y z
    poly = 1.0 * x^2 - 2.0 * x * y - 1

    model = GenericModel{Float64}()
    @variable(model, jm[1:13])

    monomap = Dict(get_basis([x,y,z],2) .=> jm)

    @test substitute_variables(poly, monomap) == 1.0*jm[5] - 2.0 * jm[6] - jm[1]
end

@testset "Constrain Moment Matrix" begin
    # Example 2.4 in Sparse Polynomial Optimization: Theory and Practice 
    model = GenericModel{Float64}()
    order = 2
    @variable(model, y[1:15])
    @polyvar x[1:2]

    moment_mtx_idcs = get_basis(x, order)
    total_basis = sort(
        unique([
            remove_zero_degree(row_idx * col_idx) for row_idx in star.(moment_mtx_idcs) for
            col_idx in moment_mtx_idcs
        ]),
    )
    monomap = Dict(total_basis .=> y)

    g1 = x[1] - x[1]^2
    local_basis = [one(x[1]),x[1],x[2]]
    localizing_matrix_constraint = constrain_moment_matrix!(model, g1, local_basis, monomap)


    myexponents = [([1, 0], [2, 0]) ([2, 0], [3, 0]) ([1, 1], [2, 1]); ([2, 0], [3, 0]) ([3, 0], [4, 0]) ([2, 1], [3, 1]); ([1, 1], [2, 1]) ([2, 1], [3, 1]) ([1, 2], [2, 2])]

    true_localizing_matrix = [mapreduce(
        a -> (monomap[prod([iszero(b[2]) ? one(x[1]) : x[b[1]]^b[2] for b in enumerate(a)])]), -, myexponents[i]) for i in eachindex(myexponents)]

    @test true_localizing_matrix == ((x->getfield(x,:func)) ∘ JuMP.constraint_object)(localizing_matrix_constraint)
end

@testset "Moment Method Example 1" begin
    order = 2
    n = 3
    @ncpolyvar x[1:n]
    f =
        x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] +
        2.0x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 +
        9x[2]^2 * x[3] +
        9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    pop = PolynomialOptimizationProblem(f, x)

    moment_problem = moment_relax(pop, order)

    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)

    # NOTE: differs from original test case value since that one is a relaxed in terms of sparsity
    # This value here is obtained by running the master branch with no sparsity relaxation
    @test is_solved_and_feasible(moment_problem.model)
    @test isapprox(objective_value(moment_problem.model), 4.372259295498716e-10, atol=1e-6)
end

@testset "Moment Method Example 2" begin
    order = 2
    n = 2
    @ncpolyvar x[1:2]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    h2 = -h1
    pop = PolynomialOptimizationProblem(f, [g, h1, h2], x)

    moment_problem = moment_relax(pop, order)

    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)
    @test is_solved_and_feasible(moment_problem.model)
    @test isapprox(objective_value(moment_problem.model), -1.0, atol=1e-6)
end

@testset "Moment Method Example 3" begin
    # NOTE: this is not doable in non-sparse case
    # order = 2
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

    # moment_problem = moment_relax(pop,order )
    # set_optimizer(moment_problem.model, Clarabel.Optimizer)
    # optimize!(moment_problem.model)

    # @test isapprox(objective_value(moment_problem.model), 0.0, atol=1e-4)
end


@testset "Moment Method Heisenberg Model on Star Graph" begin
    num_sites = 8 
    g = star_graph(num_sites)

    true_ans = -1.0

    vec_idx2ij = [(i, j) for i in 1:num_sites for j in (i + 1):num_sites]

    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

    @ncpolyvar pij[1:length(vec_idx2ij)]

    objective = sum(1.0*pij[[findvaridx(ee.src, ee.dst) for ee in edges(g)]])

    gs = [
        [(pij[findvaridx(i, j)]^2 - 1.0) for i in 1:num_sites for j in (i + 1):num_sites]
        [-(pij[findvaridx(i, j)]^2 - 1.0) for i in 1:num_sites for j in (i + 1):num_sites]
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

    moment_problem = moment_relax(pop, order )


    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)

    # FIXME: objective and dual seems to be converging why they say it's only 
    # nearly feasible? But it works for Mosek
    @test is_solved_and_feasible(moment_problem.model)
    @test isapprox(objective_value(moment_problem.model), true_ans, atol=1e-6)
end
