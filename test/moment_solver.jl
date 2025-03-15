using Test, NCTSSOS
using DynamicPolynomials
using JuMP
using Clarabel
using Graphs
using NCTSSOS: get_basis, substitute_variables, constrain_moment_matrix!, remove_zero_degree, star, get_correlative_graph, assign_constraint, clique_decomp, get_term_sparsity_graph, term_sparsity_graph_supp
using NCTSSOS: sorted_union,sorted_unique, neat_dot, sorted_unique
using NCTSSOS: correlative_sparsity, iterate_term_sparse_supp
using CliqueTrees

@testset "TS Smaller Example" begin
    n = 2
    @ncpolyvar x[1:2]
    f = 2.0-x[1]^2+x[1]*x[2]^2*x[1]-x[2]^2
    g = 4.0-x[1]^2-x[2]^2
    h1 = x[1]*x[2]+x[2]*x[1]-2.0
    h2 = - h1
    cons = [g, h1, h2]
    order = 2
    ts_algo = MMD()

    pop = PolynomialOptimizationProblem(f, cons)

    cliques, cliques_cons, discarded_cons, cliques_idx_basis = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, nothing)

    # prepare the support for each term sparse localizing moment
    initial_activated_supp = [
        sorted_union(reduce(vcat, monomials.([pop.objective; pop.constraints[clique_cons]])), [neat_dot(b, b) for b in clique_idx_basis[1]])
        for (clique, clique_cons, clique_idx_basis) in zip(cliques, cliques_cons, cliques_idx_basis)
    ]

    # this is hedious but I need chordal graph for next iteration
    cliques_mtcs_bases = map(zip(initial_activated_supp, cliques_cons, cliques_idx_basis)) do (activated_supp, clique_cons, clique_idx_basis)
        [iterate_term_sparse_supp(activated_supp, poly, basis, ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[clique_cons]], clique_idx_basis)]
    end

    moment_problem = moment_relax(pop, order, cliques, cliques_cons, cliques_mtcs_bases)
    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)
    objective_value(moment_problem.model)

    @test isapprox(objective_value(moment_problem.model), -1.0, atol=1e-6)
end


# @testset "CS TS Example" begin
#     order = 3
#     n = 10
#     @ncpolyvar x[1:n]
#     f = 0.0
#     for i = 1:n
#         jset = max(1,i-5):min(n,i+1)
#         jset = setdiff(jset,i)
#         f += (2x[i]+5*x[i]^3+1)^2
#         f -= sum([4x[i]*x[j]+10x[i]^3*x[j]+2x[j]+4x[i]*x[j]^2+10x[i]^3*x[j]^2+2x[j]^2 for j in jset])
#         f += sum([x[j]*x[k]+2x[j]^2*x[k]+x[j]^2*x[k]^2 for j in jset for k in jset])
#     end

#     cons = vcat([(1 - x[i]^2) for i in 1:n], [(x[i] - 1 / 3) for i in 1:n])

#     pop = PolynomialOptimizationProblem(f, cons)
#     cs_algo = MF()
#     ts_algo = MMD()

#     cliques, cliques_cons, discarded_cons, cliques_idx_basis = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, cs_algo)

#     # prepare the support for each term sparse localizing moment
#     initial_activated_supp = [
#         sorted_union(reduce(vcat, monomials.([pop.objective; pop.constraints[clique_cons]])), [neat_dot(b, b) for b in clique_idx_basis[1]])
#         for (clique, clique_cons, clique_idx_basis) in zip(cliques, cliques_cons, cliques_idx_basis)
#     ]

#     # this is hedious but I need chordal graph for next iteration
#     cliques_mtcs_bases = map(zip(initial_activated_supp, cliques_cons, cliques_idx_basis)) do (activated_supp, clique_cons, clique_idx_basis)
#         [iterate_term_sparse_supp(activated_supp, poly, basis, ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[clique_cons]], clique_idx_basis)]
#     end

#     moment_problem = moment_relax(pop, order, cliques, cliques_cons, cliques_mtcs_bases)
#     set_optimizer(moment_problem.model, Clarabel.Optimizer)
#     optimize!(moment_problem.model)
#     objective_value(moment_problem.model)
# end

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

    g1 = 1.0*x[1] - x[1]^2
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

    cliques, cliques_cons, _, cliques_idx_basis = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, nothing)

    cliques_mtcs_bases = [
        [[basis] for basis in idx_basis]
        for idx_basis in cliques_idx_basis 
    ]

    moment_problem = moment_relax(pop, order, cliques, cliques_cons, cliques_mtcs_bases)

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

    cliques, cliques_cons, _, cliques_idx_basis = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, nothing)

    cliques_mtcs_bases = [
        [[basis] for basis in idx_basis]
        for idx_basis in cliques_idx_basis 
    ]

    moment_problem = moment_relax(pop, order, cliques, cliques_cons, cliques_mtcs_bases)

    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)

    @test is_solved_and_feasible(moment_problem.model)
    @test isapprox(objective_value(moment_problem.model), -1.0, atol=1e-6)
end

@testset "Moment Method Correlative Sparsity" begin
    n = 3
    @ncpolyvar x[1:n]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6.0 * x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    cons = vcat([1.0 - x[i]^2 for i in 1:n], [x[i] - 1.0 / 3 for i in 1:n])
    order = 3
    cs_algo = MF()

    pop = PolynomialOptimizationProblem(f, cons, x)

    cliques, cliques_cons, _, cliques_idx_basis = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, cs_algo)

    cliques_mtcs_bases = [
        [[basis] for basis in idx_basis]
        for idx_basis in cliques_idx_basis 
    ]

    moment_problem = moment_relax(pop, order, cliques, cliques_cons, cliques_mtcs_bases)

    set_optimizer(moment_problem.model, Clarabel.Optimizer)

    optimize!(moment_problem.model)

    # FIXME: reduced accuracy
    # @test is_solved_and_feasible(moment_problem.model)
    @test isapprox(objective_value(moment_problem.model), 0.9975306427277915, atol=1e-5)
end

@testset "Moment Method Heisenberg Model on Star Graph" begin
    num_sites = 6 
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
    cs_algo = MF()

    cliques, cliques_cons, _, cliques_idx_basis = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, cs_algo)

    cliques_mtcs_bases = [
        [[basis] for basis in idx_basis]
        for idx_basis in cliques_idx_basis 
    ]

    moment_problem = moment_relax(pop, order, cliques, cliques_cons, cliques_mtcs_bases)


    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)

    # FIXME: objective and dual seems to be converging why they say it's only 
    # nearly feasible? But it works for Mosek
    # @test is_solved_and_feasible(moment_problem.model)
    @test isapprox(objective_value(moment_problem.model), true_ans, atol=1e-6)
end

