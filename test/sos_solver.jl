using Test, NCTSSOS
using DynamicPolynomials, Clarabel
using SparseArrays
using JuMP
using Graphs
using CliqueTrees
using NCTSSOS: get_Cαj, clique_decomp, correlative_sparsity, sorted_union, neat_dot, iterate_term_sparse_supp, symmetric_canonicalize, TermSparsity, moment_relax, sos_dualize

# TODO: add Broyden banded for larger test case

# NOTE: sos_dualize has performance issue have verified locally it's correct
@testset "CS TS Example" begin
    order = 3
    n = 10
    @ncpolyvar x[1:n]
    f = 0.0
    for i = 1:n
        jset = max(1, i - 5):min(n, i + 1)
        jset = setdiff(jset, i)
        f += (2x[i] + 5 * x[i]^3 + 1)^2
        f -= sum([4x[i] * x[j] + 10x[i]^3 * x[j] + 2x[j] + 4x[i] * x[j]^2 + 10x[i]^3 * x[j]^2 + 2x[j]^2 for j in jset])
        f += sum([x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset])
    end

    cons = vcat([(1 - x[i]^2) for i in 1:n], [(x[i] - 1 / 3) for i in 1:n])

    pop = PolyOpt(f, cons)
    cs_algo = MF()
    ts_algo = MMD()

    corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, cs_algo)

    cliques_objective = [reduce(+, [issubset(effective_variables(mono), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

    # prepare the support for each term sparse localizing moment
    initial_activated_supp = [
        sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, pop.constraints[cons_idx]), [neat_dot(b, b) for b in idcs_bases[1]])
        for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)
    ]

    cliques_term_sparsities = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
        [iterate_term_sparse_supp(activated_supp, poly, basis, ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
    end

    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)
    sos_problem = sos_dualize(moment_problem)

    set_optimizer(sos_problem.model, Clarabel.Optimizer)
    optimize!(sos_problem.model)
    @test isapprox(objective_value(sos_problem.model), 3.011288, atol=1e-4)
end

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

    pop = PolyOpt(f, [g1, g2, g3, g4])
    order = 2

    corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, NoElimination())

    cliques_term_sparsities = [
        [TermSparsity(Vector{Monomial{false}}(), [basis]) for basis in idx_basis]
        for idx_basis in corr_sparsity.cliques_idcs_bases
    ]

    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)

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
    pop = PolyOpt(f, [g, h1, h2])

    order = 2

    corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, NoElimination())

    @testset "Dense" begin
        cliques_term_sparsities = [
            [TermSparsity(Vector{Monomial{false}}(), [basis]) for basis in idx_basis]
            for idx_basis in corr_sparsity.cliques_idcs_bases
        ]

        moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)
        sos_problem = sos_dualize(moment_problem)

        set_optimizer(sos_problem.model, Clarabel.Optimizer)
        optimize!(sos_problem.model)

        @test is_solved_and_feasible(sos_problem.model)
        @test isapprox(objective_value(sos_problem.model), -1, atol=1e-6)
    end

    @testset "Term Sparse" begin
        ts_algo = MMD()

        cliques_objective = [reduce(+, [issubset(effective_variables(mono), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

        # prepare the support for each term sparse localizing moment
        initial_activated_supp = [
            # why does order matter here?
            sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, pop.constraints[cons_idx]), [neat_dot(b, b) for b in idcs_bases[1]])
            for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)
        ]

        cliques_term_sparsities = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
            [iterate_term_sparse_supp(activated_supp, poly, basis, ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
        end

        moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)
        sos_problem = sos_dualize(moment_problem)
        set_optimizer(sos_problem.model, Clarabel.Optimizer)
        optimize!(sos_problem.model)
        objective_value(sos_problem.model)

        @test isapprox(objective_value(sos_problem.model), -1.0, atol=1e-6)
    end
end

@testset "Dualization Trivial Example" begin
    n = 2
    true_min = 3.0
    @ncpolyvar x[1:n]

    f = x[1]^2 + x[1] * x[2] + x[2] * x[1] + x[2]^2 + true_min

    pop = PolyOpt(f, typeof(f)[])
    order = 2

    corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, NoElimination())

    cliques_term_sparsities = [
        [TermSparsity(Vector{Monomial{false}}(), [basis]) for basis in idx_basis]
        for idx_basis in corr_sparsity.cliques_idcs_bases
    ]

    moment_method = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)

    sos_method = sos_dualize(moment_method)

    set_optimizer(sos_method.model, Clarabel.Optimizer)
    optimize!(sos_method.model)

    @test is_solved_and_feasible(sos_method.model)
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

    pop = PolyOpt(f, typeof(f)[])
    order = 2

    corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, NoElimination())

    cliques_term_sparsities = [
        [TermSparsity(Vector{Monomial{false}}(), [basis]) for basis in idx_basis]
        for idx_basis in corr_sparsity.cliques_idcs_bases
    ]

    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)

    sos_problem = sos_dualize(moment_problem)

    set_optimizer(sos_problem.model, Clarabel.Optimizer)
    optimize!(sos_problem.model)

    @test is_solved_and_feasible(sos_problem.model)
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

    pop = PolyOpt(objective, gs)

    order = 1

    corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, NoElimination())

    cliques_term_sparsities = [
        [TermSparsity(Vector{Monomial{false}}(), [basis]) for basis in idx_basis]
        for idx_basis in corr_sparsity.cliques_idcs_bases
    ]

    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)

    sos_problem = sos_dualize(moment_problem)

    set_optimizer(sos_problem.model, Clarabel.Optimizer)
    optimize!(sos_problem.model)

    @test is_solved_and_feasible(sos_problem.model)
    @test isapprox(objective_value(sos_problem.model), true_ans, atol=1e-6)
end

@testset "SOS Method Correlative Sparsity" begin
    n = 3
    @ncpolyvar x[1:n]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6.0 * x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    cons = vcat([1.0 - x[i]^2 for i in 1:n], [x[i] - 1.0 / 3 for i in 1:n])
    order = 3

    cs_algo = MF()

    pop = PolyOpt(f, cons)

    corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, cs_algo)

    cliques_term_sparsities = [
        [TermSparsity(Vector{Monomial{false}}(), [basis]) for basis in idx_basis]
        for idx_basis in corr_sparsity.cliques_idcs_bases
    ]

    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)

    sos_problem = sos_dualize(moment_problem)

    set_optimizer(sos_problem.model, Clarabel.Optimizer)

    optimize!(sos_problem.model)

    # FIXME: reduced accuracy
    @test is_solved_and_feasible(sos_problem.model)
    @test isapprox(objective_value(sos_problem.model), 0.9975306427277915, atol=1e-5)
end

@testset "Benchmark" begin
    n = 3
    @ncpolyvar x[1:n]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6.0 * x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    cons = vcat([1.0 - x[i]^2 for i in 1:n], [x[i] - 1.0 / 3 for i in 1:n])
    order = 3

    cs_algo = MF()

    pop = PolyOpt(f, cons)

    corr_sparsity_s = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, cs_algo)

    cliques_term_sparsities_s = [
        [TermSparsity(Vector{Monomial{false}}(), [basis]) for basis in idx_basis]
        for idx_basis in corr_sparsity_s.cliques_idcs_bases
    ]

    corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, order, NoElimination())

    cliques_term_sparsities = [
        [TermSparsity(Vector{Monomial{false}}(), [basis]) for basis in idx_basis]
        for idx_basis in corr_sparsity.cliques_idcs_bases
    ]

    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)

    moment_problem_s = moment_relax(pop, corr_sparsity_s.cliques_cons, corr_sparsity_s.global_cons, cliques_term_sparsities_s)

    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    set_optimizer(moment_problem_s.model, Clarabel.Optimizer)

    optimize!(moment_problem.model)
    optimize!(moment_problem_s.model)

    @test isapprox(objective_value(moment_problem.model), objective_value(moment_problem_s.model), atol=1e-5)
    @test solve_time(moment_problem.model) > solve_time(moment_problem_s.model)

    sos_problem = sos_dualize(moment_problem)
    sos_problem_s = sos_dualize(moment_problem_s)

    set_optimizer(sos_problem.model, Clarabel.Optimizer)
    set_optimizer(sos_problem_s.model, Clarabel.Optimizer)

    optimize!(sos_problem.model)
    optimize!(sos_problem_s.model)

    @test isapprox(objective_value(sos_problem.model), objective_value(sos_problem_s.model), atol=1e-5)
    @test solve_time(sos_problem.model) > solve_time(sos_problem_s.model)
end
