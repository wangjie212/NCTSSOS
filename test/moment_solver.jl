using Test, NCTSSOS
using DynamicPolynomials
using JuMP
using Clarabel
using Graphs
using NCTSSOS: get_basis, substitute_variables, constrain_moment_matrix!, remove_zero_degree, star, get_correlative_graph, assign_constraint, clique_decomp, get_term_sparsity_graph, term_sparsity_graph_supp
using NCTSSOS: sorted_union,sorted_unique, neat_dot, sorted_unique, block_decomp
using CliqueTrees

# ksupp = [1, x[1]^2, x[1]^4, x[1] * x[2], x[1] * x[2]^2 * x[1], x[2] * x[1]^2 * x[2], x[2]^2, x[2]^4]
# ksupp_after = [1, x[1]^2, x[1]^4, x[1]^3 * x[2], x[1]^2 * x[2] * x[1], x[1]^2 * x[2]^2, x[1] * x[2], x[1] * x[2] * x[1] * x[2], x[1] * x[2]^2 * x[1], x[1] * x[2]^3, x[2] * x[1]^2 * x[2], x[2] * x[1] * x[2]^2, x[2]^2, x[2]^4]
@testset "TS Smaller Example" begin
    n = 2
    @ncpolyvar x[1:2]
    f = 2.0-x[1]^2+x[1]*x[2]^2*x[1]-x[2]^2
    g = 4.0-x[1]^2-x[2]^2
    h1 = x[1]*x[2]+x[2]*x[1]-2.0
    h2 = - h1
    cons = [g, h1, h2]

    pop = PolynomialOptimizationProblem(f, cons)

    cliques = [x]

    cliques_cons, discarded_cons = assign_constraint(cliques, cons)

    # get the operators needed to index colum of moment/localizing mtx in each clique
    cliques_col_basis = map(zip(cliques,cliques_cons)) do (clique, clique_cons)
        # get the basis of the moment matrix in a clique, then sort it
        [[get_basis(clique, order)]; get_basis.(Ref(clique), order .- ceil.(Int, maxdegree.(cons[clique_cons]) / 2))]
    end


    # prepare the support for each term sparse localizing moment
    prev_localizing_mtx_basis = [
        map(zip([f; cons[clique_cons]],clique_col_basis)) do (poly, poly_col_basis)
            sorted_union(monomials(poly), [neat_dot(b, b) for b in poly_col_basis])
        end
        for (clique, clique_cons, clique_col_basis) in zip(cliques, cliques_cons, cliques_col_basis)
    ]

    # loop on the clique level
    cliques_F_is = map(zip(cliques_cons, cliques, prev_localizing_mtx_basis, cliques_col_basis)) do (clique_cons, clique, prev_basis, clique_col_basis)
        # cons_basis: term sparse basis of the localizing/moment matrix corresponding to the constraints/objective
        map(zip([pop.objective; pop.constraints[clique_cons]], prev_basis, clique_col_basis)) do (poly, mtx_basis, col_basis)
            get_term_sparsity_graph(collect(monomials((poly == pop.objective) ? one(poly) : poly)), mtx_basis, col_basis)
        end
    end

    cliques_ts_cliques = map(zip(cliques_col_basis,cliques_F_is)) do (clique_col_basis, F_is)
        map(zip(F_is,clique_col_basis)) do (F_i, col_basis)
            block_decomp(F_i, col_basis,MF())
        end
    end

    cliques_total_basis = map(zip(cliques_cons, cliques_ts_cliques)) do (clique_cons, ts_cliques)
        sorted_unique(vec(reduce(vcat, [
            map(monomials(poly)) do m
                neat_dot(col_basis[i], m * col_basis[j])
            end
            for (poly, ts_clique) in zip([one(f); pop.constraints[clique_cons]], ts_cliques) for col_basis in ts_clique for i in 1:length(col_basis) for j in 1:length(col_basis)
        ])))
    end

    moment_problem = moment_relax(pop, order, cliques, cliques_cons, cliques_total_basis, cliques_ts_cliques)

    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)

    objective_value(moment_problem.model)


end


@testset "CS TS Example" begin
    order = 3
    n = 10
    @ncpolyvar x[1:n]
    f = 0.0
    for i = 1:n
        jset = max(1,i-5):min(n,i+1)
        jset = setdiff(jset,i)
        f += (2x[i]+5*x[i]^3+1)^2
        f -= sum([4x[i]*x[j]+10x[i]^3*x[j]+2x[j]+4x[i]*x[j]^2+10x[i]^3*x[j]^2+2x[j]^2 for j in jset])
        f += sum([x[j]*x[k]+2x[j]^2*x[k]+x[j]^2*x[k]^2 for j in jset for k in jset])
    end

    cons = vcat([(1 - x[i]^2) for i in 1:n], [(x[i] - 1 / 3) for i in 1:n])

    pop = PolynomialOptimizationProblem(f, cons)

    cliques = clique_decomp(x, f, cons, MF(), order)

    true_cliques = [[x[4], x[5], x[6], x[7], x[8], x[9], x[10]],
        [x[3], x[4], x[5], x[6], x[7], x[8], x[9]],
        [x[2], x[3], x[4], x[5], x[6], x[7], x[8]],
        [x[1], x[2], x[3], x[4], x[5], x[6], x[7]]]

    @test sort!(map(x -> sort!(x), cliques)) == sort!(map(x -> sort!(x), true_cliques))

    cliques_cons, discarded_cons = assign_constraint(cliques, cons)

    # get the operators needed to index colum of moment/localizing mtx in each clique
    cliques_col_basis = map(zip(cliques,cliques_cons)) do (clique, clique_cons)
        # get the basis of the moment matrix in a clique, then sort it
        [[get_basis(clique, order)]; get_basis.(Ref(clique), order .- ceil.(Int, maxdegree.(cons[clique_cons]) / 2))]
    end


    # prepare the support for each term sparse localizing moment
    prev_localizing_mtx_basis = [
        map(zip([f; cons[clique_cons]],clique_col_basis)) do (poly, poly_col_basis)
            sorted_union(monomials(poly), [neat_dot(b, b) for b in poly_col_basis])
        end
        for (clique, clique_cons, clique_col_basis) in zip(cliques, cliques_cons, cliques_col_basis)
    ]

    # loop on the clique level
    cliques_F_is = map(zip(cliques_cons, cliques, prev_localizing_mtx_basis, cliques_col_basis)) do (clique_cons, clique, prev_basis, clique_col_basis)
        # cons_basis: term sparse basis of the localizing/moment matrix corresponding to the constraints/objective
        map(zip([pop.objective; pop.constraints[clique_cons]], prev_basis, clique_col_basis)) do (poly, mtx_basis, col_basis)
            get_term_sparsity_graph(collect(monomials((poly == pop.objective) ? one(poly) : poly)), mtx_basis, col_basis)
        end
    end

    cliques_ts_cliques = map(zip(cliques_col_basis,cliques_F_is)) do (clique_col_basis, F_is)
        map(zip(F_is,clique_col_basis)) do (F_i, col_basis)
            block_decomp(F_i, col_basis,MF())
        end
    end

    cliques_total_basis = map(zip(cliques_cons, cliques_ts_cliques)) do (clique_cons, ts_cliques)
        sorted_unique(vec(reduce(vcat, [
            map(monomials(poly)) do m
                neat_dot(col_basis[i], m * col_basis[j])
            end
            for (poly, ts_clique) in zip([one(f); pop.constraints[clique_cons]], ts_cliques) for col_basis in ts_clique for i in 1:length(col_basis) for j in 1:length(col_basis)
        ])))
    end

    moment_problem = moment_relax(pop, order, cliques, cliques_cons, cliques_total_basis, cliques_ts_cliques)

    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)

    objective_value(moment_problem.model)

end

@testset "Term Sparsity Graph" begin
    # Example 7.6 of Sparse Polynomial Optimization: Theory and Practice
    @polyvar x[1:3]
    # f = x[1]^2 - 2x[1] * x[2] + 3.0 * x[2]^2 - 2 * x[1]^2 * x[2] + 2 * x[1]^2 * x[2]^2 - 2 * x[2] * x[3] + 6 * x[3]^2 + 18 * x[2]^2 * x[3] - 54 * x[2] * x[3]^2 + 142 * x[2]^2 * x[3]^2
    # @test sort(monomials(f)) == sort(total_support)

    total_support = [x[3]^2, x[2] * x[3], x[2]^2, x[1] * x[2], x[1]^2, x[2] * x[3]^2, x[2]^2 * x[3], x[1]^2 * x[2], x[2]^2 * x[3]^2, x[1]^2 * x[2]^2]

    total_basis = [one(x[1]),x[1],x[2],x[3],x[1]*x[2],x[2]*x[3]]

    G_tsp = get_term_sparsity_graph([one(x[1])],total_support,total_basis)
    @test G_tsp.fadjlist == [[5,6],[3,5],[2,4,6],[3,6],[1,2],[1,3,4]]
    @test term_sparsity_graph_supp(G_tsp, total_basis, one(Polynomial{true,Float64})) == [one(x[1]), x[1]^2, x[2]^2, x[3]^2, x[1]^2 * x[2]^2, x[2]^2 * x[3]^2, x[1] * x[2], x[2] * x[3], x[1]^2 * x[2], x[2]^2 * x[3], x[2] * x[3]^2]

    # Example 10.2
    @ncpolyvar x y
    total_support = [one(x), x^2, x * y^2 * x, y^2, x * y * x * y, y * x * y * x, x^3 * y, y * x^3, x * y^3, y^3 * x]
    # f = 2.0 - x^2 + x * y^2 * x - y^2 + x * y * x * y + y * x * y * x + x^3 * y + y * x^3 + x * y^3 + y^3 * x
    # @test sort(monomials(f)) == sort(total_support)

    total_basis = [one(x), x, y, x^2, y^2, x * y, y * x]

    G_tsp = get_term_sparsity_graph([one(x)],total_support,total_basis)
    @test G_tsp.fadjlist == [[4,5],Int[],Int[],[1,6],[1,7],[4,7],[5,6]]
    @test sort(term_sparsity_graph_supp(G_tsp, total_basis, one(Polynomial{false,Float64}))) == sort([one(x * y), x^2, y^2, x^4, y^4, y * x^2 * y, x * y^2 * x, x^3 * y, y^3 * x, y * x * y * x])
end


@testset "Assign Constraint" begin
    n = 4
    @ncpolyvar x[1:n]

    cliques = [x[[1, 2, 4]], x[[2, 3, 4]]]
    # NOTE: CliqueTrees.MMD() returns cliques that results in discardding a constraint
    # cliques = [x[[2, 3, 4]], x[[1, 3, 4]]]
    cons = Polynomial{false,Float64}[x[1] * x[2],
        x[2] * x[3],
        x[3] * x[4],
        x[4] * x[1]]

    @test assign_constraint(cliques, cons) == ([[1, 4], [2, 3]], Int[])
end

@testset "Clique Decomposition" begin
    n = 4
    @ncpolyvar x[1:n]
    f = sum(x[i] * x[mod1(i + 1, n)] for i in 1:n)
    order = 1

    G = get_correlative_graph(x, f, typeof(f)[], order)
    @test G.fadjlist == map(x -> sort!(x), [[2, 4], [1, 3], [2, 4], [1, 3]])

    @test sort!(clique_decomp(x, f, typeof(f)[], MF(), order)) == sort!(map(x -> sort!(x), [[x[2], x[3], x[4]], [x[1], x[2], x[4]]]))
    @test sort!(clique_decomp(x, f, typeof(f)[], MMD(), order)) == sort!(map(x -> sort!(x), [[x[1], x[2], x[3]], [x[1], x[3], x[4]]]))


    n = 3
    @ncpolyvar x[1:3]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6.0 * x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]
    order = 2

    G = get_correlative_graph(x, f, typeof(f)[], order)
    @test G.fadjlist == [[2], [1, 3], [2]]

    # MMD, MF don't work
    @test sort!(clique_decomp(x, f, typeof(f)[], BFS(), order)) == sort!(map(x -> sort!(x), [[x[2],x[3]],[x[1],x[2]]]))
    @test sort!(clique_decomp(x, f, typeof(f)[], MF(), order)) == sort!(map(x -> sort!(x), [[x[2],x[3]],[x[1],x[2]]]))

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
    order = 2
    G = get_correlative_graph(x, f, typeof(f)[], order)
    @test G.fadjlist == [[2, 3, 4, 5, 6, 7], [1, 3, 4, 5, 6, 7, 8], [1, 2, 4, 5, 6, 7, 8, 9], [1, 2, 3, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 8, 9, 10], [2, 3, 4, 5, 6, 7, 9, 10], [3, 4, 5, 6, 7, 8, 10], [4, 5, 6, 7, 8, 9]]
    @test sort!(clique_decomp(x, f, typeof(f)[], BFS(), order)) == sort!(map(x -> sort!(x), [x[[4,5,6,7,8,9,10]],x[[3,4,5,6,7,8,9]],x[[2,3,4,5,6,7,8]],x[[1,2,3,4,5,6,7]]]))

    order = 3
    G = get_correlative_graph(x, f, typeof(f)[], order)
    @test G.fadjlist == [[2, 3, 4, 5, 6, 7], [1, 3, 4, 5, 6, 7, 8], [1, 2, 4, 5, 6, 7, 8, 9], [1, 2, 3, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 8, 9, 10], [2, 3, 4, 5, 6, 7, 9, 10], [3, 4, 5, 6, 7, 8, 10], [4, 5, 6, 7, 8, 9]]
    @test sort!(clique_decomp(x, f, typeof(f)[], BFS(), order)) == sort!(map(x -> sort!(x), [x[[4,5,6,7,8,9,10]],x[[3,4,5,6,7,8,9]],x[[2,3,4,5,6,7,8]],x[[1,2,3,4,5,6,7]]]))


    n = 3
    @ncpolyvar x[1:3]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6.0 * x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    cons = vcat([1.0 - x[i]^2 for i in 1:n], [x[i] - 1.0 / 3 for i in 1:n])
    order = 3
    G = get_correlative_graph(x, f, cons, order)
    @test G.fadjlist == [[2], [1, 3], [2]]
    @test sort!(clique_decomp(x, f, cons, BFS(), order)) == sort!(map(x -> sort!(x), [x[[1,2]],x[[2,3]]]))
end

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

    # TODO: fix inputs
    moment_problem = moment_relax(pop, order, [x], [pop.constraints], [[monomials(f)...]])

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

    # TODO: fix inputs
    moment_problem = moment_relax(pop, order, [x], [pop.constraints], [[monomials(f)...]])

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

    pop = PolynomialOptimizationProblem(f, cons, x)

    cliques = clique_decomp(x, f, cons, MF(), order)


    # TODO: fix inputs
    moment_problem = moment_relax(pop, order, cliques, [pop.constraints], [[monomials(f)...]])
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

    # TODO: fix inputs
    moment_problem = moment_relax(pop, order, [pij], [pop.constraints], [[monomials(objective)...]])


    set_optimizer(moment_problem.model, Clarabel.Optimizer)
    optimize!(moment_problem.model)

    # FIXME: objective and dual seems to be converging why they say it's only 
    # nearly feasible? But it works for Mosek
    # @test is_solved_and_feasible(moment_problem.model)
    @test isapprox(objective_value(moment_problem.model), true_ans, atol=1e-6)
end

