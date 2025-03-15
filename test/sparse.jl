using Test, NCTSSOS
using Graphs
using CliqueTrees
using NCTSSOS: assign_constraint, get_correlative_graph, clique_decomp, get_term_sparsity_graph, term_sparsity_graph_supp

@testset "Term Sparsity Graph" begin
    # Example 7.6 of Sparse Polynomial Optimization: Theory and Practice
    @polyvar x[1:3]
    # f = x[1]^2 - 2x[1] * x[2] + 3.0 * x[2]^2 - 2 * x[1]^2 * x[2] + 2 * x[1]^2 * x[2]^2 - 2 * x[2] * x[3] + 6 * x[3]^2 + 18 * x[2]^2 * x[3] - 54 * x[2] * x[3]^2 + 142 * x[2]^2 * x[3]^2
    # @test sort(monomials(f)) == sort(total_support)

    activated_support = [x[3]^2, x[2] * x[3], x[2]^2, x[1] * x[2], x[1]^2, x[2] * x[3]^2, x[2]^2 * x[3], x[1]^2 * x[2], x[2]^2 * x[3]^2, x[1]^2 * x[2]^2]

    mtx_basis = [one(x[1]),x[1],x[2],x[3],x[1]*x[2],x[2]*x[3]]

    G_tsp = get_term_sparsity_graph([one(x[1])],activated_support,mtx_basis)
    @test G_tsp.fadjlist == [[5,6],[3,5],[2,4,6],[3,6],[1,2],[1,3,4]]
    @test term_sparsity_graph_supp(G_tsp, mtx_basis, one(Polynomial{true,Float64})) == [one(x[1]), x[1]^2, x[2]^2, x[3]^2, x[1]^2 * x[2]^2, x[2]^2 * x[3]^2, x[1] * x[2], x[2] * x[3], x[1]^2 * x[2], x[2]^2 * x[3], x[2] * x[3]^2]

    # Example 10.2
    @ncpolyvar x y
    activated_support = [one(x), x^2, x * y^2 * x, y^2, x * y * x * y, y * x * y * x, x^3 * y, y * x^3, x * y^3, y^3 * x]
    # f = 2.0 - x^2 + x * y^2 * x - y^2 + x * y * x * y + y * x * y * x + x^3 * y + y * x^3 + x * y^3 + y^3 * x
    # @test sort(monomials(f)) == sort(total_support)

    mtx_basis = [one(x), x, y, x^2, y^2, x * y, y * x]

    G_tsp = get_term_sparsity_graph([one(x)],activated_support,mtx_basis)
    @test G_tsp.fadjlist == [[4,5],Int[],Int[],[1,6],[1,7],[4,7],[5,6]]
    @test sort(term_sparsity_graph_supp(G_tsp, mtx_basis, one(Polynomial{false,Float64}))) == sort([one(x * y), x^2, y^2, x^4, y^4, y * x^2 * y, x * y^2 * x, x^3 * y, y^3 * x, y * x * y * x])
end


@testset "Get Term Sparsity Graph" begin
    n = 4
    @ncpolyvar x[1:n]
    f = sum(x[i] * x[mod1(i + 1, n)] for i in 1:n)
    order = 1
    
end

@testset "Assign Constraint" begin
    n = 4
    @ncpolyvar x[1:n]

    cliques = [x[[1, 2, 4]], x[[2, 3, 4]]]
    cons = Polynomial{false,Float64}[x[1] * x[2],
        x[2] * x[3],
        x[3] * x[4],
        x[4] * x[1]]

    @test assign_constraint(cliques, cons) == ([[1, 4], [2, 3]], Int[])

    n = 2
    @ncpolyvar x[1:2]
    g = 4.0-x[1]^2-x[2]^2
    h1 = x[1]*x[2]+x[2]*x[1]-2.0
    h2 = - h1
    cons = [g, h1, h2]

    cliques = [x]

    @test assign_constraint(cliques, cons) == ([[1,2,3]], Int[])
end

@testset "Clique Decomposition" begin
    G = SimpleGraph(Edge.([(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (4, 10), (5, 6), (5, 7), (5, 8), (5, 9), (5, 10), (6, 7), (6, 8), (6, 9), (6, 10), (7, 8), (7, 9), (7, 10), (8, 9), (8, 10), (9, 10)]))

    true_cliques = [[4, 5, 6, 7, 8, 9, 10],
        [3, 4, 5, 6, 7, 8, 9],
        [2, 3, 4, 5, 6, 7, 8],
        [1, 2, 3, 4, 5, 6, 7]]
    
    @test sort!(map(x -> sort!(x),clique_decomp(G,MF()))) == sort!(true_cliques)

    # TODO: add more test
    G = SimpleGraph(Edge.([(1, 2), (1, 4), (2, 3), (3, 4)]))

    @test sort!(map(x -> sort!(x),clique_decomp(G, MF()))) == sort!(map(x -> sort!(x), [[2, 3, 4], [1, 2, 4]]))
    @test sort!(map(x -> sort!(x),clique_decomp(G, MMD()))) == sort!(map(x -> sort!(x), [[1, 2, 3], [1, 3, 4]]))

    G = SimpleGraph(Edge.([(1, 2),(2, 3)]))

    @test sort!(map(x -> sort!(x),clique_decomp(G, BFS()))) == sort!(map(x -> sort!(x), [[2, 3], [1, 2]]))
    @test sort!(map(x -> sort!(x),clique_decomp(G, MF()))) == sort!(map(x -> sort!(x), [[2, 3], [1, 2]]))
    @test sort!(map(x -> sort!(x),clique_decomp(G, MMD()))) == sort!(map(x -> sort!(x), [[2, 3], [1, 2]]))

    G = SimpleGraph(Edge.([(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (4, 10), (5, 6), (5, 7), (5, 8), (5, 9), (5, 10), (6, 7), (6, 8), (6, 9), (6, 10), (7, 8), (7, 9), (7, 10), (8, 9), (8, 10), (9, 10)]))

    @test sort!(map(x -> sort!(x),clique_decomp(G, BFS()))) == sort!(map(x -> sort!(x), [[4, 5, 6, 7, 8, 9, 10], [3, 4, 5, 6, 7, 8, 9], [2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7]]))
    @test sort!(map(x -> sort!(x),clique_decomp(G, MF()))) == sort!(map(x -> sort!(x), [[4, 5, 6, 7, 8, 9, 10], [3, 4, 5, 6, 7, 8, 9], [2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7]]))

    G = SimpleGraph(Edge.([(1, 2),(1, 3),(1, 4),(1, 5),(1, 6),(1, 7),(2, 3),(2, 4),(2, 5),(2, 6),(2, 7),(2, 8),(3, 4),(3, 5),(3, 6),(3, 7),(3, 8),(3, 9),(4, 5),(4, 6),(4, 7),(4, 8),(4, 9),(4, 10),(5, 6),(5, 7),(5, 8),(5, 9),(5, 10),(6, 7),(6, 8),(6, 9),(6, 10),(7, 8),(7, 9),(7, 10),(8, 9),(8, 10),(9, 10)]))

    @test sort!(map(x -> sort!(x), clique_decomp(G, BFS()))) == sort!(map(x -> sort!(x), [[4, 5, 6, 7, 8, 9, 10], [3, 4, 5, 6, 7, 8, 9], [2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7]]))

    G = SimpleGraph(Edge.([(1, 2), (2, 3)]))

    @test sort!(map(x -> sort!(x), clique_decomp(G, BFS()))) == sort!(map(x -> sort!(x), [[2, 3], [1, 2]]))
end

@testset "Get Correlative Graph" begin
    n = 4
    @ncpolyvar x[1:n]
    f = sum(x[i] * x[mod1(i + 1, n)] for i in 1:n)
    order = 1

    G = get_correlative_graph(x, f, typeof(f)[], order)
    @test_throws AssertionError get_correlative_graph(sort(x), f, typeof(f)[], order)
    @test G.fadjlist == map(x -> sort!(x), [[2, 4], [1, 3], [2, 4], [1, 3]])

    n = 3
    @ncpolyvar x[1:3]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6.0 * x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]
    order = 2

    G = get_correlative_graph(x, f, typeof(f)[], order)
    @test G.fadjlist == [[2], [1, 3], [2]]

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

    order = 3
    G = get_correlative_graph(x, f, typeof(f)[], order)
    @test G.fadjlist == [[2, 3, 4, 5, 6, 7], [1, 3, 4, 5, 6, 7, 8], [1, 2, 4, 5, 6, 7, 8, 9], [1, 2, 3, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 8, 9, 10], [2, 3, 4, 5, 6, 7, 9, 10], [3, 4, 5, 6, 7, 8, 10], [4, 5, 6, 7, 8, 9]]


    n = 3
    @ncpolyvar x[1:3]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6.0 * x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    cons = vcat([1.0 - x[i]^2 for i in 1:n], [x[i] - 1.0 / 3 for i in 1:n])
    order = 3
    G = get_correlative_graph(x, f, cons, order)
    @test G.fadjlist == [[2], [1, 3], [2]]
end