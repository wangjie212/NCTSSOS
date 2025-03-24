using Test, NCTSSOS
using Clarabel
using Graphs

@testset "Jun Problem" begin
    num_sites = 10 
    g = star_graph(num_sites)

    true_ans = -1.0

    vec_idx2ij = [(i, j) for i in 1:num_sites for j in (i+1):num_sites]

    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

    @ncpolyvar pij[1:length(vec_idx2ij)]

    objective = sum(1.0 * pij[[findvaridx(ee.src, ee.dst) for ee in edges(g)]])

    gs = 
        [
            (
                pij[findvaridx(sort([i, j])...)] * pij[findvaridx(sort([j, k])...)] +
                pij[findvaridx(sort([j, k])...)] * pij[findvaridx(sort([i, j])...)] -
                pij[findvaridx(sort([i, j])...)] - pij[findvaridx(sort([j, k])...)] -
                pij[findvaridx(sort([i, k])...)] + 1.0
            ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
            (i != j && j != k && i != k)
        ]
    

    pop = PolyOpt(objective; constraints=gs, is_equality=[true for _ in gs], is_unipotent=true)

    solver_config = SolverConfig(optimizer=Clarabel.Optimizer; mom_order=1)

    result = cs_nctssos(pop, solver_config)
    @test isapprox(result.objective,true_ans; atol=1e-6)
end

@testset "Problem Creation Interface" begin
    n = 2
    @ncpolyvar x[1:n]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    # change struct name
    pop = PolyOpt(f; constraints=[g, h1], is_equality=[false, true])

    solver_config = SolverConfig(optimizer=Clarabel.Optimizer; mom_order=2, cs_algo=MF(), ts_algo=MMD())

    result = cs_nctssos(pop, solver_config)
    @test isapprox(result.objective, -1.0; atol=1e-4)

    result_higher = cs_nctssos_higher(pop, result, solver_config)
    @test isapprox(result.objective, result_higher.objective; atol=1e-4)
end


@testset "README Example Unconstrained" begin
    @ncpolyvar x[1:3]
    f = 1.0 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2] + x[2]*x[1] + x[2]*x[3] + x[3]*x[2]

    pop =  PolyOpt(f)

    solver_config_dense = SolverConfig(optimizer=Clarabel.Optimizer)

    result_dense = cs_nctssos(pop, solver_config_dense)

    result_cs = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF()))

    @test isapprox(result_dense.objective, result_cs.objective, atol=1e-4)

    result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))

    @test isapprox(result_cs.objective, result_cs_ts.objective, atol=1e-4)
end

@testset "README Example Constrained" begin
    @ncpolyvar x[1:2]
    f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1]*x[2] + x[2]*x[1] - 2.0

    pop =  PolyOpt(f; constraints=[g, h1], is_equality=[false, true])

    result_dense = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer))

    result_cs = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF()))

    @test isapprox(result_dense.objective, result_cs.objective, atol=1e-4)

    result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))

    @test isapprox(result_cs.objective, result_cs_ts.objective, atol=1e-4)

    result_cs_ts_higher = cs_nctssos_higher(pop, result_cs_ts, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))

    @test isapprox(result_dense.objective, result_cs_ts_higher.objective, atol=1e-4)
end
