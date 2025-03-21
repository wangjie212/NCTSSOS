using Test, NCTSSOS
using Clarabel

# API 
# Performance
# Correctness

@testset "Problem Creation Interface" begin
    n = 2
    @ncpolyvar x[1:n]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    h2 = -h1
    # change struct name
    pop = PolyOpt(f, [g, h1, h2])

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
    h2 = -h1

    pop =  PolyOpt(f, [g, h1, h2])

    result_dense = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer))

    result_cs = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF()))

    @test isapprox(result_dense.objective, result_cs.objective, atol=1e-4)

    result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))

    @test isapprox(result_cs.objective, result_cs_ts.objective, atol=1e-4)

    result_cs_ts_higher = cs_nctssos_higher(pop, result_cs_ts, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))

    @test isapprox(result_dense.objective, result_cs_ts_higher.objective, atol=1e-4)
end
