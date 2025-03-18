using Test, NCTSSOS
using DynamicPolynomials
using CliqueTrees
using Clarabel

@testset "Problem Creation Interface" begin
    n = 2
    @ncpolyvar x[1:n]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    h2 = -h1
    pop = PolynomialOptimizationProblem(f, [g, h1, h2])

    myans = cs_nctssos(pop; optimizer=Clarabel.Optimizer, mom_order=2, cs_algo=MF(), ts_algo=MMD())
    @test isapprox(myans.objective, -1.0, atol=1e-4)
end

@testset "README Example Unconstrained" begin
    @ncpolyvar x[1:3]
    f = 1.0 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2] + x[2]*x[1] + x[2]*x[3] + x[3]*x[2]

    pop = PolynomialOptimizationProblem(f)

    problem = cs_nctssos(pop; mom_order=2)
    myans = solve_problem(problem,Clarabel.Optimizer)

    problem = cs_nctssos(pop; mom_order=2, cs_algo=MF())

    myans_cs = solve_problem(problem, Clarabel.Optimizer)

    @test isapprox(myans.objective, myans_cs.objective, atol=1e-4)

    problem = cs_nctssos(pop; mom_order=2, ts_order=1, ts_algo=MMD())
    myans_ts = solve_problem(problem, Clarabel.Optimizer)

    @test isapprox(myans.objective, myans_ts.objective, atol=1e-4)
end

@testset "README Example Constrained" begin
    @ncpolyvar x[1:2]
    f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1]*x[2] + x[2]*x[1] - 2.0
    h2 = -h1

    pop = PolynomialOptimizationProblem(f, [g, h1, h2])

    problem = cs_nctssos(pop; mom_order=2)
    myans = solve_problem(problem,Clarabel.Optimizer)

    problem = cs_nctssos(pop; mom_order=2, cs_algo=MF())
    myans_cs = solve_problem(problem,Clarabel.Optimizer)
    @test isapprox(myans.objective, myans_cs.objective, atol=1e-4)

    problem = cs_nctssos(pop; mom_order=2, cs_algo=MF(), ts_order=1, ts_algo=MMD())
    myans_csts = solve_problem(problem, Clarabel.Optimizer)
    @test isapprox(myans.objective, myans_csts.objective, atol=1e-4)
end
