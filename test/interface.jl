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

    mom_order = 2
    ts_order = 1
    cs_algo = MF()
    ts_algo = MMD()

    problem = cs_nctssos(pop, mom_order, ts_order, cs_algo, ts_algo)

    myans = solve_problem(problem,Clarabel.Optimizer)
    @test isapprox(myans.objective, -1.0, atol=1e-4)
end
