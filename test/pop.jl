using Test, NCTSSOS
using DynamicPolynomials

@testset "PolynomialOptimizationProblem" begin
    nvars = 10
    ncons = 3
    @polyvar x[1:nvars]
    objective = sum(x .^ 2)
    constraints = [sum(i .* x) for i in 1:ncons]

    pop = PolynomialOptimizationProblem(objective, constraints, x)

    @test NCTSSOS.nvariables(pop) == nvars
    @test NCTSSOS.nconstraints(pop) == ncons
    @test NCTSSOS.iscommutative(pop) == false

    pop = PolynomialOptimizationProblem(objective, x)

    @test NCTSSOS.nvariables(pop) == nvars
    @test NCTSSOS.nconstraints(pop) == 0
    @test NCTSSOS.iscommutative(pop) == false
end
