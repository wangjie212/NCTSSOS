using Test, NCTSSOS

@testset "PolyOpt" begin
    nvars = 10
    ncons = 3
    @ncpolyvar x[1:nvars]
    objective = 1.0 * sum(x .^ 2)
    constraints = [1.0 * sum(i .* x) for i in 1:ncons]

    pop = PolyOpt(objective, constraints)

    @test NCTSSOS.nvariables(pop) == nvars
    @test NCTSSOS.nconstraints(pop) == ncons
    @test NCTSSOS.iscommutative(pop) == false

    pop = PolyOpt(objective, [])

    @test NCTSSOS.nvariables(pop) == nvars
    @test NCTSSOS.nconstraints(pop) == 0
    @test NCTSSOS.iscommutative(pop) == false

    pop = PolyOpt(objective, Set([constraints; sum(x)]))

    @test NCTSSOS.nvariables(pop) == nvars
    @test NCTSSOS.nconstraints(pop) == ncons
    @test NCTSSOS.iscommutative(pop) == false

    pop = PolyOpt(objective)

    @test NCTSSOS.nvariables(pop) == nvars
    @test NCTSSOS.nconstraints(pop) == 0
    @test NCTSSOS.iscommutative(pop) == false
end
