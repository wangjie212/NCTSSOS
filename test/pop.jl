using Test, NCTSSOS

@testset "PolyOpt Constructor" begin
    nvars = 10
    ncons = 3
    @ncpolyvar x[1:nvars]
    objective = 1.0 * sum(x .^ 2)
    constraints = [1.0 * sum(i .* x) for i in 1:ncons]

    @testset "Unconstrained" begin
        pop = PolyOpt(objective, [])

        @test NCTSSOS.nvariables(pop) == nvars
        @test NCTSSOS.nconstraints(pop) == 0
        @test NCTSSOS.iscommutative(pop) == false
        @test pop.is_equality == Bool[]

        pop = PolyOpt(objective)

        @test NCTSSOS.nvariables(pop) == nvars
        @test NCTSSOS.nconstraints(pop) == 0
        @test NCTSSOS.iscommutative(pop) == false
        @test pop.is_equality == Bool[]
    end

    @testset "Constraints All Inequality" begin
        pop = PolyOpt(objective, constraints)

        @test NCTSSOS.nvariables(pop) == nvars
        @test NCTSSOS.nconstraints(pop) == ncons
        @test NCTSSOS.iscommutative(pop) == false
        @test pop.is_equality == fill(false, ncons)

        pop = PolyOpt(objective, Set([constraints; sum(x)]))

        @test NCTSSOS.nvariables(pop) == nvars
        @test NCTSSOS.nconstraints(pop) == ncons
        @test NCTSSOS.iscommutative(pop) == false
        @test pop.is_equality == fill(false, ncons)
    end

    @testset "Constraints Mixed Equality and Inequality" begin
        is_equality = [isodd(i) ? true : false for i in 1:ncons]
        pop = PolyOpt(objective, constraints, is_equality)

        @test NCTSSOS.nvariables(pop) == nvars
        @test NCTSSOS.nconstraints(pop) == ncons
        @test NCTSSOS.iscommutative(pop) == false
        @test pop.is_equality == is_equality

        pop = PolyOpt(objective, Set([constraints; sum(x)]), is_equality)

        @test NCTSSOS.nvariables(pop) == nvars
        @test NCTSSOS.nconstraints(pop) == ncons
        @test NCTSSOS.iscommutative(pop) == false
        @test pop.is_equality == is_equality 

        @test_throws AssertionError PolyOpt(objective, constraints, fill(true, ncons + 1))
    end
end
