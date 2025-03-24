using Test, NCTSSOS
using DynamicPolynomials

@testset "PolyOpt Constructor" begin
    nvars = 10
    ncons = 3
    @ncpolyvar x[1:nvars]
    objective = 1.0 * sum(x .^ 2)
    constraints = [1.0 * sum(i .* x) for i in 1:ncons]

    @testset "Unconstrained" begin
        pop = PolyOpt(objective)

        @test pop.is_equality == Bool[]
        @test sort(pop.variables) == sort(x)
        @test pop.comm_gp == Set{PolyVar{false}}()
        @test !pop.is_unipotent 
        @test !pop.is_projective

        pop = PolyOpt(objective; comm_gp=[x[1]], obj_type=NCTSSOS.TRACE)

        @test pop.comm_gp == Set([x[1]])
        @test pop isa PolyOpt{false,Float64,NCTSSOS.TRACE} 

    end

    @testset "Constrainted Optimization Problem" begin
        pop = PolyOpt(objective; constraints=constraints)

        @test pop.constraints == constraints
        @test pop.is_equality == fill(false, ncons)

        pop = PolyOpt(objective; constraints=Set([constraints; sum(x)]))

        @test length(pop.constraints) == ncons
        @test pop.is_equality == fill(false, ncons)

        is_equality = [isodd(i) ? true : false for i in 1:ncons]
        pop = PolyOpt(objective; constraints=constraints, is_equality=is_equality,is_unipotent=false,is_projective=true)
        @test pop.is_unipotent == false
        @test pop.is_projective == true
        @test pop.is_equality == is_equality
    end

    @testset "Invalid Input" begin
        @test_throws AssertionError PolyOpt(objective; is_equality= [true])
        @test_throws AssertionError PolyOpt(objective; constraints= constraints, is_equality=fill(true, ncons + 1), is_unipotent=false, is_projective=false)
        @test_throws AssertionError PolyOpt(objective; constraints= constraints, is_unipotent=true, is_projective=true)
        @test_throws AssertionError PolyOpt(x[1]*x[2]+x[3]*x[2])
        @ncpolyvar y[1:nvars]
        @test_throws AssertionError PolyOpt(objective; comm_gp=y)
    end
end
