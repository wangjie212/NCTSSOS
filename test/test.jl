using DynamicPolynomials
using NCTSSOS
using Test

@testset "poly info" begin
    n = 3
    @ncpolyvar x[1:3]
    f = x[1]^2-x[1]*x[2]-x[2]*x[1]+3x[2]^2-2x[1]*x[2]*x[1]+2x[1]*x[2]^2*x[1]-x[2]*x[3]-x[3]*x[2]+
    6x[3]^2+9x[2]^2*x[3]+9x[3]*x[2]^2-54x[3]*x[2]*x[3]+142x[3]*x[2]^2*x[3]
    n,supp,coe = NCTSSOS.poly_info(f, x)
    @test n == 3  # number of variables
    @test supp == [[1, 2, 2, 1], [3, 2, 2, 3], [1, 2, 1], [2, 2, 3], [3, 2, 2], 
                   [3, 2, 3], [1, 1], [1, 2], [2, 2], [2, 1],
                   [2, 3], [3, 3], [3, 2]]
    @test coe == [2, 142, -2, 9, 9, -54, 1, -1, 3, -1, -1, 6, -1] 

    # reduce the basis by reflection symmetry - used in the eigen objective
    # eigen(A) = eigen(A')
    sym_supp, sym_coe = NCTSSOS.sym_canon(supp, coe)
    expected_sym_supp = [[1, 2, 2, 1], [3, 2, 2, 3], [1, 2, 1], [2, 2, 3], 
                   [3, 2, 3], [1, 1], [1, 2], [2, 2],
                   [2, 3], [3, 3]]
    expected_sym_coe = [2, 142, -2, 18,
                     -54, 1, -2, 3,
                     -2, 6, -1]
    order = sortperm(expected_sym_supp)
    @test sym_supp == expected_sym_supp[order]
    @test Int.(sym_coe) == expected_sym_coe[order]


    # reduce the basis by cyclic symmetry - used in the trace objective
    # trace(ABC) = trace(BCA) = trace(CAB), trace(A) = trace(A')
    append!(supp, [[1, 1, 2, 2], [1,2,3], [3,2,1]])
    append!(coe, [5, 5, 5])
    cyc_supp, cyc_coe = NCTSSOS.cyclic_canon(supp, coe)
    expected_cyc_supp = map(x->minimum(circshift(x, i) for i in 0:length(x)-1), [[1, 2, 2, 1], [3, 2, 2, 3], [1, 2, 1], [2, 2, 3], 
                   [3, 2, 3], [1, 1], [1, 2], [2, 2],
                   [2, 3], [3, 3], [1, 2, 3]])
    expected_cyc_coe = [7, 142, -2, 18,
                       -54, 1, -2, 3,
                       -2, 6, 10] 
    order = sortperm(expected_cyc_supp)
    @test map(x->Int.(x), cyc_supp) == expected_cyc_supp[order]
    @test Int.(cyc_coe) == expected_cyc_coe[order]
end

@testset "newton_ncbasis" begin
    supp = [[1, 2, 2, 1], [3, 2, 2, 3], [1, 2, 1], [2, 2, 3], 
                   [3, 2, 3], [1, 1], [1, 2], [2, 2],
                   [2, 3], [3, 3], [1, 2, 3]]
    supp = map(x->UInt16.(x), supp)
    # what if the basis is not even and symmetric?
    ncbasis = NCTSSOS.newton_ncbasis(supp)
    @test ncbasis == [Int[], [1], [2], [2, 1], [2, 3], [3]]

    # TODO: fix this, needs help
    # @test_broken cycbasis = NCTSSOS.newton_cyclic(supp, 3, 1)
    cycbasis = NCTSSOS.get_ncbasis(3, 2, ind=UInt16.(1:3))
    @test sort(cycbasis) == sort([Int[], [1], [2], [3], [1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]])
    # binary means x^2 = 1
    cycbasis = NCTSSOS.get_ncbasis(3, 2, ind=UInt16.(1:3), binary=true)
    @test sort(cycbasis) == sort([Int[], [1], [2], [3], [1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2]])

    cycbasis = NCTSSOS.get_ncbasis(3, 2, ind=UInt16[7, 8, 9])
    @test sort(cycbasis) == sort([Int[], [7], [8], [9], [7, 7], [7, 8], [7, 9], [8, 7], [8, 8], [8, 9], [9, 7], [9, 8], [9, 9]])
    # Q: does the order of the basis matter?
end

@testset "get_blocks" begin
    supp = [[1, 2, 2, 1], [3, 2, 2, 3], [1, 2, 1], [2, 2, 3], 
                   [3, 2, 3], [1, 1], [1, 2], [2, 2],
                   [2, 3], [3, 3], [1, 2, 3]]
    supp = map(x->UInt16.(x), supp)
    coe = ones(length(supp))
    supp, coe = NCTSSOS.sym_canon(supp, coe)
    basis = NCTSSOS.newton_ncbasis(supp)
    blocks,cl,blocksize = NCTSSOS.get_blocks(supp, basis, TS="MD", QUIET=true, merge=false, md=3, obj="eigen", partition=0, constraint=nothing)
    @test blocks == [[0x0003, 0x0004, 0x0005, 0x0006], [0x0002, 0x0003, 0x0004, 0x0005], [0x0001, 0x0004, 0x0005]]
end

@testset "case 1" begin
    n = 3
    @ncpolyvar x[1:3]
    f = x[1]^2-x[1]*x[2]-x[2]*x[1]+3x[2]^2-2x[1]*x[2]*x[1]+2x[1]*x[2]^2*x[1]-x[2]*x[3]-x[3]*x[2]+
    6x[3]^2+9x[2]^2*x[3]+9x[3]*x[2]^2-54x[3]*x[2]*x[3]+142x[3]*x[2]^2*x[3]

    opt,data = nctssos_first(f, x, newton=true, reducebasis=true, TS="MD", obj="eigen", QUIET=true, solver="COSMO")
    @test isapprox(opt, -0.0035512, atol=1e-6)

end

@testset "case 2" begin
    n = 2
    @ncpolyvar x[1:2]
    f = 2-x[1]^2+x[1]*x[2]^2*x[1]-x[2]^2
    g = 4-x[1]^2-x[2]^2
    h = x[1]*x[2]+x[2]*x[1]-2
    pop = [f, g, h]

    opt,data = nctssos_first(pop, x, 2, numeq=1, TS="MD", obj="eigen", QUIET=true, solver="COSMO")
    @test isapprox(opt, -1, atol=1e-6)

    opt,data = nctssos_higher!(data, TS="MD", QUIET=true, solver="COSMO")
    @test isapprox(opt, -1, atol=1e-6)

    opt,data = nctssos_first(pop, x, 2, numeq=1, TS="MD", obj="trace", QUIET=true, solver="COSMO")
    @test isapprox(opt, -1, atol=1e-6)


    opt,data = nctssos_higher!(data, TS="MD", QUIET=true, solver="COSMO")
    @test isapprox(opt, -1, atol=1e-6)
end

@testset "case 3" begin 
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

    opt,data = cs_nctssos_first([f], x, 3, TS="MD", obj="trace", QUIET=true, solver="COSMO")
    # TODO: verify the optimum result is 0.0
    @test isapprox(opt, 0.0, atol=1e-4)
    # optimum = 6.5e-8

    # NOTE: following blocks are commented out because they are too slow

    # opt,data = cs_nctssos_higher!(data, TS="MD", QUIET=true, solver="COSMO")
    # @test isapprox(opt, 6.5e-7, atol=1e-10)
    # optimum = 6.5e-7

    # pop = [f]
    # for i = 1:n
    #     push!(pop, 1-x[i]^2)
    #     push!(pop, x[i]-1/3)
    # end

    # opt,data = cs_nctssos_first(pop, x, 3, TS="MD", obj="eigen", QUIET=true, solver="COSMO")
    # @test isapprox(opt, 3.011288, atol=1e-10)
    # # optimum = 3.011288

    # opt,data = cs_nctssos_first(pop, x, 3, TS="MD", obj="trace", QUIET=true, solver="COSMO")
    # @test isapprox(opt, 3.011288, atol=1e-10)
    # optimum = 3.011288
end
