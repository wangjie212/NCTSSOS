using Test, NCTSSOS
using DynamicPolynomials


@testset "Remove Variables with zero degree" begin
    # TODO: should I put variable initialization in outer testset scope?
    @ncpolyvar x y z

    mono1 = Monomial{false}([x, y, z], [2, 0, 1])
    mono1_simp = remove_zero_degree(mono1)
    @test mono1_simp.vars == [x, z]
    @test mono1_simp.z == [2,1]

    mono2 = Monomial{false}([x, y, z], [2, 1, 1])
    mono2_simp = remove_zero_degree(mono2)
    @test mono2_simp.vars == [x, y, z]
    @test mono2_simp.z == [2,1,1]

    mono3 = Monomial{false}([x, y, z], [0, 0, 0])
    mono3_simp = remove_zero_degree(mono3)
    @test mono3_simp == 1
end


@testset "Star Operation" begin
    @ncpolyvar x y z

    mono1 = Monomial{false}([x, y, z], [2, 0, 1])

    # TODO: I am assuming all variables are Hermitian
    mono1_star = star(mono1)

    @test mono1_star.vars == [z, x]
    @test mono1_star.z == [1, 2]

    mono2 = Monomial{false}([x, y, z], [0, 0, 0])
    mono2_star = star(mono2)

    @test mono2_star == 1

    mono3 = x*y*z
    mono3_star = star(mono3)
end

@testset "Symmetric Canonical Form" begin
    @ncpolyvar x y z
    mono1 = x^2*y*z
    mono1_sym = symmetric_canonicalize(mono1)
    @test mono1_sym.vars == [z, y, x]
    @test mono1_sym.z == [1, 1, 2]

    mono2 = z*y*x^2
    mono2_sym = symmetric_canonicalize(mono2)
    @test mono2_sym.vars == [z, y, x]
    @test mono2_sym.z == [1, 1, 2]

    poly1 = 0.1*x^2*y*z + 0.2*z*y*x^2 + 0.3
    poly1_sym = symmetric_canonicalize(poly1)

    @test coefficients(poly1_sym) ≈ [0.3, 0.3]
    @test monomials(poly1_sym) == [z*y*x^2, 1]
end

@testset "Get basis" begin
    @ncpolyvar x y z

    nc_basis_deg2 = get_basis([x,y,z],2)

    @test sort(nc_basis_deg2) == sort([one(x), x, y, z, x^2, y^2, z^2, x * y, x * z, y * z, z * x, z * y, y * x])

    @polyvar x y z

    comm_basis_deg2 = get_basis([x,y,z],2)

    @test sort(comm_basis_deg2) == sort([one(x), x, y, z, x^2, y^2, z^2, x * y, x * z, y * z])
end

@testset "get support" begin

    @ncpolyvar x y z

    poly = 0.1 * x^2 * y + 0.2 * x - 0.1

    @test sort(support(poly, identity)) == sort([1, x, x^2 * y])
end