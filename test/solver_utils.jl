using Test, NCTSSOS
using DynamicPolynomials
using JuMP
using NCTSSOS: remove_zero_degree, star, symmetric_canonicalize, get_basis, support, neat_dot, get_dim, _comm, _unipotent, _projective, reduce!, cyclic_canonicalize

@testset "Utilities" begin
    @ncpolyvar x y z

    @testset "Remove Variables with zero degree" begin
        mono1 = Monomial{false}([x, y, z], [2, 0, 1])
        mono1_simp = remove_zero_degree(mono1)
        @test mono1_simp.vars == [x, z]
        @test mono1_simp.z == [2, 1]

        mono2 = Monomial{false}([x, y, z], [2, 1, 1])
        mono2_simp = remove_zero_degree(mono2)
        @test mono2_simp.vars == [x, y, z]
        @test mono2_simp.z == [2, 1, 1]

        mono3 = Monomial{false}([x, y, z], [0, 0, 0])
        mono3_simp = remove_zero_degree(mono3)
        @test mono3_simp == 1
    end

    @testset "Star Operation" begin
        mono1 = Monomial{false}([x, y, z], [2, 0, 1])

        # NOTE: I am assuming all variables are Hermitian
        mono1_star = star(mono1)

        @test mono1_star == z * x^2

        mono2 = Monomial{false}([x, y, z], [0, 0, 0])
        mono2_star = star(mono2)

        @test mono2_star == 1

        mono3 = x * y * z
        mono3_star = star(mono3)
        @test mono3_star == z * y * x
    end

    @testset "Symmetric Canonical Form" begin
        mono1 = x^2 * y * z
        mono1_sym = symmetric_canonicalize(mono1)
        @test mono1_sym.vars == [z, y, x]
        @test mono1_sym.z == [1, 1, 2]

        mono2 = z * y * x^2
        mono2_sym = symmetric_canonicalize(mono2)
        @test mono2_sym.vars == [z, y, x]
        @test mono2_sym.z == [1, 1, 2]

        poly1 = 0.1 * x^2 * y * z + 0.2 * z * y * x^2 + 0.3
        poly1_sym = symmetric_canonicalize(poly1)

        @test coefficients(poly1_sym) ≈ [0.3, 0.3]
        @test monomials(poly1_sym) == [z * y * x^2, 1]

        n = 3
        @ncpolyvar a[1:n]
        poly3 =
            a[1]^2 - a[1] * a[2] - a[2] * a[1] + 3a[2]^2 - 2a[1] * a[2] * a[1] +
            2a[1] * a[2]^2 * a[1] - a[2] * a[3] - a[3] * a[2] +
            6a[3]^2 +
            9a[2]^2 * a[3] +
            9a[3] * a[2]^2 - 54a[3] * a[2] * a[3] + 142a[3] * a[2]^2 * a[3]

        supp = [
            [1, 2, 2, 1],
            [3, 2, 2, 3],
            [1, 2, 1],
            [3, 2, 2],
            [3, 2, 3],
            [1, 1],
            [2, 1],
            [2, 2],
            [3, 2],
            [3, 3],
        ]

        coe = [2, 142, -2, 18, -54, 1, -2, 3, -2, 6]

        poly3_sym = symmetric_canonicalize(poly3)

        @test poly3_sym == mapreduce(v -> v[1] * reduce(*, a[v[2]]), +, zip(coe, supp))
    end

    @testset "Cyclic Canonical Form" begin
        @test cyclic_canonicalize(x * y^2 * x) == y^2 * x^2
        @test cyclic_canonicalize(y*y*z) == z*y^2
        @test cyclic_canonicalize(z*z) == z^2
        @test cyclic_canonicalize(one(x)) == 1

        n = 3
        @ncpolyvar a[1:n]
        f = a[1]^2 - a[1] * a[2] - a[2] * a[1] + 3a[2]^2 - 2a[1] * a[2] * a[1] + 2a[1] * a[2]^2 * a[1] - a[2] * a[3] - a[3] * a[2] +
            6a[3]^2 + 9a[2]^2 * a[3] + 9a[3] * a[2]^2 - 54a[3] * a[2] * a[3] + 142a[3] * a[2]^2 * a[3] + 5 * a[1]^2 * a[2]^2 + 5 * a[1] * a[2] * a[3] + 5 * a[3] * a[2] * a[1]
        @test cyclic_canonicalize(f) == cyclic_canonicalize(7 * a[1] * a[2]^2 * a[1] + 142 * a[3] * a[2]^2 * a[3] - 2 * a[1] * a[2] * a[1] + 18 * a[2]^2 * a[3] - 54 * a[3] * a[2] * a[3] + 1 * a[1]^2 - 2 * a[1] * a[2] + 3 * a[2]^2 - 2 * a[2] * a[3] + 6 * a[3]^2 + 10 * a[1] * a[2] * a[3])
    end

    @testset "Get basis" begin
        nc_basis_deg2 = get_basis([x, y, z], 2)

        @test sort(nc_basis_deg2) ==
              sort([one(x), x, y, z, x^2, y^2, z^2, x * y, x * z, y * z, z * x, z * y, y * x])

        @polyvar a b c

        comm_basis_deg2 = get_basis([a, b, c], 2)

        @test sort(comm_basis_deg2) ==
              sort([one(a), b, c, a, a^2, b^2, c^2, a * b, a * c, b * c])
    end

    @testset "get support" begin
        poly = 0.1 * x^2 * y + 0.2 * x - 0.1

        @test sort(support(poly, identity)) == sort([one(x), x, x^2 * y])
    end

    @testset "neat_dot" begin
        mono1 = Monomial{false}([x, y], [1, 0])

        mono2 = Monomial{false}([x, y], [1, 1])

        @test neat_dot(mono1, mono2) == Monomial{false}([x, y], [2, 1])
    end

    @testset "VectorConstraint Dim" begin
        model = Model()
        n = 5
        var1 = @variable(model, [1:n, 1:n])
        var2 = @variable(model, [1:2*n, 1:2*n])

        cons1 = @constraint(model, var1 in PSDCone())
        cons2 = @constraint(model, var2 in Zeros())

        @test get_dim(constraint_object(cons1)) == n
        @test get_dim(constraint_object(cons2)) == 2 * n
    end

    @testset "_comm" begin
        @ncpolyvar a[1:3]
        @ncpolyvar b[1:3]

        mono = a[1]^2*b[2]^2*a[2]*b[1]^3
        @test _comm(mono,a) == a[1]^2*a[2]*b[2]^2*b[1]^3
        mono = a[1]^3*a[3]
        @test _comm(mono,a) == a[1]^3*a[3]
    end

    @testset "_projective" begin
        mono = y * x^3 * y*z^3
        @test _projective(mono) == y*x*y*z
    end

    @testset "_unipotent" begin
        mono = z*x*y*z^2*y*x*z
        @test _unipotent(mono) == one(mono)
    end

    @testset "reduce!" begin
        basis = get_basis([x,y,z],3)
        for b in basis
            @show b
        end

        reduce!(basis, [x], _unipotent)
        @test sort(basis) == sort([one(x*y),z,y,x,z*y,x*y,y*z,x*z,x*z*y,z*y*z,y*z*y,x*y*z])

        reduce!(basis, [x], _projective)
        @test sort(basis) == sort([one(x*y*z),z,y,x,z*y,x*z,x*y,y*z,x*z*y,z*y*z,x*y*z,y*z*y])
    end
end
