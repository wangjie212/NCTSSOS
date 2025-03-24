using BenchmarkTools, Profile, PProf
using NCTSSOS, Test

using DynamicPolynomials
@testset "PolyOpt Construction" begin
    n = 300
    @ncpolyvar x[1:n]
    f = 0.0
    for i = 1:n
        jset = max(1, i - 5):min(n, i + 1)
        jset = setdiff(jset, i)
        g = sum(x[j] + x[j]^2 for j in jset)
        f += (2 * x[i] + 5 * x[i]^3 + 1 - g)^2
    end

    cons = Polynomial{false,Float64}[1.0*x[i]*x[mod1(i+5,n)] for i in 1:n]
    is_equality = fill(false, length(cons))

    @benchmark PolyOpt($f, $cons, $is_equality)

    Profile.Allocs.clear()
    Profile.Allocs.@profile PolyOpt(f,cons, is_equality);

    prof = Profile.Allocs.fetch();

    PProf.Allocs.pprof(prof,from_c=false)

    t = @benchmark PolyOpt($f, $cons, $is_equality)

    @test t.allocs == 0.0
end
