# NOTE: for ncpolyvar generating monomials from `monomials` may give you x^0*y^2*z^0*x^1 terms
# which is not equal to y^2*x^1
# plue this bug: https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/118#issue-1412618512
# I will keep this function here
function remove_zero_degree(m::Monomial)
    isconstant(m) && return m
    return prod([x[1]^x[2] for x in filter(!(iszero ∘ last), collect(zip(m.vars, m.z)))])
end

function star(m::Monomial)
    return prod([x[1]^x[2] for x in filter(!(iszero ∘ last), reverse(collect(zip(m.vars, m.z))))])
end

function symmetric_canonicalize(monomial::Monomial)
    return min(monomial, star(monomial))
end

function symmetric_canonicalize(poly::Polynomial)
    return mapreduce(p -> coefficient(p)' * symmetric_canonicalize(monomial(p)), +, terms(poly); init=zero(poly))
end

function cyclic_canonicalize(monomial::Monomial)
    isconstant(monomial) && return monomial
    return minimum(mapreduce(vcat, 0:length(monomial.vars)-1) do shift
        shifted_mono = prod([var^expo for (var, expo) in circshift!(collect(zip(monomial.vars, monomial.z)), shift)])
        [shifted_mono, star(shifted_mono)]
    end)
end

function cyclic_canonicalize(poly::Polynomial)
    return mapreduce(p -> coefficient(p) * cyclic_canonicalize(monomial(p)), +, terms(poly); init=zero(poly))
end

function get_basis(vars::Vector{PolyVar{C}}, d::Int) where {C}
    # need to remove zero degree other wise sortting fails
    # return mapreduce(cur_d -> remove_zero_degree.(monomials(vars, cur_d)), vcat, 0:d)
    return remove_zero_degree.(sort(monomials(vars, 0:d)))
end

function support(poly::Polynomial{C,T}, canonicalize::Function) where {C,T}
    return canonicalize.(monomials(poly))
end

function neat_dot(x::Monomial{C}, y::Monomial{C}) where {C}
    # NOTE: the `*` in DynamicPolynomials sometimes creates monomials with degree 0, which we don't want
    return remove_zero_degree(star(x) * y)
end

sorted_unique(xs) = sort(unique(xs))
sorted_union(xs...) = sort(union(xs...))

get_dim(cons::VectorConstraint) = cons.set isa MOI.PositiveSemidefiniteConeSquare ? JuMP.shape(cons).side_dimension : JuMP.shape(cons).dims[1]

function _comm(mono::Monomial{C}, comm_gp::Set{PolyVar{C}}) where {C}
    return prod(zip(mono.vars, mono.z)) do (var, expo)
        var in comm_gp ? var^expo : var^(zero(expo))
    end *
           prod(zip(mono.vars, mono.z)) do (var, expo)
        !(var in comm_gp) ? var^expo : var^(zero(expo))
    end
end

function _unipotent(mono::Monomial)
    prev_mono = mono
    local cur_mono
    while true
        cur_mono = prod(zip(prev_mono.vars, prev_mono.z)) do (var, expo)
            var^(expo % 2)
        end
        cur_mono == prev_mono && break
        prev_mono = cur_mono
    end
    return cur_mono
end

_projective(mono::Monomial) =
    prod(zip(mono.vars, mono.z)) do (var, expo)
        var^(iszero(expo) ? expo : one(expo))
    end

function reducer(pop::PolyOpt)
    function (x)
        cx = _comm(x, pop.comm_gp)
        return pop.is_unipotent ? _unipotent(cx) : (pop.is_projective ? _projective(cx) : cx)
    end
end
