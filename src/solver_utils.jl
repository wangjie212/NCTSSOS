function remove_zero_degree(m::Monomial{C}) where {C}
    isconstant(m) && return m
    return prod([x[1]^x[2] for x in filter(!(iszero ∘ last), collect(zip(m.vars, m.z)))])
end

function star(m::Monomial{C}) where {C}
    return prod([
        x[1]^x[2] for x in filter(!(iszero ∘ last), reverse(collect(zip(m.vars, m.z))))
    ])
end

function symmetric_canonicalize(monomial::Monomial{C}) where {C}
    return min(remove_zero_degree(monomial), star(monomial))
end

function symmetric_canonicalize(poly::Polynomial{C,T}) where {C,T}
    return mapreduce(
        p -> coefficient(p)' * symmetric_canonicalize(monomial(p)),
        +,
        terms(poly);
        init=zero(poly),
    )
end

#NOTE: I did not consider binary variable but it's easy to extend, just filter out in monomial z==2 && vars in binary set
function get_basis(vars::Vector{PolyVar{C}}, d::Int) where {C}
    # need to remove zero degree other wise sortting fails
    return mapreduce(cur_d -> remove_zero_degree.(monomials(vars, cur_d)), vcat, 0:d)
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
