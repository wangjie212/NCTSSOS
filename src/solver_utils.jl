function remove_zero_degree(m::Monomial{C}) where {C}
    isconstant(m) && return m
    return prod([x[1]^x[2] for x in filter(!(iszero ∘ last), collect(zip(m.vars, m.z)))])
end

function star(m::Monomial{C}) where {C}
    return prod([
        x[1]^x[2] for
        x in filter(!(iszero ∘ last), collect(zip(reverse(m.vars), reverse(m.z))))
    ])
end

function symmetric_canonicalize(monomial::Monomial{C}) where {C}
    return min(remove_zero_degree(monomial), star(monomial))
end

function symmetric_canonicalize(poly::Polynomial{C,T}) where {C,T}
    return mapreduce(
        p -> DP.coefficient(p)' * symmetric_canonicalize(DP.monomial(p)),
        +,
        DP.terms(poly);
        init=zero(poly),
    )
end

#TODO: I did not consider binary variable but it's easy to extend, just filter out in monomial z==2 && vars in binary set 
get_basis(vars::Vector{V}, d::Int) where {V<:AbstractVariable} = mapreduce(cur_d -> remove_zero_degree.(monomials(vars, cur_d)), vcat, 0:d)
# need to remove zero degree other wise sortting fails

support(poly::Polynomial{C,T}, canonicalize::Function) where {C,T} = canonicalize ∘ monomials(poly)
