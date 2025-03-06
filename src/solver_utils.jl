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

#TODO: I did not consider binary variable but it's easy to extend, just filter out in monomial z==2 && vars in binary set 
function get_basis(vars::Vector{PolyVar{C}}, d::Int) where {C}
    return mapreduce(cur_d -> remove_zero_degree.(monomials(vars, cur_d)), vcat, 0:d)
end
# need to remove zero degree other wise sortting fails

function support(poly::Polynomial{C,T}, canonicalize::Function) where {C,T}
    return canonicalize.(monomials(poly))
end
