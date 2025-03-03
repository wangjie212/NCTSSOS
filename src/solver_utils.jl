function remove_zero_degree(m::M) where {M<:AbstractMonomialLike}
    isconstant(m) && return m
    return prod([x[1]^x[2] for x in filter(!(iszero ∘ last), collect(zip(m.vars, m.z)))])
end

star(m::MD) where {MD<:AbstractMonomialLike} = prod([x[1]^x[2] for x in filter(!(iszero ∘ last), collect(zip(reverse(m.vars), reverse(m.z))))])

symmetric_canonicalize(monomial::MD) where {MD<:AbstractMonomialLike} = min(monomial, star(monomial))

symmetric_canonicalize(poly::PD) where {PD<:AbstractPolynomialLike{<:Real}} = mapreduce(p -> coefficient(p)' * symmetric_canonicalize(monomial(p)), +, terms(poly); init=zero(poly))