function remove_zero_degree(m::M) where {M<:AbstractMonomialLike}
    isconstant(m) && return m
    return prod([x[1]^x[2] for x in filter(!(iszero ∘ last), collect(zip(m.vars, m.z)))])
end

star(m::MD) where {MD<:AbstractMonomialLike} = prod([x[1]^x[2] for x in filter(!(iszero ∘ last), collect(zip(reverse(m.vars), reverse(m.z))))])

symmetric_canonicalize(monomial::MD) where {MD<:AbstractMonomialLike} = min(monomial, star(monomial))

symmetric_canonicalize(poly::PD) where {PD<:AbstractPolynomialLike{<:Real}} = mapreduce(p -> coefficient(p)' * symmetric_canonicalize(monomial(p)), +, terms(poly); init=zero(poly))

#TODO: I did not consider binary variable but it's easy to extend, just filter out in monomial z==2 && vars in binary set 
get_basis(vars::VV, d::T) where {VV<:AbstractVector{<:AbstractVariable},T<:Integer} =
    mapreduce(cur_d -> remove_zero_degree.(monomials(vars, cur_d)), vcat, 0:d) # need to remove zero degree other wise sortting fails