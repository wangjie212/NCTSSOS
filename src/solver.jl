abstract type AbstratMethod end

struct MomentMethod <: AbstratMethod 
	order::Int
	clique_func #function to get clique
end

get_order(method::MomentMethod) = method.order

get_clique_func(method::MomentMethod) = method.clique_func

struct SOSMethod <: AbstractMethod

end

function init_moment_vector!(model::Model, basis::VM ) where {VM<:AbstractVector{<:AbstractMonomialLike}}
	moment_vector = @variable(model, y[1:length(basis)])
	return Dict(basis .=> moment_vector)
end

function make_moment_matrix(model::Model, poly::PD, basis2var_dict::Dict{VM,VRef}) where {PD<:AbstractPolynomialLike{<:Real},VM<:AbstractVector{<:AbstractMonomialLike},VRef<:AbstractVariableRef}

end

function make_sdp(method::MomentMethod, pop::PolynomialOptimizationProblem)
    objective = pop |> symmetric_canonicalize ∘ get_objective

    model = Model() # TODO: I should leave option to user for which sovler to use right?

	basis2var_dict = init_moment_vector!(model, basis)

	# get support of obj and constriants
	# store in supp

	# get basis for moment matrix and localizing matrices
	# store in basis a vector of basis monomials

    # ksupp should contain support of obj and all constraints + support of moment matrix and localizing matrices

	return model
end



