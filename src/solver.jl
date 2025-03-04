abstract type AbstractMethod end

struct MomentMethod <: AbstractMethod 
	order::Int
	clique_func #function to get clique
end

get_order(method::MomentMethod) = method.order

get_clique_func(method::MomentMethod) = method.clique_func

struct SOSMethod <: AbstractMethod

end

function init_moment_vector!(model::Model, basis::VM) where {VM<:AbstractVector{<:AbstractMonomialLike}}
	moment_vector = @variable(model, y[1:length(basis)])
	@constraint(model, first(y) == 1)
	return Dict(basis .=> moment_vector)
end

replace_variable_with_jump_variables(poly::PD, total_basis2var_dict::Dict{M,VRef}) where {PD<:AbstractPolynomialLike{<:Real},M<:AbstractMonomialLike,VRef<:AbstractVariableRef} = mapreduce(x -> DP.coefficient(x) * total_basis2var_dict[remove_zero_degree(DP.monomial(x))], +, terms(poly))


function constrain_moment_matrix!(model::Model, poly::PD, local_basis::VM, basis2var_dict::Dict{M,VRef}) where {PD<:AbstractPolynomialLike{<:Real},VM<:AbstractVector{<:AbstractMonomialLike},M<:AbstractMonomialLike,VRef<:AbstractVariableRef}
    moment_mtx = [replace_variable_with_jump_variables(poly * row_idx * col_idx, basis2var_dict) for row_idx in star.(local_basis), col_idx in local_basis]
	return @constraint(model, moment_mtx in PSDCone())
end

function make_sdp(method::MomentMethod, pop::PolynomialOptimizationProblem)
    objective = pop |> symmetric_canonicalize ∘ get_objective

    model = Model() # TODO: I should leave option to user for which sovler to use right?

	# total_basis that is the union of support of all constraints and objective and moment matrix and localizing matrices
    total_basis = begin
        moment_mtx_idcs = get_basis(get_variables(pop), get_order(method))
        sort(unique([remove_zero_degree(row_idx * col_idx) for row_idx in star.(moment_mtx_idcs) for col_idx in moment_mtx_idcs]))
    end

	total_basis2var_dict = init_moment_vector!(model, total_basis)

    constraint_matrices = [constrain_moment_matrix!(model, cur_poly, get_basis(get_variables(pop), get_order(method) - ceil(Int, maxdegree(cur_poly) / 2)), total_basis2var_dict) for cur_poly in vcat(get_constraints(pop)..., one(objective))]


	# get support of obj and constriants
	# store in supp

	# get basis for moment matrix and localizing matrices
	# store in basis a vector of basis monomials

    # ksupp should contain support of obj and all constraints + support of moment matrix and localizing matrices

    @objective(model, Min, replace_variable_with_jump_variables(objective, total_basis2var_dict))
	return model
end



