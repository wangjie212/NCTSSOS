#TODO: add user interface
# learn from OMEinsum.jl
# consider adding Solver Interface
# consider obtaining enough information on Moment matrix etc to check if problem solved correctly

function cs_nctssos(pop::PolynomialOptimizationProblem{C,T}, mom_order::Int; ts_order::Int=0, cs_algo::Union{Nothing,EliminationAlgorithm}=nothing, ts_algo::Union{Nothing,EliminationAlgorithm}=nothing) where {C,T}
    isnothing(ts_algo) && (@assert iszero(ts_order) "No Chordal Decomposition Algorithm provided for Term Sparsity, ts_order needs to be zero")
    @assert mom_order >= maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.constraints]]) "Moment Matrix Order Must be large enough to include all monomials in problem"

    corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, mom_order, cs_algo)

    if isnothing(ts_algo)
        cliques_term_sparsities = [
            [TermSparsity(Vector{Monomial{C}}(), [basis]) for basis in idx_basis]
            for idx_basis in corr_sparsity.cliques_idcs_bases
        ]
    else
        cliques_objective = [reduce(+, [issubset(effective_variables(mono), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

        # prepare the support for each term sparse localizing moment
        initial_activated_supp = [
            # why does order matter here?
            sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, pop.constraints[cons_idx]), [neat_dot(b, b) for b in idcs_bases[1]])
            for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)
        ]

        # need to iterate
        cliques_term_sparsities = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
            [iterate_term_sparse_supp(activated_supp, poly, basis, ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
        end
    end

    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, cliques_term_sparsities)
    sos_problem = sos_dualize(moment_problem)
    return sos_problem
end

struct Answer{F}
    objective::F # support for high precision solution
    # should contain objective and moment matrix or gram matrix for manually checking what happened
end


# FIXME: should it just be `sovle` or too general a name?
function solve_problem(problem::Union{MomentProblem,SOSProblem}, optimizer)
    set_optimizer(problem.model, optimizer)
    optimize!(problem.model)
    return Answer(objective_value(problem.model))
end

function cs_nctssos_higher!()

end
