abstract type OptimizationProblem end

@enum Objective EIGEN TRACE STATE
# C: false if the variables are not commutative
# T: type of the coefficients
# OBJ: type of the objective function, one of `Objective`
struct PolyOpt{C,T,OBJ} <: OptimizationProblem
    objective::Polynomial{C,T}
    constraints::Vector{Polynomial{C,T}} # NOTE: assuming constraints are all simplified using comm_gp, is_unipotent, and is_projective
    is_equality::Vector{Bool}
    variables::Vector{PolyVar{C}}
    comm_gp::Set{PolyVar{C}} # Set of variables that commutes with variables not in the set
    is_unipotent::Bool # square to 1. Examples: Pauli Operators, SWAP Operators
    is_projective::Bool # X^2 = X. Is projective.
end

# FIXME: is this default parameters make sense?
function PolyOpt(objective::Polynomial{C,T}; constraints=Any[], is_equality=fill(false, length(constraints)), comm_gp=PolyVar{C}[], is_unipotent::Bool=false, is_projective::Bool=false, obj_type::Objective=EIGEN) where {C,T}
    @assert !(T <: Integer) "The polynomial coefficients can not be integers (not supported by JuMP solvers)."
    cons = collect(Polynomial{C,T}, constraints)
    is_eq = collect(Bool, is_equality)
    @assert length(is_eq) == length(cons) "The number of constraints must be the same as the number of equality conditions."
    vars = sorted_union(variables(objective), [variables(c) for c in cons]...)
    @assert issubset(comm_gp, vars) "The commutative variables must be a subset of the variables."
    @assert !(is_unipotent && is_projective) "The problem cannot be both unipotent and projective."
    return PolyOpt{C,T,obj_type}(objective, cons, is_eq, vars, Set(comm_gp), is_unipotent, is_projective)
end
