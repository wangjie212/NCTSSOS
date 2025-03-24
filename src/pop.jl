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

# for user convenience
# be as general as possible
function PolyOpt(objective::Polynomial{C,T}, constraints, is_equality, comm_gp, is_unipotent::Bool, is_projective::Bool, obj_type::Objective) where {C,T}
    @assert !(T <: Integer) "The polynomial coefficients can not be integers (not supported by JuMP solvers)."
    cons = collect(Polynomial{C,T}, constraints)
    comm_vars = collect(PolyVar{C}, comm_gp)
    @assert length(is_equality) == length(cons) "The number of constraints must be the same as the number of equality conditions."
    vars = sorted_union(variables(objective), [variables(c) for c in cons]...)
    @assert issubset(comm_vars, vars) "The commutative variables must be a subset of the variables."
    @assert !(is_unipotent && is_projective) "The problem cannot be both unipotent and projective."
    return PolyOpt{C,T,obj_type}(objective, cons, is_equality, vars, Set(comm_vars), is_unipotent, is_projective)
end

function PolyOpt(objective::Polynomial{C,T}, constraints) where {C,T}
    return PolyOpt(objective, constraints, fill(false, length(constraints)), PolyVar{C}[], false, false, EIGEN)
end

# for user, unconstrained pop
function PolyOpt(objective::Polynomial{C,T}) where {C,T}
    return PolyOpt(objective, Polynomial{C,T}[], Bool[], PolyVar{C}[], false, false, EIGEN)
end
