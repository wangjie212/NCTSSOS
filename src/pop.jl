abstract type OptimizationProblem end

# C: false if the variables are not commutative
# T: type of the coefficients
struct PolyOpt{C,T} <: OptimizationProblem
    objective::Polynomial{C,T}
    constraints::Vector{Polynomial{C,T}} # NOTE: assuming constraints are all simplified using comm_gp, is_unipotent, and is_projective
    is_equality::Vector{Bool}
    variables::Vector{PolyVar{C}}
    comm_gp::Vector{PolyVar{C}}
    is_unipotent::Bool
    is_projective::Bool
end

# for user convenience
# be as general as possible
function PolyOpt(objective::Polynomial{C,T}, constraints, is_equality, comm_gp::Vector{PolyVar{C}}, is_unipotent::Bool, is_projective::Bool) where {C,T}
    @assert !(T <: Integer) "The polynomial coefficients can not be integers (not supported by JuMP solvers)."
    cons = collect(Polynomial{C,T}, constraints)
    @assert length(is_equality) == length(cons) "The number of constraints must be the same as the number of equality conditions."
    vars = sorted_union(variables(objective), [variables(c) for c in cons]...)
    @assert issubset(comm_gp, vars) "The commutative variables must be a subset of the variables."
    @assert !(is_unipotent && is_projective) "The problem cannot be both unipotent and projective."
    return PolyOpt(objective, cons, is_equality, vars, comm_gp, is_unipotent, is_projective)
end

function PolyOpt(objective::Polynomial{C,T}, constraints) where {C,T}
    return PolyOpt(objective, constraints, fill(false, length(constraints)), PolyVar{C}[], false, false)
end

# for user, unconstrained pop
function PolyOpt(objective::Polynomial{C,T}) where {C,T}
    return PolyOpt(objective, Polynomial{C,T}[], Bool[], PolyVar{C}[], false, false)
end

nvariables(pop::PolyOpt) = length(pop.variables)
nconstraints(pop::PolyOpt) = length(pop.constraints)
iscommutative(::PolyOpt{C}) where {C} = C
