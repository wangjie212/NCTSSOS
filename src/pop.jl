abstract type OptimizationProblem end

# C: false if the variables are not commutative
# T: type of the coefficients
struct PolyOpt{C,T} <: OptimizationProblem
    objective::Polynomial{C,T}
    constraints::Vector{Polynomial{C,T}}
    is_equality::Vector{Bool}
    variables::Vector{PolyVar{C}}
end

# for user convenience
# be as general as possible
function PolyOpt(objective::Polynomial{C,T}, constraints, is_equality) where {C,T}
    @assert !(T isa Integer) "The polynomial coefficients can not be integers (not supported by JuMP solvers)."
    cons = collect(Polynomial{C,T}, constraints)
    @assert length(is_equality) == length(cons) "The number of constraints must be the same as the number of equality conditions."
    vars = sorted_union(variables(objective), [variables(c) for c in cons]...)
    return PolyOpt(objective, cons, is_equality, vars)
end

function PolyOpt(objective::Polynomial{C,T}, constraints) where {C,T}
    return PolyOpt(objective, constraints, fill(false, length(constraints)))
end

# for user, unconstrained pop
function PolyOpt(objective::Polynomial{C,T}) where {C,T}
    return PolyOpt(objective, Polynomial{C,T}[], Bool[])
end

nvariables(pop::PolyOpt) = length(pop.variables)
nconstraints(pop::PolyOpt) = length(pop.constraints)
iscommutative(::PolyOpt{C}) where {C} = C
