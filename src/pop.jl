abstract type OptimizationProblem end

# C: false if the variables are not commutative
# T: type of the coefficients
struct PolynomialOptimizationProblem{C,T} <: OptimizationProblem
    objective::Polynomial{C,T}
    constraints::Vector{Polynomial{C,T}}
    variables::Vector{PolyVar{C}}
end

# for user convenience
# be as general as possible
function PolynomialOptimizationProblem(objective::Polynomial{C,T}, constraints) where {C,T}
    cons = collect(Polynomial{C,T}, constraints)
    variables = union(variables(objective), [variables(c) for c in cons]...)
    return PolynomialOptimizationProblem(objective, cons, variables)
end

nvariables(pop::PolynomialOptimizationProblem) = length(pop.variables)
nconstraints(pop::PolynomialOptimizationProblem) = length(pop.constraints)
iscommutative(::PolynomialOptimizationProblem{C}) where {C} = C
