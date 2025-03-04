abstract type OptimizationProblem end

struct PolynomialOptimizationProblem{C,N,VV<:AbstractVector{<:DP.AbstractVariable},PD<:AbstractPolynomialLike} <: OptimizationProblem
    objective::PD
    constraints::NTuple{N,PD}
    variables::VV
end

function PolynomialOptimizationProblem(objective::PD, constraints::AbstractVector{PD}, variables::VV) where {VV<:AbstractVector{<:DP.AbstractVariable},PD<:AbstractPolynomialLike}
    N = length(constraints)
	C = any(x -> !(DP.iscomm(typeof(x))), variables)
    PolynomialOptimizationProblem{C,N,VV,PD}(objective, ntuple(i -> constraints[i], N), variables)
end

PolynomialOptimizationProblem(objective::PD, variables::VV) where {VV<:AbstractVector{<:DP.AbstractVariable},PD<:AbstractPolynomialLike} = PolynomialOptimizationProblem(objective, PD[], variables)

get_objective(pop::PolynomialOptimizationProblem) = pop.objective

get_constraints(pop::PolynomialOptimizationProblem) = pop.constraints

get_variables(pop::PolynomialOptimizationProblem) = pop.variables

nvariables(pop::PolynomialOptimizationProblem) = length(pop.variables)

nconstraints(::PolynomialOptimizationProblem{C,N,VV,PD}) where {C,N,VV,PD} = N

iscommutative(::PolynomialOptimizationProblem{C,N,VV,PD}) where {C,N,VV,PD} = C
