abstract type OptimizationProblem end

struct PolynomialOptimizationProblem{C,N,VV<:AbstractVector{<:AbstractVariable},PD<:AbstractPolynomialLike{<:Real}} <: OptimizationProblem
    objective::PD
    constraints::NTuple{N,PD}
    variables::VV
end

function PolynomialOptimizationProblem(objective::PD, constraints::AbstractVector{PD}, variables::VV) where {VV<:AbstractVector{<:AbstractVariable},PD<:AbstractPolynomialLike{<:Real}}
    N = length(constraints)
	C = any(x -> !(DP.iscomm(typeof(x))), variables)
    PolynomialOptimizationProblem{C,N,VV,PD}(objective, ntuple(i -> constraints[i], N), variables)
end

get_objective(pop::PolynomialOptimizationProblem) = pop.objective

nvariables(pop::PolynomialOptimizationProblem) = length(pop.variables)

nconstraints(::PolynomialOptimizationProblem{C,N,VV,PD}) where {C,N,VV,PD} = N

iscommutative(::PolynomialOptimizationProblem{C,N,VV,PD}) where {C,N,VV,PD} = C
