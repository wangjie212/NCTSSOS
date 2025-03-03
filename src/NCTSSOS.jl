module NCTSSOS


using DynamicPolynomials
const DP = DynamicPolynomials


include("pop.jl")
export PolynomialOptimizationProblem

include("solver.jl")


include("solver_utils.jl")
export remove_zero_degree, star, symmetric_canonicalize, get_basis


end
