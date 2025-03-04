module NCTSSOS


using DynamicPolynomials
const DP = DynamicPolynomials
using JuMP


include("pop.jl")
export PolynomialOptimizationProblem

include("solver.jl")
export MomentMethod, SDSOSMethod, init_moment_vector!, constrain_moment_matrix!, make_sdp


include("solver_utils.jl")
export remove_zero_degree, star, symmetric_canonicalize, get_basis
export support


end
