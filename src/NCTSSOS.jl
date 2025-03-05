module NCTSSOS


using DynamicPolynomials
using SparseArrays
const DP = DynamicPolynomials
using JuMP
using LinearAlgebra


include("pop.jl")
export PolynomialOptimizationProblem

include("moment_solver.jl")
export MomentMethod, init_moment_vector!, constrain_moment_matrix!, make_sdp, set_total_basis2var_dict!, get_total_basis2var_dict

include("sos_solver.jl")
export dualize, get_C_α_j


include("solver_utils.jl")
export remove_zero_degree, star, symmetric_canonicalize, get_basis
export support


end
