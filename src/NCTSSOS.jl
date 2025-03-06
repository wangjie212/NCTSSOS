module NCTSSOS

using DynamicPolynomials
using DynamicPolynomials: AbstractVariable
using SparseArrays
const DP = DynamicPolynomials
using JuMP
using LinearAlgebra

export PolynomialOptimizationProblem
export MomentMethod, relax
export dualize
export remove_zero_degree, star, symmetric_canonicalize, get_basis
export support

include("pop.jl")

include("moment_solver.jl")

include("sos_solver.jl")

include("solver_utils.jl")

end
