module NCTSSOS

using DynamicPolynomials
using DynamicPolynomials: AbstractVariable, variables, coefficient, monomial, terms
using SparseArrays
using JuMP
using LinearAlgebra

export PolynomialOptimizationProblem
export MomentProblem, relax
export dualize
export remove_zero_degree, star, symmetric_canonicalize, get_basis
export support

include("pop.jl")

include("moment_solver.jl")

include("sos_solver.jl")

include("solver_utils.jl")

end
