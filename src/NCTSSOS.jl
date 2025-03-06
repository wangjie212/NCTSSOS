module NCTSSOS

using DynamicPolynomials
using DynamicPolynomials: AbstractVariable, variables, coefficient, monomial, terms
using SparseArrays
using JuMP
using LinearAlgebra

export PolynomialOptimizationProblem
export moment_relax
export dualize

include("pop.jl")

include("moment_solver.jl")

include("sos_solver.jl")

include("solver_utils.jl")

end
