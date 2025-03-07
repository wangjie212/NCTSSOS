module NCTSSOS

using DynamicPolynomials
using DynamicPolynomials:
    AbstractVariable, variables, coefficient, monomial, terms, isconstant
using SparseArrays, LinearAlgebra
using JuMP

export PolynomialOptimizationProblem
export moment_relax
export sos_dualize

include("pop.jl")

include("moment_solver.jl")

include("sos_solver.jl")

include("solver_utils.jl")

end
