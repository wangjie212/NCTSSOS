module NCTSSOS

using DynamicPolynomials
using MultivariatePolynomials
using JuMP
using MosekTools
using Graphs
using ChordalGraph
using LinearAlgebra
using SparseArrays
import DynamicPolynomials as DP
import MultivariatePolynomials as MP

export ncpop, statepop, add_psatz!, get_moment_matrix, add_SOHS!, add_poly!, arrange
export star, symmetrize, rational_certificate, rational_certificate_sparse

mutable struct mosek_para
    tol_pfeas::Float64
    tol_dfeas::Float64
    tol_relgap::Float64
    time_limit::Int64
    num_threads::Int64
end

mosek_para() = mosek_para(1e-8, 1e-8, 1e-8, -1, 0)

include("polynomial.jl")
include("utils.jl")
include("ncpop.jl")
include("add_psatz.jl")
include("statepop.jl")

# --- Certification routines  ---
include(joinpath(@__DIR__, "certification", "caches.jl"))
include(joinpath(@__DIR__, "certification", "helpers.jl"))
include(joinpath(@__DIR__, "certification", "dense.jl"))
include(joinpath(@__DIR__, "certification", "sparse.jl"))

end
