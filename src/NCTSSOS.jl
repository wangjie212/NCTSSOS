module NCTSSOS

using DynamicPolynomials
using MultivariatePolynomials
using JuMP
using MosekTools
using Graphs
using ChordalGraph
using MetaGraphs
using LinearAlgebra
using SparseArrays
using COSMO
import DynamicPolynomials as DP

export nctssos_first, nctssos_higher!, cs_nctssos_first, cs_nctssos_higher!, ptraceopt_first, ptraceopt_higher!, 
pstateopt_first, pstateopt_higher!, Werner_witness_first, Werner_witness_higher!, cosmo_para, mixword, traceopt_first,
traceopt_higher!, cpstateopt_first, stateopt_first, stateopt_higher!, cpstateopt_higher!, add_psatz!, get_moment_matrix,
add_SOHS!, add_poly!, arrange, star, symmetrize, rational_certificate, rational_certificate_sparse

const Mono = DP.Monomial{DP.NonCommutative{DP.CreationOrder}, Graded{LexOrder}}
const Poly{T} = DP.Polynomial{DP.NonCommutative{DP.CreationOrder}, Graded{LexOrder}, T}

include("clique_merge.jl")
include("ncupop.jl")
include("nccpop.jl")
include("mixncpop.jl")
include("trace.jl")
include("state.jl")
include("complex.jl")
include("add_psatz.jl")
include("utils.jl")

# --- Certification routines  ---
include(joinpath(@__DIR__, "certification", "caches.jl"))
include(joinpath(@__DIR__, "certification", "helpers.jl"))
include(joinpath(@__DIR__, "certification", "dense.jl"))
include(joinpath(@__DIR__, "certification", "sparse.jl"))

end
