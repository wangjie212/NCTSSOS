global _basis_cache     = IdDict{Tuple{Vector,Vector{Vector{UInt16}}},Vector}()
global _poly_index_cache = IdDict{Any,Dict{String,Int}}()
global _array_cache      = IdDict{Tuple{Vector,Vector},Vector{Vector{Int}}}()

global _decomp_cache = IdDict{Tuple{Ptr{Nothing},UInt}, DynamicPolynomials.Polynomial}()
global _rhs_cache    = IdDict{Tuple{Ptr{Nothing},UInt}, DynamicPolynomials.Polynomial}()
global _nf_cache        = IdDict{Vector{Int},Any}()
global _concat_nf_cache =
    IdDict{Tuple{Vector{Int},Vector{Int},Int,String},Any}()

global _gram_rat =
    IdDict{Ptr{Nothing}, Matrix{Rational{BigInt}}}()

global _pair_table_cache        = IdDict{Tuple,Matrix{Tuple{Vector{Int},Vector{Int},Rational{BigInt}}}}()
global _pair_table_ideal_cache = IdDict{Tuple,
     Tuple{
         Matrix{Tuple{Vector{Int},Vector{Int},Rational{BigInt}}},
         Dict{Any,Vector{Tuple{Int,Int}}}
     }
 }()

global _lhs_coeff_cache      = IdDict{Tuple, Matrix{Rational{BigInt}}}()
global _lhs_mod_coeff_cache  = IdDict{Tuple, Matrix{Rational{BigInt}}}()
    
function clear_caches!()
    global _basis_cache      = Dict{Any,Any}()
    global _poly_index_cache = Dict{Any,Any}()
    global _array_cache      = Dict{Any,Any}()
    global _decomp_cache     = Dict{Any,Any}()   # if you keep it
    global _rhs_cache        = Dict{Any,Any}()   # if you keep it
    global _nf_cache         = Dict{Any,Any}()
    global _concat_nf_cache  = Dict{Any,Any}()
    global _gram_rat         = Dict{Any,Any}()
    global _pair_table_cache       = Dict{Any,Any}()
    global _pair_table_ideal_cache = Dict{Any,Any}()
    global _lhs_coeff_cache        = Dict{Any,Any}()
    global _lhs_mod_coeff_cache    = Dict{Any,Any}()
end