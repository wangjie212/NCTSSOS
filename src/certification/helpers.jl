using Arblib

function basis_to_monovec(xvars::Vector, data::Vector{Vector{UInt16}})
    key = (xvars, data)
    if haskey(_basis_cache, key)
        return _basis_cache[key]
    end
    basis = [ var_from_array(xvars, Int.(sub)) for sub in data ]
    _basis_cache[key] = basis
    return basis
end

normal_form(arr, vars, partition, constraint) =
    get!(_nf_cache, arr) do
        NCTSSOS.reduce(
            var_from_array(vars, arr),
            vars; obj = "eigen", partition = partition, constraint = constraint)
    end

concat_nf(u,v,vars,partition,constraint) =
    get!(_concat_nf_cache,
         (u, v,
          partition === nothing ? 0 : partition,
          String(constraint))) do
        normal_form(vcat(reverse(u),v), vars, partition, constraint)
    end

arrays_of_basis(vars, basis) = get!(_array_cache, (vars,basis)) do
    [array_from_var(b, vars) for b in basis]
end

function conj_poly(p::DynamicPolynomials.Polynomial)
    return p.a'*conj_monovec(Vector(p.x))
end

function symmetrize(p::DynamicPolynomials.Polynomial)
    return (p + conj_poly(p))/2
end

function rationalize_mat(mat; tol=1e-10)
    A = mat isa Matrix{Rational{BigInt}} ? mat : rationalize.(mat; tol=tol)
    return Rational{BigInt}.(A)
end

function rationalize_poly(poly::DynamicPolynomials.Polynomial; tol=1e-15)
    poly.a isa Vector{Rational{BigInt}} && return poly

    if eltype(poly.a) <: Real
         coeffs = Rational{BigInt}.(rationalize.(poly.a; tol=tol))
         return DynamicPolynomials.Polynomial(coeffs, poly.x)
    end

    return poly
end

function g0_to_decomp(mat, basis)
    return conj_monovec(basis)'*mat*basis
end

@inline function g0_to_decomp_cached(mat, basis)
    key = (pointer_from_objref(mat), objectid(basis))
    get!(_decomp_cache, key) do
        g0_to_decomp(mat, basis)    
    end
end

function gi_to_decomp_rat(vars, data, gram, i, ineq; tol=1e-15)
    basis      = basis_to_monovec(vars, data.basis[i])
    L          = rationalize_mat(psd_decomposition(gram), tol=tol)
    basis_conj = conj_monovec(basis)'
    ineq_rat   = rationalize_poly(ineq, tol=tol)

    decomp = zero(ineq_rat)
    @inbounds for j in 1:size(L,2)
        ℓ = view(L, :, j)
        v = ℓ' * basis
        v̄ = basis_conj * ℓ
        decomp += v̄ * ineq_rat * v
    end
    return decomp
end

function gi_to_decomp_rat_sparse(vars, basis_cell::Vector, gram, ineq; tol = 1e-15)
    basis      = basis_to_monovec(vars, basis_cell)
    L          = rationalize_mat(psd_decomposition(gram), tol = tol)
    basis_conj = conj_monovec(basis)'
    ineq_rat   = rationalize_poly(ineq, tol = tol)

    decomp = zero(ineq_rat)
    @inbounds for j in 1:size(L, 2)
        ℓ = view(L, :, j)
        decomp += (basis_conj * ℓ) * ineq_rat * (ℓ' * basis)
    end
    return decomp
end

@inline function rhs_from_bucket_cached(bucket, G)
    key = (pointer_from_objref(G), objectid(bucket))
    get!(_rhs_cache, key) do
        _rhs_from_bucket(bucket, G)  
    end
end

function gi_to_decomp_rat_eq_sparse(
        vars,      
        basis_cell::Vector,  
        gram,           
        eq;             
        tol = 1e-15)

    basis      = basis_to_monovec(vars, basis_cell)
    conjbasis  = conj_monovec(basis)'

    gram_rat   = rationalize_mat(gram; tol = tol)
    eq_rat     = rationalize_poly(eq;   tol = tol)

    coef1      = one(typeof(first(eq_rat.a)))            
    pbasis     = [ DynamicPolynomials.Polynomial([coef1], [m]) for m in basis ]
    pconj      = [ DynamicPolynomials.Polynomial([coef1], [m]) for m in vec(conjbasis) ]

    return pconj' * gram_rat * (eq_rat .* pbasis)
end

@inline function rationalize_mat_cached(G::AbstractMatrix; tol = 1e-10)
    get!(_gram_rat, pointer_from_objref(G)) do
        rationalize_mat(G; tol = tol)    
    end
end

@inline _pair_key(arr_basis, r::Int) = (Tuple(arr_basis), r)

function _coeff_matrix(vars, arr_basis, LHS)
    key = (vars, arr_basis, objectid(LHS))
    get!(_lhs_coeff_cache, key) do
        dim   = length(arr_basis)
        idx   = _poly_index(LHS)
        A     = zeros(Rational{BigInt}, dim, dim)
        for i in 1:dim, j in i:dim
            mkey = string(var_from_arrays(vars, arr_basis[i], arr_basis[j]))
            k    = get(idx, mkey, 0)
            c    = k==0 ? zero(LHS.a[1]) : LHS.a[k]
            A[i,j] = c
            A[j,i] = c
        end
        A
    end
end

function _coeff_matrix_mod(vars, arr_basis, LHS; partition=nothing, constraint=false)
    key = (Tuple(arr_basis), objectid(LHS),
       objectid(partition), objectid(constraint))
    get!(_lhs_mod_coeff_cache, key) do
        dim = length(arr_basis)
        A   = zeros(Rational{BigInt}, dim, dim)
        for i in 1:dim, j in i:dim
            A[i,j] = p_entry_mod(vars, arr_basis[i], arr_basis[j], LHS;
                                 partition=partition, constraint=constraint)
            A[j,i] = A[i,j]
        end
        A
    end
end

function gi_to_decomp_rat_eq(vars, data, gram, i, eq; tol=1e-15)
    basis   = basis_to_monovec(vars, data.basis[i])
    monoms  = conj_monovec(basis)

    eq_rat  = rationalize_poly(eq, tol=tol)
    gram_rat = rationalize_mat_cached(gram; tol = tol)

    coef1   = one(typeof(first(eq_rat.a)))

    pbasis = [ DynamicPolynomials.Polynomial([coef1], [m]) for m in basis   ]
    pconj  = [ DynamicPolynomials.Polynomial([coef1], [m]) for m in monoms   ]

    return pconj' * gram_rat * (eq_rat .* pbasis)
end


function psd_decomposition(mat)
    vals, vecs = eigen(Symmetric(mat))             
    valpos        = max.(vals, 0.0)                  
    return vecs * Diagonal(sqrt.(valpos))          
end

function make_psd(mat)                             
    vals, vecs = eigen(Symmetric(mat))
    valpos        = max.(vals, 0.0)
    return vecs * Diagonal(valpos) * vecs'           
end

function reverse_monomial(m::DynamicPolynomials.Monomial)
    word = Vector{eltype(m.vars)}()
    for (v, p) in zip(m.vars, m.z)
        for _ in 1:p
            push!(word, v)
        end
    end
    reverse!(word)

    new_vars = Vector{eltype(m.vars)}()
    new_pows = Int[]
    if !isempty(word)
        current_var = word[1]
        count = 1
        for i in 2:length(word)
            if word[i] == current_var
                count += 1
            else
                push!(new_vars, current_var)
                push!(new_pows, count)
                current_var = word[i]
                count = 1
            end
        end
        push!(new_vars, current_var)
        push!(new_pows, count)
    end

    return DynamicPolynomials.Monomial(new_vars, new_pows)
end

function conj_monovec(ms::Vector{<:DynamicPolynomials.Monomial})
    return map(reverse_monomial, ms)
end

function pairs_up_to_length(u, v, maxlen::Int, allowed::Set{Vector{Int}})
    X       = vcat(reverse(u), v)
    results = Tuple{Vector{Int},Vector{Int}}[]
    @inbounds for i in 0:length(X)
        vp = X[1:i]; wp = X[i+1:end]
        (length(vp) > maxlen || length(wp) > maxlen) && continue
        (vp in allowed && wp in allowed) && push!(results,(reverse(vp),wp))
    end
    return results, length(results)
end

function var_from_array(var, arr)
    result = DynamicPolynomials.Monomial(var, zeros(Int, length(var)))
    if(length(arr) == 0)
        return result
    end
    for i in 1:length(arr)
        result*=var[arr[i]]
    end
return result
end

function var_from_arrays(var, u, v)
    return var_from_array(var, reverse(u))*var_from_array(var, v)
end

_poly_index(p) = get!(_poly_index_cache, p) do
    Dict(string(m)=>i for (i,m) in enumerate(p.x))
end

function p_entry(var, u, v, p)
    idx = get(_poly_index(p), string(var_from_arrays(var,u,v)), 0)
    return idx==0 ? 0 : p.a[idx]
end

function p_entry_mod(var, u, v, p; partition=nothing, constraint=false)
    nf = concat_nf(u, v, var, partition, constraint)
    idx = get(_poly_index(p), string(nf), 0)
    return idx == 0 ? 0 : p.a[idx]
end

function array_from_var(M, var)
    if M == 1
        return []
    end
    return [findfirst(==(string(v)), string.(var)) for (v, e) in zip(M.vars, M.z) for _ in 1:e]
end

function compute_LHS(f, lambda, vars, r, data, ineqs, eqs; sparse=false, tol=1e-20)
    if sparse
        LHS = rationalize_poly(f, tol=tol) - rationalize(lambda; tol=tol)

        next_block = fill(2, length(data.cliques))

        for (j, g) in enumerate(ineqs)
            k   = poly_clique(g, data.cliques, vars)
            idx = next_block[k]                   
            next_block[k] += 1
            gram = data.GramMat[k][idx][1]
            if r == 1
                LHS -= rationalize(gram; tol=tol) * rationalize_poly(g, tol=tol)
            else
                basis_cell = data.basis[k][idx]
                LHS -= gi_to_decomp_rat_sparse(vars, basis_cell, gram, g; tol=tol)
            end
        end

        for (j, e) in enumerate(eqs)
            k   = poly_clique(e, data.cliques, vars)
            idx = next_block[k]
            next_block[k] += 1
            gram = data.GramMat[k][idx][1]
            if r == 1
                LHS -= rationalize(gram; tol=tol) * rationalize_poly(e, tol=tol)
            else
                basis_cell = data.basis[k][idx]
                LHS -= gi_to_decomp_rat_eq_sparse(vars, basis_cell, gram, e; tol=tol)
            end
        end

        return LHS
    end

    LHS = rationalize_poly(f, tol=tol) - rationalize(lambda; tol=tol)
    grammats = data.GramMat
    bases    = data.basis

    for i in 2:length(ineqs)+1
        if r == 1
            raw = grammats[i][1]
            gram = raw isa Matrix ? raw : raw[1]
            if gram > 0
                LHS -= rationalize(gram, tol=tol) * rationalize_poly(ineqs[i-1], tol=tol)
            end
        else
            basis = basis_to_monovec(vars, bases[i])
            gram  = grammats[i][1]
            decomp = gi_to_decomp_rat(vars, data, make_psd(gram), i, ineqs[i-1]; tol=tol)
            LHS -= decomp
        end
    end

    for i in length(ineqs)+2:length(grammats)
        if r == 1
            raw = grammats[i][1]
            gram = raw isa Matrix ? raw : raw[1]
            LHS -= rationalize(gram, tol=tol) * rationalize_poly(eqs[i-length(ineqs)-1], tol=tol)
        else
            basis = basis_to_monovec(vars, bases[i])
            gram  = grammats[i][1]
            decomp = gi_to_decomp_rat_eq(vars, data, gram, i, eqs[i-length(ineqs)-1]; tol=tol)
            LHS -= decomp
        end
    end

    return LHS
end

function project_gram(gram, basis, LHS, r, vars)
    dim        = length(basis)
    arr_basis  = arrays_of_basis(vars, basis)
    allowed    = Set(arr_basis)           

    key = _pair_key(arr_basis, r)
    pair_data = get!(_pair_table_cache, key) do       
        arr_to_idx = Dict(arr_basis[k] => k for k in 1:dim)
        tbl = Matrix{Tuple{Vector{Int},Vector{Int},Rational{BigInt}}}(undef, dim, dim)
        row_idx = Int[]; col_idx = Int[]
        for i in 1:dim, j in i:dim
            raw_pairs, n_pairs = pairs_up_to_length(arr_basis[i], arr_basis[j], r, allowed)
            inv_n = BigInt(1)//BigInt(n_pairs)
            empty!(row_idx); empty!(col_idx)
            for (up,vp) in raw_pairs
                push!(row_idx, arr_to_idx[up])
                push!(col_idx, arr_to_idx[vp])
            end
            tbl[i,j] = (copy(row_idx), copy(col_idx), inv_n)
            tbl[j,i] = tbl[i,j]
        end
        tbl
    end

    C = _coeff_matrix(vars, arr_basis, LHS)

    g_proj = zeros(Rational{BigInt}, dim, dim)
    for i in 1:dim, j in i:dim
        ui_vec, vj_vec, inv_n = pair_data[i,j]

        acc = gram[i,j] + inv_n * C[i,j]
        @inbounds for k in 1:length(ui_vec)
            acc -= inv_n * gram[ui_vec[k], vj_vec[k]]
        end
        g_proj[i,j] = acc
        g_proj[j,i] = acc
    end
    return g_proj
end

function round_project_gram(gram, basis, LHS, r, vars; tol=10e-20)
    gram_rat = rationalize_mat_cached(gram; tol = tol)
    g_proj = project_gram(gram_rat, basis, LHS, r, vars)
    return g_proj
end

function _build_pair_table_ideal(arr_basis, vars, r; partition=nothing, constraint=false)
    dim   = length(arr_basis)
    lens  = map(length, arr_basis)

    nf_map = Dict{Any, Vector{Tuple{Int,Int}}}()
    for u in 1:dim, v in 1:dim
        max(lens[u], lens[v]) > r && continue
        nf = concat_nf(arr_basis[u], arr_basis[v], vars, partition, constraint)
        push!(get!(nf_map, nf, Vector{Tuple{Int,Int}}()), (u,v))
    end

    tbl = Matrix{Tuple{Vector{Int},Vector{Int},Rational{BigInt}}}(undef, dim, dim)
    for i in 1:dim, j in i:dim
        pairs   = nf_map[ concat_nf(arr_basis[i], arr_basis[j], vars,
                                    partition, constraint) ]
        row_idx = [p[1] for p in pairs]
        col_idx = [p[2] for p in pairs]
        inv_n   = BigInt(1) // BigInt(length(pairs))
        tbl[i,j] = (row_idx, col_idx, inv_n)
        tbl[j,i] = tbl[i,j]
    end
    return tbl, nf_map
end

function _rhs_from_bucket(bucket::Dict{Any,Vector{Tuple{Int,Int}}}, G)
    it = collect(bucket)
    nf₀, pairs₀ = it[1]
    s₀ = sum(G[i,j] for (i,j) in pairs₀)
    term₀ = DynamicPolynomials.Polynomial([one(Rational{BigInt})], [nf₀])
    poly = s₀ * term₀
    for (nf, pairs) in it[2:end]
        s = sum(G[i,j] for (i,j) in pairs)
        term = DynamicPolynomials.Polynomial([one(Rational{BigInt})], [nf])
        poly += s * term
    end
    return poly
end


function project_gram_ideal(gram, basis, LHS, r, vars; partition=nothing, constraint=false)
    dim        = length(basis)
    arr_basis  = arrays_of_basis(vars, basis)

    key = (_pair_key(arr_basis,r)..., objectid(partition), objectid(constraint))
    pair_data, bucket = get!(_pair_table_ideal_cache, key) do
         _build_pair_table_ideal(arr_basis, vars, r;
                                 partition=partition, constraint=constraint)
     end

    C = _coeff_matrix_mod(vars, arr_basis, LHS;
                          partition=partition, constraint=constraint)

    g_proj = zeros(Rational{BigInt}, dim, dim)
    for i in 1:dim, j in i:dim
        ui_vec, vj_vec, inv_n = pair_data[i,j]

        acc = gram[i,j] + inv_n * C[i,j]
        @inbounds for k in 1:length(ui_vec)
            acc -= inv_n * gram[ui_vec[k], vj_vec[k]]
        end
        g_proj[i,j] = acc
        g_proj[j,i] = acc
    end
    return g_proj, bucket
end

function coeff_map_float(p::DynamicPolynomials.Polynomial)
    d = Dict{String,Float64}()
    @inbounds for (c, m) in zip(p.a, p.x)
        d[string(m)] = get(d, string(m), 0.0) + Float64(c)
    end
    return d
end

function raw_lhs_rhs_diff_from_bucket(bucket::Dict{Any,Vector{Tuple{Int,Int}}},
                                      G::AbstractMatrix{<:Real},
                                      LHS::DynamicPolynomials.Polynomial)
    coeffL = coeff_map_float(LHS) 

    err = 0.0
    @inbounds for (nf, pairs) in bucket
        s = 0.0
        for (i,j) in pairs
            s += Float64(G[i,j])  
        end
        cL = get(coeffL, string(nf), 0.0)
        err += abs(s - cL)
    end

    return err
end


function round_project_gram_ideal(gram, basis, LHS, r, vars; partition=nothing, constraint=false, tol=10e-10)
    gram_rat = rationalize_mat_cached(gram; tol = tol)
    g_proj, bucket = project_gram_ideal(gram_rat, basis, LHS, r, vars; partition=partition, constraint=constraint)
    return g_proj, bucket
end

function rigorous_min_eig(m::AbstractMatrix; prec::Int = 128)
    mc = AcbMatrix(m; prec=64)

    ev_approx, R_approx = Arblib.approx_eig_qr(mc; prec=prec)

    eps = Arblib.eig_global_enclosure(mc, ev_approx, R_approx; prec=prec)

    min_v = BigFloat(minimum(Arb.(real.(ev_approx); prec=prec)) - Arb(eps; prec=prec)) 
end

function poly_clique(p,
                     cliques::AbstractVector{<:AbstractVector},
                     vars::AbstractVector) :: Int

    vset = Set{Int}()
    for mon in p.x
        for v in mon.vars
            push!(vset, findfirst(==(v), vars))
        end
    end

    isempty(vset) && return 1          

    for k in eachindex(cliques)
        clique_set = Set(Int.(cliques[k]))
        all(i in clique_set for i in vset) && return k
    end
    error("Polynomial involves variables from several cliques.")
end

function merged_coeffs(p::DynamicPolynomials.Polynomial)
    d = Dict{String,Rational{BigInt}}()
    for (c,m) in zip(p.a, p.x)
        d[string(m)] = get(d, string(m), 0//1) + c
    end
    return d
end