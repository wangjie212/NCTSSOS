using Arblib

function rational_certificate_sparse(
    f,                    
    ineq, eq,           
    vars,                  
    r;                    
    partition  = nothing,  
    constraint = nothing,
    QUIET      = false,
    tol        = 1e-12,
    eigprec    = 256
)

    pop = [f]
    opt, data = ncpop(pop, vars, r;
        partition  = partition,
        constraint = constraint,
        TS         = false,
        numeq      = 0,
        QUIET      = true,
        Gram       = true)

    K       = length(data.cliques)
    lambda  = opt
    ideal   = !(partition === nothing && constraint === nothing)

    if !QUIET
        println("\nSparse numerical certificate computed.  Lower bound = ", lambda)
    end

    obj_poly = rationalize_poly(f; tol=tol) - rationalize(lambda; tol=tol)

    LHS_nf = obj_poly
    begin
        monL, coeL = arrange(LHS_nf, vars;
                             partition  = partition,
                             constraint = constraint)
        LHS_nf = sum(coeL[i]*monL[i] for i in eachindex(coeL))
        LHS_nf = rationalize_poly(LHS_nf; tol=tol)
    end

    buckets = Vector{Dict{Any,Vector{Tuple{Int,Int}}}}(undef, K)
    G0_rat  = Vector{Matrix{Rational{BigInt}}}(undef, K)
    basisKs = Vector{Vector}(undef, K)

    RHS_raw = Vector{DynamicPolynomials.Polynomial}(undef, K)

    for k in 1:K
        G0k = data.GramMat[k][1][1]
        G0_rat[k]  = rationalize_mat(G0k; tol=tol)
        basisKs[k] = basis_to_monovec(vars, data.basis[k][1])

        arr_basis = arrays_of_basis(vars, basisKs[k])
        _, bucket = _build_pair_table_ideal(arr_basis, vars, r;
                                            partition  = partition,
                                            constraint = constraint)
        buckets[k] = bucket

        RHS_raw[k] = rhs_from_bucket_cached(bucket, G0_rat[k])
    end

    for k in 1:K
        monR, coeR = arrange(RHS_raw[k], vars;
                             partition  = partition,
                             constraint = constraint)
        RHS_raw[k] = sum(coeR[i]*monR[i] for i in eachindex(coeR))
        RHS_raw[k] = rationalize_poly(RHS_raw[k]; tol=tol)
    end

    RHS_sum = Base.reduce(+, RHS_raw; init=zero(LHS_nf))
    diff    = merged_coeffs(LHS_nf - RHS_sum) 

    if !QUIET
        err_glob = sum(abs.(values(diff)))
        println("|(Sum raw RHS_k) − LHS| = ", Float64(err_glob))
    end

    G_proj = [copy(G0_rat[k]) for k in 1:K]

    r_map = diff 

    S_counts = Dict{String, Int}()

    for k in 1:K
        for (nf, pairs) in buckets[k]
            s = string(nf)
            S_counts[s] = get(S_counts, s, 0) + length(pairs)
        end
    end

    for (k, Gk) in enumerate(G_proj)
        bucket = buckets[k]
        for (nf, pairs) in bucket
            s = string(nf)
            r_t = get(r_map, s, 0//1)
            S_t = get(S_counts, s, 0)

            (r_t == 0//1 || S_t == 0) && continue
            lt = r_t / Rational{BigInt}(S_t) 

            @inbounds for (i,j) in pairs
                Gk[i,j] += lt
            end
        end
    end

    RHS_proj   = Vector{DynamicPolynomials.Polynomial}(undef, K)
    total_shift = zero(lambda)

    for k in 1:K
        RHS_k = rhs_from_bucket_cached(buckets[k], G_proj[k])

        monR, coeR = arrange(RHS_k, vars;
                             partition  = partition,
                             constraint = constraint)
        RHS_k = sum(coeR[i]*monR[i] for i in eachindex(coeR))
        RHS_k = rationalize_poly(RHS_k; tol=tol)
        RHS_proj[k] = RHS_k

        emin    = rigorous_min_eig(Matrix{Rational{BigInt}}(G_proj[k]); prec=eigprec)
        glength = size(G_proj[k], 1)
        total_shift += (-emin) * glength

        if !QUIET
            println("\nClique $k after global projection: min eig = $(emin)")
        end
    end

    new_bound = lambda - total_shift

    if !QUIET
        println("\nOld bound  = ", lambda)
        println("Total shift = ", total_shift)
        println("New bound   = ", new_bound)

        global_RHS  = Base.reduce(+, RHS_proj; init=zero(LHS_nf))
        monR, coeR  = arrange(global_RHS, vars;
                              partition  = partition,
                              constraint = constraint)
        global_RHS  = sum(coeR[i]*monR[i] for i in eachindex(coeR))
        global_RHS  = rationalize_poly(global_RHS; tol=tol)

        diff_coeff = merged_coeffs(global_RHS - LHS_nf)
        err_glob2  = sum(abs.(values(diff_coeff)))
    end

    return (newbound = new_bound,
            oldbound = lambda,
            totalshift = total_shift,
            RHS_per_clique = RHS_proj,
            Gproj_per_clique = G_proj,
            nCliques = K)
end