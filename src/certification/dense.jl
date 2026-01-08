using Arblib

function rational_certificate(f, ineq, eq, vars, r; partition=nothing, constraint=nothing, QUIET=false, TS = false, QUIETTS = true, tol=10e-10, eigprec=256)

    println("Building POP and computing numerical certificate...")
    clear_caches!()
    GC.gc()    

    if isempty(ineq) && isempty(eq)
        pop = [f]
    elseif isempty(ineq)
        pop = [f; eq]
    elseif isempty(eq)
        pop = [f; ineq]
    else
        pop = [f; ineq; eq]
    end

    time_sdp = @elapsed begin
        if partition == nothing && constraint == nothing
            opt,data = nctssos_first(pop, vars ,r, TS=TS,numeq=length(eq), QUIET=QUIETTS, Gram=true)
        elseif partition == nothing
            opt,data = nctssos_first(pop, vars ,r, TS=TS,numeq=length(eq), constraint=constraint, QUIET=QUIETTS, Gram=true)
        elseif constraint == nothing
            opt,data = nctssos_first(pop, vars ,r, TS=TS,numeq=length(eq), partition=partition, QUIET=QUIETTS, Gram=true)
        else
            opt,data = nctssos_first(pop, vars ,r, TS=TS,numeq=length(eq), partition=partition, constraint=constraint, QUIET=QUIETTS, Gram=true)
        end
    end

    println("Numerical certificate computed in $time_sdp seconds.")
    println("Computing LHS...")
    t_lhs = @elapsed begin
        g0 = data.GramMat[1][1]
        basis0 = basis_to_monovec(vars, data.basis[1])
        LHS = compute_LHS(f, opt, vars, r, data, ineq, eq, tol=tol)
    end
    println("LHS computed in $t_lhs seconds.")

    if partition == nothing && constraint == nothing
        bucket = nothing                     
    else
        arr_basis0   = arrays_of_basis(vars, basis0)
        _, bucket    = _build_pair_table_ideal(arr_basis0, vars, r;
                                               partition  = partition,
                                               constraint = constraint)
    end

    println("Rounding and projecting gram matrix...")

    t_rp = @elapsed begin
        if partition == nothing && constraint == nothing
            g_proj = round_project_gram(g0, basis0, LHS, r, vars; tol=tol)
            bucket = nothing
        else
            
            LHS_mon, LHS_coe = arrange(LHS, vars, partition=partition, constraint=constraint)
            LHS = sum(LHS_coe[i] * LHS_mon[i] for i in eachindex(LHS_coe))
            g_proj, bucket = round_project_gram_ideal(g0, basis0, LHS, r, vars; partition=partition, constraint=constraint, tol=tol)
        end
    end

    println("Gram matrix projected in $t_rp seconds.")

    println("Computing minimum eigenvalue of the projected gram matrix...")
    t_eig = @elapsed begin
        eig_min = rigorous_min_eig(Matrix{Rational{BigInt}}(g_proj); prec=eigprec)
    end
    println("Minimum eigenvalue computed in $t_eig seconds.")

    println("Computing right-hand side...")
    t_rhs = @elapsed begin
        RHS = bucket === nothing ? g0_to_decomp_cached(g_proj, basis0) : rhs_from_bucket_cached(bucket, g_proj)
    end
    
    println("Right-hand side computed in $t_rhs seconds.")
    println("Computing minimal eigenvalue of raw gram matrix...")
    eig_min_old = minimum(eigen(Matrix{Float64}(g0)).values)
    println("\nMinimum eigenvalue of the raw matrix")
    println(eig_min_old)

    diff_raw = NaN
    if bucket === nothing
        println("\nNo bucket (no ideal). Skipping raw LHS/RHS diagnostic.")
    else
        println("\nComputing raw LHS/RHS difference from bucket")
        t_raw = @elapsed begin
            diff_raw = raw_lhs_rhs_diff_from_bucket(bucket, g0, LHS)
        end
        println("Raw LHS/RHS difference computed in $t_raw seconds.")
        println("\nDifference between raw LHS and RHS:")
        println(diff_raw)
    end
    if !QUIET
        println("\nProjected matrix rational?")
        println(typeof(g_proj) == Matrix{Rational{BigInt}})

        diff = LHS - RHS
        coeff_map = Dict{String,Rational{BigInt}}()
        for (c,m) in zip(diff.a, diff.x)
            coeff_map[string(m)] = get(coeff_map, string(m), 0//1) + c
        end

        println("\nDifference between projected LHS and RHS:")
        println(isempty(coeff_map) ? 0.0 :
                Float64(sum(abs.(v) for v in values(coeff_map))))


        if(eig_min > 0)
            println("\nMatrix is already PSD, no need for further computation.")
            return opt
        else
            println("\nMatrix is not PSD, we need to compute a new bound.")

            println("\nMinimum eigenvalue of projected matrix:")
            println(eig_min)

            println("\nOld bound:")
            println(opt)

            println("\nNew bound after round and project procedure:")
            m = (-eig_min)*length(g_proj[1,:])
            new_bound = opt-m
            println(new_bound)

            println("\nDifference between bounds:")
            bound_diff = new_bound - opt
            println(bound_diff)
    
            return (newbound = new_bound, oldbound = opt, bdiff = bound_diff, eigproj = eig_min, eigraw =  eig_min_old, diffraw = diff_raw, glength = length(g_proj[1,:]))
        end
    end
    if(eig_min > 0)
        return opt
    else
        bound_diff = (-eig_min)*length(g_proj[1,:])
        new_bound = opt-bound_diff
        println("\n Bound lowered by $(bound_diff) to $new_bound.")
        return (newbound = new_bound, oldbound = opt, bdiff = bound_diff, eigproj = eig_min, eigraw = eig_min_old, diffraw = diff_raw, glength = length(g_proj[1,:]))
    end
end