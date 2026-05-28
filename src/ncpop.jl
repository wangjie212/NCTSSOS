mutable struct ncpop_data
    pop # polynomial optimiztion problem
    obj # "eigen" or "trace"
    d # relaxation order
    n # number of variables
    numeq # number of equality constraints
    econs_type # type of equality constraints: 1 - symmetric, -1 - antisymmetric, 0 - others
    partition # the first 'partition' variables commute with the remaining variables
    comm_var # the first 'comm_var' variables commute each other
    constraint # nothing or "projection" or "unipotent"
    ksupp # extended support at the k-th TS step
    basis # monomial bases
    cliquesize # sizes of cliques
    cliques # cliques of variables
    I # index sets of inequality constraints
    blocksize # sizes of blocks
    blocks # block structure
    moment # Moment matrix
    GramMat # Gram matrices
end

"""
    opt,data = ncpop(pop, x, d; numeq=0, CS="MF", TS="block", solver="Mosek", writetofile=false, QUIET=false, 
    obj="eigen", partition=0, comm_var=0, constraint=nothing, solve=true, Gram=false)

Solve the `d`th level relaxation for noncommutative polynomial optimization problems. Return the optimum and other auxiliary data.
If `newton=true`, then compute a monomial basis by the Newton chip method.
If `reducebasis=true`, then remove monomials from the monomial basis by diagonal inconsistency.
If `TS="block"`, use maximal chordal extensions; if `TS="MD"`, use approximately smallest chordal extensions. 
If obj="eigen", minimize the eigenvalue; if obj="trace", minimize the trace.

# Arguments
- `pop`: vector of the objective function, inequality constraints, and equality constraints
- `x`: noncommuting variables
- `d`: relaxation order
- `numeq`: number of equality constraints
- `partition`: the first 'partition' variables commutes with the remaining variables
- `comm_var`: the first 'comm_var' variables commutes each other
- `constraint`: nothing or "projection" or "unipotent"
"""
function ncpop(pop::Vector{NCPoly{T}}, x, d; numeq=0, CS="MF", TS="block", newton=false, soc=false, normality=false, 
    QUIET=false, obj="eigen", solve=true, Gram=false, partition=0, comm_var=0, constraint=nothing, model=nothing, 
    mosek_setting=mosek_para(), writetofile=false) where {T<:Number}
    println("********************************** NCTSSOS **********************************")
    println("NCTSSOS is launching...")
    cost = ncpoly(pop[1], x)
    pop_cons = [ncpoly{T}([Int[]], [1]); [ncpoly(p, x) for p in pop[2:end]]]
    econs_type = nothing
    if numeq > 0
        econs_type = zeros(Int, numeq)
        for i = length(pop) - numeq + 1 : length(pop)
            temp = star(pop[i])
            if temp == pop[i]
                econs_type[i - (length(pop) - numeq)] = 1
            elseif temp == -pop[i]
                econs_type[i - (length(pop) - numeq)] = -1
            end
        end
    end
    if obj == "trace"
        cost = cyclic_canon(cost)
    else
        cost = sym_canon(cost)
    end
    time = @elapsed begin
    cliques,cql,cliquesize = clique_decomp(cost, pop_cons, length(x), alg=CS, QUIET=QUIET)
    end
    if CS != false && QUIET == false
        println("Obtained the variable cliques in $time seconds.\nThe maximal size of cliques is $(maximum(cliquesize)).")
    end
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    I = assign_constraint(pop_cons, cliques, cql)
    time = @elapsed begin
    basis = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    for i = 1:cql
        basis[i] = Vector{Vector{Vector{Int}}}(undef, length(I[i]))
        if newton == false
            for j = 1:length(I[i])
                basis[i][j] = get_ncbasis(cliquesize[i], d-ceil(Int, maxdeg(pop_cons[I[i][j]])/2), ind=cliques[i], binary=constraint!==nothing)
            end
            for ba in basis[i]
                ind = [_comm(item, partition, comm_var) == item for item in ba]
                ba = ba[ind]
            end
        else
            if obj == "trace"
                basis[1][1] = newton_cyclic(cost, length(x), Int(maxdeg(cost)/2))
            else
                basis[1][1] = newton_ncbasis(cost)
            end
       end
    end
    end
    if QUIET == false
        println("Obtained the monomial basis in $time seconds.")
    end
    ksupp = nothing
    if TS != false
        ksupp = reduce(vcat, [p.supp for p in [cost; pop_cons[2:end]]])
        for i = 1:cql
            if obj == "trace"
                append!(ksupp, [_cyclic_canon([item[end:-1:1]; item]) for item in basis[i][1]])
            else
                append!(ksupp, [[item[end:-1:1]; item] for item in basis[i][1]])
            end
        end
        ksupp = reduce!.(ksupp, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
        sort!(ksupp)
        unique!(ksupp)
    end
    time = @elapsed begin
    blocks,cl,blocksize = get_blocks(pop_cons, I, ksupp, basis, cliques, cql, TS=TS, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.([maximum.(bs) for bs in blocksize]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,GramMat = solvesdp(cost, pop_cons, basis, cql, I, blocks, cl, blocksize, numeq=numeq, econs_type=econs_type, TS=TS, soc=soc, 
    d=d, n=length(x), normality=normality, QUIET=QUIET, obj=obj, solve=solve, writetofile=writetofile, Gram=Gram, partition=partition, comm_var=comm_var, 
    constraint=constraint, model=model, mosek_setting=mosek_setting)
    data = ncpop_data([cost; pop_cons[2:end]], obj, d, length(x), numeq, econs_type, partition, comm_var, constraint, ksupp, basis, cliquesize, cliques, I, blocksize, blocks, 
    moment, GramMat)
    return opt,data
end

"""
    opt,data = ncpop(data; TS="block", QUIET=false, solve=true, Gram=false)

Compute higher TS steps. Return the optimum and other auxiliary data.
"""
function ncpop(data::ncpop_data; TS="block", soc=false, normality=false, QUIET=false, solve=true, Gram=false, model=nothing, mosek_setting=mosek_para(), writetofile=false)
    numeq = data.numeq
    econs_type = data.econs_type
    partition = data.partition
    comm_var = data.comm_var
    constraint = data.constraint
    obj = data.obj
    basis = data.basis
    cliques = data.cliques
    cql = length(cliques)
    pop_cons = [ncpoly{Float64}([Int[]], [1]); data.pop[2:end]] 
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize = get_blocks(pop_cons, data.I, data.ksupp, basis, cliques, cql, TS=TS, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    end
    if blocksize == data.blocksize
        println("No higher TS step!")
        opt = nothing
    else
        if QUIET == false
            mb = maximum(maximum.([maximum.(bs) for bs in blocksize]))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,GramMat = solvesdp(data.pop[1], pop_cons, basis, cql, data.I, blocks, cl, blocksize, numeq=numeq, econs_type=econs_type, 
        TS=TS, soc=soc, d=data.d, n=data.n, normality=normality, QUIET=QUIET, obj=obj, solve=solve, Gram=Gram, partition=partition, comm_var=comm_var, 
        constraint=constraint, writetofile=writetofile, model=model, mosek_setting=mosek_setting)
        data.moment = moment
        data.GramMat = GramMat
        data.ksupp = ksupp
        data.blocks = blocks
        data.blocksize = blocksize
    end
    return opt,data
end

function get_graph(tsupp::Vector{Vector{Int}}, supp, basis; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    lb = length(basis)
    G = SimpleGraph(lb)
    for i = 1:lb, j = i+1:lb
        r = 1
        while r <= length(supp)
            bi = [basis[i][end:-1:1]; supp[r]; basis[j]]
            bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
            if bfind(tsupp, bi) !== nothing
               break
            else
               r += 1
            end
        end
        if r <= length(supp)
           add_edge!(G, i, j)
        end
    end
    return G
end

function get_blocks(pop_cons::Vector{T}, tsupp, basis; TS="block", obj="eigen", QUIET=true, partition=0, comm_var=0, constraint=nothing) where {T<:ncpoly}
    blocks = Vector{Vector{Vector{Int}}}(undef, length(pop_cons))
    blocksize = Vector{Vector{Int}}(undef, length(pop_cons))
    cl = Vector{Int}(undef, length(pop_cons))
    if TS == false
        for k = 1:length(pop_cons)
            blocks[k],blocksize[k],cl[k] = [Vector(1:length(basis[k]))],[length(basis[k])],1
        end
    else
        for k = 1:length(pop_cons)
            G = get_graph(tsupp, pop_cons[k].supp, basis[k], obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
            if TS == "block"
                blocks[k] = connected_components(G)
                blocksize[k] = length.(blocks[k])
                cl[k] = length(blocksize[k])
            else
                blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS, minimize=false)
            end
            if QUIET == false
                sb = sort(Int.(unique(blocksize[k])), rev=true)
                numb = [sum(blocksize[k].== i) for i in sb]
                println("-----------------------------------------------------------------------------")
                println("The sizes of PSD blocks for the $k-th multiplier:\n$sb\n$numb")
                println("-----------------------------------------------------------------------------")
            end
        end
    end
    return blocks,cl,blocksize
end

function get_blocks(pop_cons, I, tsupp, basis, cliques, cql; TS="block", obj="eigen", partition=0, comm_var=0, constraint=nothing)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    for i = 1:cql
        ksupp = TS == false ? nothing : tsupp[[issubset(item, cliques[i]) for item in tsupp]]
        blocks[i],cl[i],blocksize[i] = get_blocks(pop_cons[I[i]], ksupp, basis[i], QUIET=true, TS=TS, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    end
    return blocks,cl,blocksize
end

function assign_constraint(pop_cons::Vector{T}, cliques, cql) where {T<:ncpoly}
    I = [Int[] for i = 1:cql]
    for (i, p) in enumerate(pop_cons)
        ind = findall(k->issubset(unique(reduce(vcat, p.supp)), cliques[k]), 1:cql)
        push!.(I[ind], i)
    end
    return I
end

function clique_decomp(cost::T, pop_cons::Vector{T}, n; alg="MF", QUIET=false) where {T<:ncpoly}
    if alg == false
        cliques,cql,cliquesize = [Vector(1:n)],1,[n]
    else
        G = SimpleGraph(n)
        foreach(item -> add_clique!(G, unique(item)), cost.supp)
        for p in pop_cons[2:end]
            add_clique!(G, unique(reduce(vcat, p.supp)))
        end
        if alg == "NC"
            cliques,cql,cliquesize = max_cliques(G)
        else
            cliques,cql,cliquesize = chordal_cliques!(G, method=alg, minimize=true)
        end
    end
    if QUIET == false
        uc = unique(cliquesize)
        sizes = [sum(cliquesize.== i) for i in uc]
        println("-----------------------------------------------------------------------------")
        println("The clique sizes of varibles:\n$uc\n$sizes")
        println("-----------------------------------------------------------------------------")
    end
    return cliques,cql,cliquesize
end

function solvesdp(cost, pop_cons, basis, cql, I, blocks, cl, blocksize; numeq=0, econs_type=nothing, TS="block", QUIET=false, obj="eigen", solve=true, 
    Gram=false, partition=0, comm_var=0, constraint=nothing, soc=false, d=0, n=0, normality=false, model=nothing, mosek_setting=mosek_para(), writetofile=false)
    ksupp = Vector{Int}[]
    for i = 1:cql, j = 1:cl[i][1], k = 1:blocksize[i][1][j], r = k:blocksize[i][1][j]
        @inbounds bi = [basis[i][1][blocks[i][1][j][k]][end:-1:1]; basis[i][1][blocks[i][1][j][r]]]
        push!(ksupp, bi)
    end
    if TS != false
        for i = 1:cql, (j, w) in enumerate(I[i][2:end]), l = 1:cl[i][j+1], t = 1:blocksize[i][j+1][l], r = t:blocksize[i][j+1][l], item in pop_cons[w].supp
            @inbounds bi = [basis[i][j+1][blocks[i][j+1][l][t]][end:-1:1]; item; basis[i][j+1][blocks[i][j+1][l][r]]]
            push!(ksupp, bi)
        end
    end
    if soc == true && cql == 1
        slb = basis[1][1][length.(basis[1][1]) .<= 2d-maxdeg(cost)]
        spb = basis[1][1][length.(basis[1][1]) .<= d-Int(ceil(maxdeg(cost)/2))]
        lspb = length(spb)
        if TS != false
            for item in slb, nitem in cost.supp
                push!(ksupp, [nitem; item], [item; nitem])
            end
            for i = 1:lspb, j = i:lspb, item in cost.supp
                push!(ksupp, [spb[i][end:-1:1]; item; spb[j]], [item; spb[i][end:-1:1]; spb[j]], [spb[i][end:-1:1]; spb[j]; item])
            end
        end
    end
    if normality == true && cql == 1
        wbasis = basis[1][1]
        wbs = length(wbasis)
        for i = 1:n, j = 1:wbs, k = j:wbs
            push!(ksupp, [wbasis[j][end:-1:1]; i; wbasis[k]], [wbasis[j][end:-1:1]; i; i; wbasis[k]])
        end
    end
    ksupp = reduce!.(ksupp, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    sort!(ksupp)
    unique!(ksupp)
    if QUIET == false
        println("There are $(length(ksupp)) affine constraints.")
    end
    objv = moment = GramMat = nothing
    if solve == true
        if QUIET == false
            println("Assembling the SDP...")
        end
        if model === nothing
            model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = [AffExpr(0) for i = 1:length(ksupp)]
        pos = Vector{Vector{Vector{Union{VariableRef,Vector{VariableRef},Symmetric{VariableRef},Matrix{VariableRef}}}}}(undef, cql)
        for i = 1:cql
            pos[i] = Vector{Vector{Union{VariableRef,Vector{VariableRef},Symmetric{VariableRef},Matrix{VariableRef}}}}(undef, length(I[i]))
            for (j, w) in enumerate(I[i])
                pos[i][j] = Vector{Union{VariableRef,Vector{VariableRef},Symmetric{VariableRef},Matrix{VariableRef}}}(undef, cl[i][j])
                for l = 1:cl[i][j]
                    bs = blocksize[i][j][l]
                    if w <= length(pop_cons) - numeq || econs_type[w - (length(pop_cons) - numeq)] == 1
                        if w <= length(pop_cons) - numeq
                            pos[i][j][l] = @variable(model, [1:bs, 1:bs], PSD)
                        else
                            pos[i][j][l] = @variable(model, [1:bs, 1:bs], Symmetric)
                        end
                        for t = 1:bs, r = t:bs, (s, item) in enumerate(pop_cons[w].supp)
                            @inbounds bi = [basis[i][j][blocks[i][j][l][t]][end:-1:1]; item; basis[i][j][blocks[i][j][l][r]]]
                            bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                            Locb = bfind(ksupp, bi)
                            if t == r
                                @inbounds add_to_expression!(cons[Locb], pop_cons[w].coe[s], pos[i][j][l][t,r])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*pop_cons[w].coe[s], pos[i][j][l][t,r])
                            end
                        end
                    elseif econs_type[w - (length(pop_cons) - numeq)] == -1
                        pos[i][j][l] = @variable(model, [1:Int(bs*(bs-1)/2)])
                        for t = 1:bs-1, r = t+1:bs, (s, item) in enumerate(pop_cons[w].supp)
                            @inbounds bi = [basis[i][j][blocks[i][j][l][t]][end:-1:1]; item; basis[i][j][blocks[i][j][l][r]]]
                            bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                            Locb = bfind(ksupp, bi)
                            @inbounds add_to_expression!(cons[Locb], 2*pop_cons[w].coe[s], pos[i][j][l][Int((2bs-t)*(t-1)/2)+r-t])
                        end
                    else
                        pos[i][j][l] = @variable(model, [1:bs, 1:bs])
                        for t = 1:bs, r = 1:bs, (s, item) in enumerate(pop_cons[w].supp)
                            @inbounds bi = [basis[i][j][blocks[i][j][l][t]][end:-1:1]; item; basis[i][j][blocks[i][j][l][r]]]
                            bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                            Locb = bfind(ksupp, bi)
                            @inbounds add_to_expression!(cons[Locb], 2*pop_cons[w].coe[s], pos[i][j][l][t,r])
                        end
                    end
                end
            end
        end
        if soc == true && cql == 1
            fr = @variable(model, [1:length(slb)])
            for (i, item) in enumerate(slb), (s, nitem) in enumerate(cost.supp)
                bi = reduce!([nitem; item], obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                Locb = bfind(ksupp, bi)
                add_to_expression!(cons[Locb], cost.coe[s], fr[i])
                bi = reduce!([item; nitem], obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                Locb = bfind(ksupp, bi)
                add_to_expression!(cons[Locb], -cost.coe[s], fr[i])
            end
            npos = @variable(model, [1:lspb, 1:lspb], PSD)
            for i = 1:lspb, j = i:lspb, (s, item) in enumerate(cost.supp)
                @inbounds bi = [spb[i][end:-1:1]; item; spb[j]]
                bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                Locb = bfind(ksupp, bi)
                if i == j
                    @inbounds add_to_expression!(cons[Locb], cost.coe[s], npos[i,j])
                else
                    @inbounds add_to_expression!(cons[Locb], 2*cost.coe[s], npos[i,j])
                end
                @inbounds bi = [item; spb[i][end:-1:1]; spb[j]]
                bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                Locb = bfind(ksupp, bi)
                if i == j
                    @inbounds add_to_expression!(cons[Locb], -0.5*cost.coe[s], npos[i,j])
                else
                    @inbounds add_to_expression!(cons[Locb], -cost.coe[s], npos[i,j])
                end
                @inbounds bi = [spb[i][end:-1:1]; spb[j]; item]
                bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                if i == j
                    @inbounds add_to_expression!(cons[Locb], -0.5*cost.coe[s], npos[i,j])
                else
                    @inbounds add_to_expression!(cons[Locb], -cost.coe[s], npos[i,j])
                end
            end
        end
        if normality == true && cql == 1
            for i = 1:n
                hnom = @variable(model, [1:2wbs, 1:2wbs], PSD)
                for j = 1:wbs, k = j:wbs
                    @inbounds bi = [wbasis[j][end:-1:1]; wbasis[k]]
                    bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                    Locb = bfind(ksupp, bi)
                    if j == k
                        @inbounds add_to_expression!(cons[Locb], hnom[j,k])
                    else
                        @inbounds add_to_expression!(cons[Locb], 2, hnom[j,k])
                    end
                    @inbounds bi = [wbasis[j][end:-1:1]; i; i; wbasis[k]]
                    bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                    Locb = bfind(ksupp, bi)
                    if j == k
                        @inbounds add_to_expression!(cons[Locb], hnom[j+wbs,k+wbs])
                    else
                        @inbounds add_to_expression!(cons[Locb], 2, hnom[j+wbs,k+wbs])
                    end    
                    @inbounds bi = [wbasis[j][end:-1:1]; i; wbasis[k]]
                    bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                    Locb = bfind(ksupp, bi)
                    if j == k
                        @inbounds add_to_expression!(cons[Locb], 2, hnom[j,k+wbs])
                    else
                        @inbounds add_to_expression!(cons[Locb], 2, hnom[j,k+wbs]+hnom[k,j+wbs])
                    end     
                end
            end
        end
        for i = 1:length(cost.supp)
            Locb = bfind(ksupp, cost.supp[i])
            if Locb === nothing
                @error "The monomial basis is not enough!"
                return nothing,ksupp,nothing,nothing
            else
                cons[Locb] -= cost.coe[i]
            end
        end
        @variable(model, lower)
        cons[1] += lower
        @constraint(model, con, cons==zeros(length(ksupp)))
        @objective(model, Max, lower)
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
        end
        if writetofile != false
            write_to_file(dualize(model), writetofile)
        else
            if QUIET == false
                println("Solving the SDP...")
            end
            time = @elapsed begin
            optimize!(model)
            end
            if QUIET == false
                println("SDP solving time: $time seconds.")
            end
            status = termination_status(model)
            objv = objective_value(model)
            if status != MOI.OPTIMAL
                println("termination status: $status")
                status = primal_status(model)
                println("solution status: $status")
            end
            println("optimum = $objv")
            if Gram == true
                GramMat = Vector{Vector{Vector{Union{Float64,Matrix{Float64}}}}}(undef, cql)
                for i = 1:cql
                    GramMat[i] = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, length(I[i]))
                    for (j, w) in enumerate(I[i])
                        if w > length(pop_cons) - numeq && econs_type[w - (length(pop_cons) - numeq)] == -1
                            GramMat[i][j] = Vector{Union{Float64,Matrix{Float64}}}(undef, cl[i][j])
                            for l = 1:cl[i][j]
                                GramMat[i][j][l] = zeros(blocksize[i][j][l], blocksize[i][j][l])
                                for t = 1:blocksize[i][j][l]-1, r = t+1:blocksize[i][j][l]
                                    GramMat[i][j][l][t,r] = value(pos[i][j][l][Int((2*blocksize[i][j][l]-t)*(t-1)/2)+r-t])
                                    GramMat[i][j][l][r,t] = -GramMat[i][j][l][t,r]
                                end
                            end
                        else
                            GramMat[i][j] = [value.(pos[i][j][l]) for l = 1:cl[i][j]]
                        end
                    end
                end
            end
            dual_var = -dual(con)
            moment = Vector{Vector{Matrix{Float64}}}(undef, cql)
            for i = 1:cql
                moment[i] = Vector{Matrix{Float64}}(undef, cl[i][1])
                for k = 1:cl[i][1]
                    moment[i][k] = zeros(blocksize[i][1][k], blocksize[i][1][k])
                    for j = 1:blocksize[i][1][k], r = j:blocksize[i][1][k]
                        @inbounds bi = [basis[i][1][blocks[i][1][k][j]][end:-1:1]; basis[i][1][blocks[i][1][k][r]]]
                        bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                        Locb = bfind(ksupp, bi)
                        moment[i][k][j,r] = dual_var[Locb]
                    end
                    moment[i][k] = Symmetric(moment[i][k], :U)
                end
            end
        end
    end
    return objv,ksupp,moment,GramMat
end
