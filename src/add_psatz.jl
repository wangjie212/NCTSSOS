mutable struct sohs_data
    cliquesize # size of cliques
    cliques # cliques
    basis # monomial basis
    blocksize # size of blocks
    blocks # block structrue
    tsupp # total support
    I # index sets of constraints
    gram # Gram variables
    constrs # constraint name
end

"""
    info = add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order; newton=false, obj="eigen", CS=false, TS="block", 
    SO=1, partition=0, comm_var=0, constraint=nothing, QUIET=false, constrs=nothing)

Add a Putinar's style SOHS representation of the nc polynomial `nonneg` to the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `nonneg`: a nonnegative polynomial constrained to be a Putinar's style SOHS on a nc semialgebraic set
- `vars`: the set of NCPOP variables
- `ineq_cons`: inequality constraints
- `eq_cons`: equality constraints
- `order`: relaxation order
- `CS`: method of chordal extension for correlative sparsity (`"MF"`, `"MD"`, `"NC"`, `false`)
- `TS`: type of term sparsity (`"block"`, `"MD"`, `"MF"`, `false`)
- `SO`: sparse order
- `partition`: the first 'partition' variables commutes with the remaining variables
- `comm_var`: the first 'comm_var' variables commutes each other
- `constraint`: nothing or "projection" or "unipotent"
- `QUIET`: run in the quiet mode (`true`, `false`)
- `constrs`: the constraint name used in the JuMP model

# Output arguments
- `info`: auxiliary data
"""

function add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order; newton=false, obj="eigen", partition=0, comm_var=0, constraint=nothing, CS=false, TS="block", SO=1, QUIET=false, constrs=nothing)
    cost = ncpoly(nonneg, vars)
    pop_cons = [ncpoly{Float64}([Int[]], [1]); [ncpoly(p, vars) for p in [ineq_cons; eq_cons]]]
    numeq = length(eq_cons)
    econs_type = nothing
    if numeq > 0
        econs_type = zeros(Int, numeq)
        for (i, p) in enumerate(eq_cons)
            temp = star(p)
            if temp == p
                econs_type[i] = 1
            elseif temp == -p
                econs_type[i] = -1
            end
        end
    end
    if obj == "trace"
        cost = cyclic_canon(cost)
    else
        cost = sym_canon(cost)
    end
    time = @elapsed begin
    cliques,cql,cliquesize = clique_decomp(cost, pop_cons, length(vars), alg=CS, QUIET=QUIET)
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
                basis[i][j] = get_ncbasis(cliquesize[i], order-ceil(Int, maxdeg(pop_cons[I[i][j]])/2), ind=cliques[i], binary=constraint!==nothing)
            end
            for ba in basis[i]
                ind = [_comm(item, partition, comm_var) == item for item in ba]
                ba = ba[ind]
            end
        else
            if obj == "trace"
                basis[1][1] = newton_cyclic(cost, length(vars), Int(maxdeg(cost)/2))
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
    blocks,cl,blocksize = get_pblocks(pop_cons, I, ksupp, basis, cliques, cql, TS=TS, SO=SO, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.([maximum.(bs) for bs in blocksize]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
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
    ksupp = reduce!.(ksupp, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    sort!(ksupp)
    unique!(ksupp)
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
    for i = 1:length(cost.supp)
        Locb = bfind(ksupp, cost.supp[i])
        if Locb === nothing
            @error "The monomial basis is not enough!"
            return info
        else
            cons[Locb] -= cost.coe[i]
        end
    end
    if constrs !== nothing
        @constraint(model, cons==zeros(length(ksupp)), base_name=constrs)
    else
        @constraint(model, cons .== 0)
    end
    info = sohs_data(cliquesize, cliques, basis, blocksize, blocks, ksupp, I, pos, constrs)
    return info
end

function get_pblocks(pop_cons, I, tsupp, basis, cliques, cql; TS="block", SO=1, obj="eigen", partition=0, comm_var=0, constraint=nothing)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    status = ones(Int, cql)
    for i = 1:cql
        ksupp = TS == false ? nothing : tsupp[[issubset(item, cliques[i]) for item in tsupp]]
        blocks[i],cl[i],blocksize[i],status[i] = get_pblocks(pop_cons[I[i]], ksupp, basis[i], TS=TS, SO=SO, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    end
    if minimum(status) == 1
        println("No higher TS step!")
    end
    return blocks,cl,blocksize
end

function get_pblocks(pop_cons, tsupp, basis; TS="block", SO=1, obj="eigen", partition=0, comm_var=0, constraint=nothing)
    blocks = Vector{Vector{Vector{Int}}}(undef, length(pop_cons))
    blocksize = Vector{Vector{Int}}(undef, length(pop_cons))
    cl = Vector{Int}(undef, length(pop_cons))
    status = 0
    if TS == false
        for k = 1:length(pop_cons)
            blocks[k],blocksize[k],cl[k] = [Vector(1:length(basis[k]))],[length(basis[k])],1
        end
    else
        for i = 1:SO
            if i > 1
                oblocksize = deepcopy(blocksize)
            end
            for k = 1:length(pop_cons)
                G = get_graph(tsupp, pop_cons[k].supp, basis[k], obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                if TS == "block"
                    blocks[k] = connected_components(G)
                    blocksize[k] = length.(blocks[k])
                    cl[k] = length(blocksize[k])
                else
                    blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS, minimize=false)
                end
            end
            if i > 1 && blocksize == oblocksize
                status = 1
                break
            end
            if i < SO
                tsupp = Vector{Int}[]
                for j = 1:length(blocks[1]), k = 1:blocksize[1][j], l = k:blocksize[1][j]
                    @inbounds bi = [basis[1][blocks[1][j][k]][end:-1:1]; basis[1][blocks[1][j][l]]]
                    push!(tsupp, bi)
                end
                tsupp = reduce!.(tsupp, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                sort!(tsupp)
                unique!(tsupp)
            end
        end
    end
    return blocks,cl,blocksize,status
end

function get_moment_matrix(moment, info; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    MomMat = Vector{Union{Float64, Symmetric{Float64}, Array{Float64,2}}}(undef, length(info.basis))
    for i = 1:length(info.basis)
        lb = length(info.basis[i][1])
        MomMat[i] = zeros(Float64, lb, lb)
        for j = 1:lb, k = j:lb
            bi = [info.basis[i][1][j][end:-1:1]; info.basis[i][1][k]]
            bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
            Locb = bfind(info.tsupp, bi)
            if Locb !== nothing
                MomMat[i][j,k] = moment[Locb]
            end
        end
        MomMat[i] = Symmetric(MomMat[i], :U)
    end
    return MomMat
end
