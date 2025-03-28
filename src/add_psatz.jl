mutable struct struct_data
    cql # number of cliques
    cliquesize # size of cliques
    cliques # cliques
    basis # monomial basis
    cl # number of blocks
    blocksize # size of blocks
    blocks # blocks
    tsupp # total support
    I # index sets of constraints
    gram # Gram variables
    constrs # constraint name
end

"""
    info = add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order; obj="eigen", CS=false, cliques=[], TS="block", 
    SO=1, partition=0, constraint=nothing, QUIET=false, constrs=nothing)

Add a Putinar's style SOHS representation of the nc polynomial `nonneg` to the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `nonneg`: a nonnegative polynomial constrained to be a Putinar's style SOHS on a nc semialgebraic set
- `vars`: the set of NCPOP variables
- `ineq_cons`: inequality constraints
- `eq_cons`: equality constraints
- `order`: relaxation order
- `CS`: method of chordal extension for correlative sparsity (`"MF"`, `"MD"`, `"NC"`, `false`)
- `cliques`: the set of cliques used in correlative sparsity
- `TS`: type of term sparsity (`"block"`, `"signsymmetry"`, `"MD"`, `"MF"`, `false`)
- `SO`: sparse order
- `QUIET`: run in the quiet mode (`true`, `false`)
- `constrs`: the constraint name used in the JuMP model

# Output arguments
- `model`: the modified JuMP model
- `info`: other auxiliary data
"""

function add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order; obj="eigen", partition=0, constraint=nothing, CS=false, cliques=[], TS="block", SO=1, QUIET=false, constrs=nothing)
    n = length(vars)
    _,fsupp,fcoe = poly_info(nonneg, vars)
    m = length(ineq_cons) + length(eq_cons)
    numeq = length(eq_cons)
    if numeq > 0
        eq_constraint_type = zeros(Int, numeq)
        for i = 1:numeq
            temp = star(eq_cons[i])
            if temp == eq_cons[i]
                eq_constraint_type[i] = 1
            elseif temp == -eq_cons[i]
                eq_constraint_type[i] = -1
            end
        end
    end
    if m > 0
        _,gsupp,gcoe = polys_info([ineq_cons; eq_cons], vars)
        dg = [maximum(length.(item)) for item in gsupp]
    else
        gsupp = Vector{Vector{UInt16}}[]
    end
    if obj == "trace"
        fsupp,fcoe = cyclic_canon(fsupp, fcoe, type=AffExpr)
    else
        fsupp,fcoe = sym_canon(fsupp, fcoe, type=AffExpr)
    end
    if CS != false
        if cliques == []
            cliques,cql,cliquesize = clique_decomp(n, m, fsupp, gsupp, alg=CS, QUIET=false)
        else
            cql = length(cliques)
            cliquesize = length.(cliques)
        end
    else
        cliques,cql,cliquesize = [Vector(1:n)],1,[n]
    end
    I = assign_constraint(m, gsupp, cliques, cql)
    basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
    for t = 1:cql
        basis[t] = Vector{Vector{Vector{UInt16}}}(undef, length(I[t])+1)
        basis[t][1] = get_ncbasis(cliquesize[t], order, ind=cliques[t], binary=constraint!==nothing)
        for s = 1:length(I[t])
            basis[t][s+1] = get_ncbasis(cliquesize[t], order-ceil(Int, dg[I[t][s]]/2), ind=cliques[t], binary=constraint!==nothing)
        end
        if partition > 0
            for ba in basis[t]
                ind = [_comm(item, partition) == item for item in ba]
                ba = ba[ind]
            end
        end
    end
    blocks,cl,blocksize = get_blocks(I, m, fsupp, gsupp, basis, cliques, cql, TS=TS, SO=SO, QUIET=QUIET, obj=obj, partition=partition, constraint=constraint)
    tsupp = Vector{UInt16}[]
    for i = 1:cql
        for j = 1:cl[i][1], k = 1:blocksize[i][1][j], r = k:blocksize[i][1][j]
            @inbounds bi = [basis[i][1][blocks[i][1][j][k]][end:-1:1]; basis[i][1][blocks[i][1][j][r]]]
            push!(tsupp, bi)
        end
        if TS != false
            for (j, w) in enumerate(I[i])
                for l = 1:cl[i][j+1], t = 1:blocksize[i][j+1][l], r = t:blocksize[i][j+1][l], item in gsupp[w]
                    ind1 = blocks[i][j+1][l][t]
                    ind2 = blocks[i][j+1][l][r]
                    @inbounds bi = [basis[i][j+1][ind1][end:-1:1]; item; basis[i][j+1][ind2]]
                    push!(tsupp, bi)
                end
            end
        end
    end
    tsupp = reduce!.(tsupp, obj=obj, partition=partition, constraint=constraint)
    sort!(tsupp)
    unique!(tsupp)
    ltsupp = length(tsupp)
    cons = [AffExpr(0) for i=1:ltsupp]
    pos = Vector{Vector{Vector{Union{VariableRef,Vector{VariableRef},Symmetric{VariableRef},Matrix{VariableRef}}}}}(undef, cql)
    for t = 1:cql
        pos[t] = Vector{Vector{Union{VariableRef,Vector{VariableRef},Symmetric{VariableRef},Matrix{VariableRef}}}}(undef, 1+length(I[t]))
        pos[t][1] = Vector{Union{VariableRef,Vector{VariableRef},Symmetric{VariableRef},Matrix{VariableRef}}}(undef, cl[t][1])
        for i = 1:cl[t][1]
            bs = blocksize[t][1][i]
            if bs == 1
                pos[t][1][i] = @variable(model, lower_bound=0)
                @inbounds bi = [basis[t][1][blocks[t][1][i][1]][end:-1:1]; basis[t][1][blocks[t][1][i][1]]]
                bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint)
                Locb = bfind(tsupp, ltsupp, bi)
                @inbounds add_to_expression!(cons[Locb], pos[t][1][i])
            else
               pos[t][1][i] = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:bs, r = j:bs
                    @inbounds bi = [basis[t][1][blocks[t][1][i][j]][end:-1:1]; basis[t][1][blocks[t][1][i][r]]]
                    bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint)
                    Locb = bfind(tsupp, ltsupp, bi)
                    if j == r
                       @inbounds add_to_expression!(cons[Locb], pos[t][1][i][j,r])
                    else
                       @inbounds add_to_expression!(cons[Locb], 2, pos[t][1][i][j,r])
                    end
               end
            end
        end
        for k = 1:length(I[t])
            pos[t][k+1] = Vector{Union{VariableRef,Vector{VariableRef},Symmetric{VariableRef},Matrix{VariableRef}}}(undef, cl[t][k+1])
            for i = 1:cl[t][k+1]
                bs = blocksize[t][k+1][i]
                if bs == 1
                    if k <= m-numeq
                        pos[t][k+1][i] = @variable(model, lower_bound=0)
                    elseif eq_constraint_type[k-(m-numeq)] != -1
                        pos[t][k+1][i] = @variable(model)
                    end
                    if k <= m-numeq || eq_constraint_type[k-(m-numeq)] != -1
                        for (s, item) in enumerate(gsupp[I[t][k]])
                            @inbounds bi = [basis[t][k+1][blocks[t][k+1][i][1]][end:-1:1]; item; basis[t][k+1][blocks[t][k+1][i][1]]]
                            bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint)
                            Locb = bfind(tsupp, ltsupp, bi)
                            if k <= m-numeq || eq_constraint_type[k-(m-numeq)] == 1
                                @inbounds add_to_expression!(cons[Locb], gcoe[I[t][k]][s], pos[t][k+1][i])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*gcoe[I[t][k]][s], pos[t][k+1][i])
                            end
                        end
                    end
                else
                    if k <= m-numeq || eq_constraint_type[k-(m-numeq)] == 1
                        if k <= m-numeq
                            pos[t][k+1][i] = @variable(model, [1:bs, 1:bs], PSD)
                        else
                            pos[t][k+1][i] = @variable(model, [1:bs, 1:bs], Symmetric)
                        end
                        for j = 1:bs, r = j:bs, (s, item) in enumerate(gsupp[I[t][k]])
                            @inbounds bi = [basis[t][k+1][blocks[t][k+1][i][j]][end:-1:1]; item; basis[t][k+1][blocks[t][k+1][i][r]]]
                            bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint)
                            Locb = bfind(tsupp, ltsupp, bi)
                            if j == r
                                @inbounds add_to_expression!(cons[Locb], gcoe[I[t][k]][s], pos[t][k+1][i][j,r])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*gcoe[I[t][k]][s], pos[t][k+1][i][j,r])
                            end
                        end
                    elseif eq_constraint_type[k-(m-numeq)] == -1
                        pos[t][k+1][i] = @variable(model, [1:Int(bs*(bs-1)/2)])
                        for j = 1:bs-1, r = j+1:bs, (s, item) in enumerate(gsupp[I[t][k]])
                            @inbounds bi = [basis[t][k+1][blocks[t][k+1][i][j]][end:-1:1]; item; basis[t][k+1][blocks[t][k+1][i][r]]]
                            bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint)
                            Locb = bfind(tsupp, ltsupp, bi)
                            @inbounds add_to_expression!(cons[Locb], 2*gcoe[I[t][k]][s], pos[t][k+1][i][Int((2bs-j)*(j-1)/2)+r-j])
                        end
                    else
                        pos[t][k+1][i] = @variable(model, [1:bs, 1:bs])
                        for j = 1:bs, r = 1:bs, (s, item) in enumerate(gsupp[I[t][k]])
                            @inbounds bi = [basis[t][k+1][blocks[t][k+1][i][j]][end:-1:1]; item; basis[t][k+1][blocks[t][k+1][i][r]]]
                            bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint)
                            Locb = bfind(tsupp, ltsupp, bi)
                            @inbounds add_to_expression!(cons[Locb], 2*gcoe[I[t][k]][s], pos[t][k+1][i][j,r])
                        end
                    end
                end
            end
        end
    end
    bc = [AffExpr(0) for i=1:ltsupp]
    for (i, item) in enumerate(fsupp)
        Locb = bfind(tsupp, ltsupp, item)
        if Locb === nothing
            @error "The monomial basis is not enough!"
            return model,info
        else
            bc[Locb] = fcoe[i]
        end
    end
    if constrs !== nothing
        @constraint(model, cons==bc, base_name=constrs)
    else
        @constraint(model, cons==bc)
    end
    info = struct_data(cql,cliquesize,cliques,basis,cl,blocksize,blocks,tsupp,I,pos,constrs)
    return info
end

function get_blocks(I, m, fsupp::Vector{Vector{UInt16}}, gsupp::Vector{Vector{Vector{UInt16}}}, basis, cliques, cql; tsupp=[], TS="block", SO=1, QUIET=false, obj="eigen", partition=0, constraint=nothing)
    blocks = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    status = ones(Int, cql)
    if isempty(tsupp)
        tsupp = copy(fsupp)
        for i = 1:m
            append!(tsupp, gsupp[i])
        end
        tsupp = reduce!.(tsupp, obj=obj, partition=partition, constraint=constraint)
        sort!(tsupp)
        unique!(tsupp)
    end
    for i = 1:cql
        ind = [issubset(unique(item), cliques[i]) for item in tsupp]
        supp = copy(tsupp[ind])
        if obj == "trace"
            append!(supp, [_cyclic_canon([item[end:-1:1]; item]) for item in basis[i][1]])
        else
            append!(supp, [[item[end:-1:1]; item] for item in basis[i][1]])
        end
        supp = reduce!.(supp, obj=obj, partition=partition, constraint=constraint)
        sort!(supp)
        unique!(supp)
        blocks[i] = Vector{Vector{Vector{UInt16}}}(undef, length(I[i])+1)
        cl[i] = Vector{Int}(undef, length(I[i])+1)
        blocksize[i] = Vector{Vector{Int}}(undef, length(I[i])+1)
        blocks[i],cl[i],blocksize[i],status[i] = get_blocks(length(I[i]), supp, gsupp[I[i]], basis[i], TS=TS, SO=SO, obj=obj, partition=partition, constraint=constraint)
    end
    if minimum(status) == 1
        println("No higher TS step of the CS-TSSOS hierarchy!")
    end
    return blocks,cl,blocksize
end

function get_blocks(m::Int, tsupp, supp::Vector{Vector{Vector{UInt16}}}, basis::Vector{Vector{Vector{UInt16}}}; TS="block", SO=1, merge=false, md=3, obj="eigen", partition=0, constraint=nothing)
    blocks = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    blocksize = Vector{Vector{Int}}(undef, m+1)
    cl = Vector{Int}(undef, m+1)
    status = 0
    if TS == false
        for k = 1:m+1
            blocks[k],blocksize[k],cl[k] = [[i for i=1:length(basis[k])]],[length(basis[k])],1       
        end
    else       
        for i = 1:SO
            if i > 1
                oblocksize = deepcopy(blocksize)
            end
            for k = 1:m+1
                if k == 1
                    G = get_graph(tsupp, basis[1], obj=obj, partition=partition, constraint=constraint)
                else
                    G = get_graph(tsupp, supp[k-1], basis[k], obj=obj, partition=partition, constraint=constraint)
                end
                if TS == "block"
                    blocks[k] = connected_components(G)
                    blocksize[k] = length.(blocks[k])
                    cl[k] = length(blocksize[k])
                else
                    blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS)
                    if merge == true
                        blocks[k],cl[k],blocksize[k] = clique_merge!(blocks[k], d=md, QUIET=true)
                    end
                end
            end
            if i > 1 && blocksize == oblocksize
                status = 1
                break
            end
            if i < SO
                tsupp = Vector{UInt16}[]
                for t = 1:length(blocks[1]), j = 1:blocksize[1][t], r = j:blocksize[1][t]
                    ind1 = blocks[1][t][j]
                    ind2 = blocks[1][t][r]
                    @inbounds bi = [basis[1][ind1][end:-1:1]; basis[1][ind2]]
                    push!(tsupp, bi)
                end
                tsupp = reduce!.(tsupp, obj=obj, partition=partition, constraint=constraint)
                sort!(tsupp)
                unique!(tsupp)
            end
        end
    end
    return blocks,cl,blocksize,status
end

function assign_constraint(m, gsupp::Vector{Vector{Vector{UInt16}}}, cliques, cql)
    I = [UInt16[] for i=1:cql]
    for i = 1:m
        temp = UInt16[]
        for item in gsupp[i]
            append!(temp, item)
        end
        unique!(temp)
        ind = findall(k->issubset(temp, cliques[k]), 1:cql)
        push!.(I[ind], i)
    end
    return I
end

function clique_decomp(n, m, fsupp::Vector{Vector{UInt16}}, gsupp::Vector{Vector{Vector{UInt16}}}; alg="MF", QUIET=false)
    G = SimpleGraph(n)
    for item in fsupp
        add_clique!(G, unique(item))
    end
    for i = 1:m
        add_clique!(G, unique(Base.reduce(vcat, gsupp[i])))
    end
    if alg == "NC"
        cliques,cql,cliquesize = max_cliques(G)
    else
        cliques,cql,cliquesize = chordal_cliques!(G, method=alg, minimize=true)
    end
    uc = unique(cliquesize)
    sizes = [sum(cliquesize.== i) for i in uc]
    if QUIET == false
        println("-----------------------------------------------------------------------------")
        println("The clique sizes of varibles:\n$uc\n$sizes")
        println("-----------------------------------------------------------------------------")
    end
    return cliques,cql,cliquesize
end

function get_moment_matrix(moment, info; obj="eigen", partition=0, constraint=nothing)
    MomMat = Vector{Union{Float64, Symmetric{Float64}, Array{Float64,2}}}(undef, info.cql)
    ltsupp = length(info.tsupp)
    for i = 1:info.cql
        lb = length(info.basis[i][1])
        MomMat[i] = zeros(Float64, lb, lb)
        for j = 1:lb, k = j:lb
            bi = [info.basis[i][1][j][end:-1:1]; info.basis[i][1][k]]
            bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint)
            Locb = bfind(info.tsupp, ltsupp, bi)
            if Locb !== nothing
                MomMat[i][j,k] = moment[Locb]
            end
        end
        MomMat[i] = Symmetric(MomMat[i], :U)
    end
    return MomMat
end
