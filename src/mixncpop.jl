mutable struct ncmpop_data
    n::Int # number of all variables
    m::Int # number of all constraints
    numeq::Int # number of equality constraints
    supp # support data
    coe # coefficient data
    partition # the first 'partition' variables commute with the remaining variables
    comm_var # the first 'comm_var' variables commute each other
    constraint # nothing or "projection" or "unipotent"
    obj # "eigen" or "trace"
    ksupp # extended support at the k-th step
    basis # monomial bses
    cql # number of cliques
    cliques # cliques of variables
    cliquesize # numbers of cliques
    I # constraints associated to each clique
    blocks # block structure
    cl # numbers of blocks
    blocksize # sizes of blocks
    moment # moment matrix
    GramMat # Gram matrix
end

"""
    opt,data = cs_nctssos_first(pop, x, d; numeq=0, CS="MF", TS="block", solver="Mosek", writetofile=false, QUIET=false, 
    obj="eigen", partition=0, comm_var=0, constraint=nothing, solve=true, Gram=false)

Compute the first step of the CS-NCTSSOS hierarchy for constrained noncommutative polynomial
optimization with relaxation order `d`. Return the optimum and other auxiliary data.

# Arguments
- `pop`: the vector of the objective function, inequality constraints, and equality constraints
- `x`: the set of noncommuting variables
- `d`: the relaxation order of the moment-SOHS hierarchy
- `partition`: the first 'partition' variables commutes with the remaining variables
- `comm_var`: the first 'comm_var' variables commutes each other
- `constraint`: nothing or "projection" or "unipotent"
- `numeq`: the number of equality constraints
"""
function cs_nctssos_first(pop, x, d; numeq=0, CS="MF", minimize=false, TS="block", QUIET=false, obj="eigen", 
    solve=true, Gram=false, partition=0, comm_var=0, constraint=nothing, solver="Mosek", writetofile=false, cosmo_setting=cosmo_para())
    n,supp,coe = polys_info(pop, x)
    opt,data = cs_nctssos_first(supp, coe, n, d, numeq=numeq, CS=CS, minimize=minimize, TS=TS, QUIET=QUIET, obj=obj,
    solve=solve, solver=solver, writetofile=writetofile, Gram=Gram, partition=partition, comm_var=comm_var, constraint=constraint, cosmo_setting=cosmo_setting)
    return opt,data
end

"""
    opt,data = cs_nctssos_first(supp::Vector{Vector{Vector{UInt16}}}, coe, n::Int, d::Int; numeq=0,
    CS="MF", TS="block", QUIET=false, obj="eigen", partition=0, comm_var=0, constraint=nothing, solve=true, Gram=false)

Compute the first step of the CS-NCTSSOS hierarchy for constrained noncommutative polynomial optimization
with relaxation order `d`. Here the polynomial optimization problem is defined by `supp` and `coe`,
corresponding to the supports and coeffients of `pop` respectively. Return the optimum and other auxiliary data.

# Arguments
- `supp`: the supports of the polynomial optimization problem
- `coe`: the coeffients of the polynomial optimization problem
- `d`: the relaxation order of the moment-SOHS hierarchy
- `partition`: the first 'partition' variables commutes with the remaining variables
- `comm_var`: the first 'comm_var' variables commutes each other
- `constraint`: nothing or "projection" or "unipotent"
- `numeq`: the number of equality constraints
"""
function cs_nctssos_first(supp::Vector{Vector{Vector{UInt16}}}, coe, n::Int, d::Int; numeq=0, CS="MF", solver="Mosek", writetofile=false,
    minimize=false, TS="block", QUIET=false, obj="eigen", solve=true, Gram=false, partition=0, comm_var=0, constraint=nothing, cosmo_setting=cosmo_para())
    println("********************************** NCTSSOS **********************************")
    println("NCTSSOS is launching...")
    m = length(supp)-1
    dg = [maximum(length.(supp[i])) for i=2:m+1]
    if obj == "trace"
        supp[1],coe[1] = cyclic_canon(supp[1], coe[1])
    else
        supp[1],coe[1] = sym_canon(supp[1], coe[1])
    end
    time = @elapsed begin
    cliques,cql,cliquesize = clique_decomp(n, m, supp, alg=CS, minimize=minimize)
    end
    if CS != false && QUIET == false
        mc = maximum(cliquesize)
        println("Obtained the variable cliques in $time seconds.\nThe maximal size of cliques is $mc.")
    end
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    I = assign_constraint(m, supp, cliques, cql)
    basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
    for t = 1:cql
        basis[t] = Vector{Vector{Vector{UInt16}}}(undef, length(I[t])+1)
        basis[t][1] = get_ncbasis(cliquesize[t], d, ind=cliques[t], binary=constraint!==nothing)
        for s = 1:length(I[t])
            basis[t][s+1] = get_ncbasis(cliquesize[t], d-ceil(Int, dg[I[t][s]]/2), ind=cliques[t], binary=constraint!==nothing)
        end
        for ba in basis[t]
            ind = [_comm(item, partition, comm_var) == item for item in ba]
            ba = ba[ind]
        end
    end
    tsupp = nothing
    if TS != false
        tsupp = copy(supp[1])
        for i = 1:m
            append!(tsupp, supp[i+1])
        end
        tsupp = reduce!.(tsupp, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
        sort!(tsupp)
        unique!(tsupp)
    end
    time = @elapsed begin
    blocks,cl,blocksize = get_blocks(I, supp, tsupp, basis, cliques, cql, TS=TS, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,GramMat = solvesdp(m, supp, coe, basis, cql, I, blocks, cl, blocksize, numeq=numeq, QUIET=QUIET, obj=obj, solve=solve, 
    solver=solver, writetofile=writetofile, Gram=Gram, partition=partition, comm_var=comm_var, constraint=constraint, cosmo_setting=cosmo_setting)
    data = ncmpop_data(n, m, numeq, supp, coe, partition, comm_var, constraint, obj, ksupp, basis, cql, cliques, cliquesize, I, blocks, cl, blocksize, 
    moment,GramMat)
    return opt,data
end

"""
    opt,data = cs_nctssos_higher!(data; TS="block", QUIET=false, solve=true, Gram=false)

Compute higher steps of the CS-NCTSSOS hierarchy.
Return the optimum and other auxiliary data.
"""
function cs_nctssos_higher!(data::ncmpop_data; TS="block", QUIET=false, solve=true, Gram=false, solver="Mosek", writetofile=false, cosmo_setting=cosmo_para())
    n = data.n
    m = data.m
    numeq = data.numeq
    supp = data.supp
    coe = data.coe
    partition = data.partition
    comm_var = data.comm_var
    constraint = data.constraint
    obj = data.obj
    ksupp = data.ksupp
    basis = data.basis
    cql = data.cql
    cliques = data.cliques
    I = data.I
    blocks = data.blocks
    cl = data.cl
    blocksize = data.blocksize
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    oblocksize = deepcopy(data.blocksize)
    time = @elapsed begin
    blocks,cl,blocksize = get_blocks(I, supp, ksupp, basis, cliques, cql, blocks=blocks, cl=cl, 
    blocksize=blocksize, TS=TS, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    end
    if blocksize == oblocksize
        println("No higher TS step of the CS-NCTSSOS hierarchy!")
        opt = nothing
    else
        if QUIET == false
            mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,GramMat = solvesdp(m, supp, coe, basis, cql, I, blocks, cl, blocksize, numeq=numeq, QUIET=QUIET, obj=obj, solve=solve, 
        Gram=Gram, partition=partition, comm_var=comm_var, constraint=constraint, solver=solver, writetofile=writetofile, cosmo_setting=cosmo_setting)
        data.moment = moment
        data.GramMat = GramMat
        data.ksupp = ksupp
        data.blocks = blocks
        data.cl = cl
        data.blocksize = blocksize
    end
    return opt,data
end

function solvesdp(m::Int, supp::Vector{Vector{Vector{UInt16}}}, coe, basis, cql, I, blocks, cl, blocksize; numeq=0, QUIET=false, obj="eigen", solve=true, Gram=false, 
    partition=0, comm_var=0, constraint=nothing, solver="Mosek", writetofile=false, cosmo_setting=cosmo_para())
    ksupp = Vector{UInt16}[]
    for i = 1:cql
        for j = 1:cl[i][1], k = 1:blocksize[i][1][j], r = k:blocksize[i][1][j]
            @inbounds bi = [basis[i][1][blocks[i][1][j][k]][end:-1:1]; basis[i][1][blocks[i][1][j][r]]]
            push!(ksupp, bi)
        end
        for (j, w) in enumerate(I[i])
            for l = 1:cl[i][j+1], t = 1:blocksize[i][j+1][l], r = t:blocksize[i][j+1][l], s = 1:length(supp[w+1])
                ind1 = blocks[i][j+1][l][t]
                ind2 = blocks[i][j+1][l][r]
                @inbounds bi = [basis[i][j+1][ind1][end:-1:1]; supp[w+1][s]; basis[i][j+1][ind2]]
                push!(ksupp, bi)
            end
        end
    end
    ksupp = reduce!.(ksupp, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    sort!(ksupp)
    unique!(ksupp)
    lksupp = length(ksupp)
    if QUIET == false
        println("There are $lksupp affine constraints.")
    end
    objv = moment = GramMat = nothing
    if solve == true
        if QUIET == false
            println("Assembling the SDP...")
        end
        if solver == "Mosek"
            model = Model(optimizer_with_attributes(Mosek.Optimizer))
        elseif solver == "COSMO"
            model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter))
        else
            @error "The solver is currently not supported!"
            return nothing,nothing,nothing,nothing
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = [AffExpr(0) for i=1:lksupp]
        for i = 1:cql, l = 1:cl[i][1]
            if blocksize[i][1][l] == 1
               @inbounds pos = @variable(model, lower_bound=0)
               @inbounds bi = [basis[i][1][blocks[i][1][l][1]][end:-1:1]; basis[i][1][blocks[i][1][l][1]]]
               bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
               Locb = bfind(ksupp, lksupp,bi)
               @inbounds add_to_expression!(cons[Locb], pos)
            else
               @inbounds bs = blocksize[i][1][l]
               @inbounds pos = @variable(model, [1:bs, 1:bs], PSD)
               for t = 1:bs, r = t:bs
                   @inbounds ind1 = blocks[i][1][l][t]
                   @inbounds ind2 = blocks[i][1][l][r]
                   @inbounds bi = [basis[i][1][ind1][end:-1:1]; basis[i][1][ind2]]
                   bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                   Locb = bfind(ksupp, lksupp, bi)
                   if t == r
                      @inbounds add_to_expression!(cons[Locb], pos[t,r])
                   else
                      @inbounds add_to_expression!(cons[Locb], 2, pos[t,r])
                   end
               end
            end
        end
        for i = 1:cql, (j, w) in enumerate(I[i])
            for l = 1:cl[i][j+1]
                bs = blocksize[i][j+1][l]
                if bs == 1
                    if j <= m-numeq
                        pos = @variable(model, lower_bound=0)
                    else
                        pos = @variable(model)
                    end
                    ind1 = blocks[i][j+1][l][1]
                    for s = 1:length(supp[w+1])
                        @inbounds bi = [basis[i][j+1][ind1][end:-1:1]; supp[w+1][s]; basis[i][j+1][ind1]]
                        bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                        Locb = bfind(ksupp, lksupp, bi)
                        @inbounds add_to_expression!(cons[Locb], coe[w+1][s], pos)
                    end
                else
                    if j <= m-numeq
                        pos = @variable(model, [1:bs, 1:bs], PSD)
                    else
                        pos = @variable(model, [1:bs, 1:bs], Symmetric)
                    end
                    for t = 1:bs, r = t:bs
                        ind1 = blocks[i][j+1][l][t]
                        ind2 = blocks[i][j+1][l][r]
                        for s = 1:length(supp[w+1])
                            @inbounds bi = [basis[i][j+1][ind1][end:-1:1]; supp[w+1][s]; basis[i][j+1][ind2]]
                            bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                            Locb = bfind(ksupp, lksupp, bi)
                            if t == r
                                @inbounds add_to_expression!(cons[Locb], coe[w+1][s], pos[t,r])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*coe[w+1][s], pos[t,r])
                            end
                        end
                    end
                end
            end
        end
        bc = zeros(lksupp)
        for i = 1:length(supp[1])
            Locb = bfind(ksupp, lksupp, supp[1][i])
            if Locb === nothing
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing,nothing
            else
               bc[Locb] = coe[1][i]
            end
        end
        @variable(model, lower)
        cons[1] += lower
        @constraint(model, con, cons.==bc)
        @objective(model, Max, lower)
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
            println("Solving the SDP...")
        end
        time=@elapsed begin
        optimize!(model)
        end
        if QUIET == false
            println("SDP solving time: $time seconds.")
        end
        if writetofile != false
            write_to_file(dualize(model), writetofile)
        end
        status = termination_status(model)
        objv = objective_value(model)
        if status != MOI.OPTIMAL
           println("termination status: $status")
           status = primal_status(model)
           println("solution status: $status")
        end
        println("optimum = $objv")
        dual_var = -dual.(con)
        moment = Vector{Vector{Matrix{Float64}}}(undef, cql)
        for i = 1:cql
            moment[i] = Vector{Matrix{Float64}}(undef, cl[i][1])
            for k = 1:cl[i][1]
                moment[i][k] = zeros(blocksize[i][1][k],blocksize[i][1][k])
                for j = 1:blocksize[i][1][k], r = j:blocksize[i][1][k]
                    @inbounds bi = [basis[i][1][blocks[i][1][k][j]][end:-1:1]; basis[i][1][blocks[i][1][k][r]]]
                    bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                    Locb = bfind(ksupp, lksupp, bi)
                    moment[i][k][j,r] = dual_var[Locb]
                end
                moment[i][k] = Symmetric(moment[i][k],:U)
            end
        end
    end
    return objv,ksupp,moment,GramMat
end

function get_blocks(I, supp::Vector{Vector{Vector{UInt16}}}, tsupp, basis, cliques, cql; TS="block", blocks=[], cl=[], blocksize=[], obj="eigen", partition=0, comm_var=0, constraint=nothing)
    if isempty(blocks)
        blocks = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        cl = Vector{Vector{Int}}(undef, cql)
        blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
        ksupp = nothing
        for i = 1:cql
            if TS != false
                ind = [issubset(unique(item), cliques[i]) for item in tsupp]
                ksupp = copy(tsupp[ind])
                if obj == "trace"
                    append!(ksupp, [_cyclic_canon([item[end:-1:1]; item]) for item in basis[i][1]])
                else
                    append!(ksupp, [[item[end:-1:1]; item] for item in basis[i][1]])
                end
                ksupp = reduce!.(ksupp, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                sort!(ksupp)
                unique!(ksupp)
            end
            blocks[i] = Vector{Vector{Vector{UInt16}}}(undef, length(I[i])+1)
            cl[i] = Vector{Int}(undef, length(I[i])+1)
            blocksize[i] = Vector{Vector{Int}}(undef, length(I[i])+1)
            blocks[i],cl[i],blocksize[i] = get_cblocks(length(I[i]), ksupp, supp[I[i].+1], basis[i], QUIET=true, TS=TS, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
        end
    else
        for i = 1:cql
            ind = [issubset(unique(item), cliques[i]) for item in tsupp]
            blocks[i],cl[i],blocksize[i] = get_cblocks(length(I[i]), tsupp[ind], supp[I[i].+1], basis[i], QUIET=true, TS=TS, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
        end
    end
    return blocks,cl,blocksize
end

function assign_constraint(m, supp, cliques, cql)
    I = [UInt16[] for i=1:cql]
    for i = 1:m
        ind = findall(k->issubset(unique(Base.reduce(vcat, supp[i+1])), cliques[k]), 1:cql)
        push!.(I[ind], i)
    end
    return I
end

function clique_decomp(n::Int, m::Int, supp::Vector{Vector{Vector{UInt16}}}; alg="MF", minimize=false)
    if alg == false
        cliques = [UInt16[i for i=1:n]]
        cql = 1
        cliquesize=[n]
    else
        G = SimpleGraph(n)
        for item in supp[1]
            add_clique!(G, unique(item))
        end
        for i = 2:m+1
            add_clique!(G, unique(Base.reduce(vcat, supp[i])))
        end
        if alg == "NC"
            cliques,cql,cliquesize = max_cliques(G)
        else
            cliques,cql,cliquesize = chordal_cliques!(G, method=alg, minimize=minimize)
        end
    end
    uc = unique(cliquesize)
    sizes = [sum(cliquesize.== i) for i in uc]
    println("-----------------------------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("-----------------------------------------------------------------------------")
    return cliques,cql,cliquesize
end
