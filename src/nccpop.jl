mutable struct nccpop_data
    n # number of variables
    m # number of constraints
    order # relaxation order
    numeq # number of equality constraints
    supp # support data
    coe # coefficient data
    partition # the first 'partition' variables commutes with the remaining variables
    constraint # nothing or "projection" or "unipotent"
    obj # "eigen" or "trace"
    basis # monomial bases
    ksupp # extended support at the k-th step
    blocks # block structure
    cl # numbers of blocks
    blocksize # sizes of blocks
    moment # moment matrix
    GramMat # Gram matrix
end

"""
    opt,data = nctssos_first(pop::Vector{Polynomial{false, T}} where T<:Number, x::Vector{PolyVar{false}},
        order::Int; numeq=0, reducebasis=false, TS="block", normality=false, obj="eigen", merge=false, md=3, 
        solve=true, Gram=false, QUIET=false)

Compute the first step of the NCTSSOS hierarchy for constrained noncommutative polynomial optimization with
relaxation order `order`.
Return the optimum and other auxiliary data.

# Arguments
- `pop`: the vector of the objective function, inequality constraints, and equality constraints.
- `x`: the set of noncommuting variables.
- `order`: the relaxation order of the moment-SOHS hierarchy.
- `numeq`: the number of equality constraints.
"""

function nctssos_first(
    pop::Vector{Polynomial{false,T}} where {T<:Number},
    x::Vector{PolyVar{false}},
    order::Int;
    numeq=0,
    reducebasis=false,
    TS="block",
    soc=false,
    obj="eigen",
    merge=false,
    md=3,
    solve=true,
    Gram=false,
    QUIET=false,
    solver="Mosek",
    partition=0,
    constraint=nothing,
    normality=false,
    cosmo_setting=cosmo_para(),
)
    n, supp, coe = polys_info(pop, x)
    opt, data = nctssos_first(
        supp,
        coe,
        n,
        order;
        numeq=numeq,
        reducebasis=reducebasis,
        TS=TS,
        obj=obj,
        merge=merge,
        soc=soc,
        md=md,
        QUIET=QUIET,
        solve=solve,
        solver=solver,
        Gram=Gram,
        partition=partition,
        constraint=constraint,
        normality=normality,
        cosmo_setting=cosmo_setting,
    )
    return opt, data
end

function nctssos_first(
    supp::Vector{Vector{Vector{UInt16}}},
    coe,
    n::Int64,
    order::Int64;
    numeq=0,
    reducebasis=false,
    TS="block",
    soc=false,
    obj="eigen",
    merge=false,
    md=3,
    solve=true,
    solver="Mosek",
    Gram=false,
    QUIET=false,
    partition=0,
    constraint=nothing,
    normality=false,
    cosmo_setting=cosmo_para(),
)
    println("********************************** NCTSSOS **********************************")
    println("NCTSSOS is launching...")
    m = length(supp) - 1
    dg = [maximum(length.(supp[i])) for i in 2:(m + 1)]
    if obj == "trace"
        supp[1], coe[1] = cyclic_canon(supp[1], coe[1])
    else
        supp[1], coe[1] = sym_canon(supp[1], coe[1])
    end
    basis = Vector{Vector{Vector{UInt16}}}(undef, m + 1)

    basis[1] = get_ncbasis(n, order; binary=constraint !== nothing)
    # TODO: I did not consider some variables are commutative
    if partition > 0
        ind = [_comm(basis[1][i], partition) == basis[1][i] for i in 1:length(basis[1])]
        basis[1] = basis[1][ind]
    end

    # FIXME: I am here
    ksupp = copy(supp[1])
    for i in 1:m
        basis[i + 1] = get_ncbasis(
            n, order - Int(ceil(dg[i] / 2)); binary=constraint !== nothing
        )
        if partition > 0
            ind = [_comm(item, partition) == item for item in basis[i + 1]]
            basis[i + 1] = basis[i + 1][ind]
        end
        if obj == "trace"
            append!(
                ksupp,
                [
                    min(_cyclic_canon(word), _cyclic_canon(reverse(word))) for
                    word in supp[i + 1]
                ],
            )
        else
            append!(ksupp, _sym_canon.(supp[i + 1]))
        end
    end
    if obj == "trace"
        append!(ksupp, [_cyclic_canon([item[end:-1:1]; item]) for item in basis[1]])
    else
        append!(ksupp, [[item[end:-1:1]; item] for item in basis[1]])
    end
    ksupp = reduce!.(ksupp, obj=obj, partition=partition, constraint=constraint)
    sort!(ksupp)
    unique!(ksupp)
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
        blocks, cl, blocksize = get_blocks(
            m,
            ksupp,
            supp[2:end],
            basis;
            TS=TS,
            obj=obj,
            QUIET=QUIET,
            merge=merge,
            md=md,
            partition=partition,
            constraint=constraint,
        )
        if reducebasis == true && obj == "eigen" && constraint === nothing
            gsupp = get_gsupp(
                m, supp, basis[2:end], blocks[2:end], cl[2:end], blocksize[2:end]
            )
            psupp = copy(supp[1])
            push!(psupp, UInt16[])
            append!(psupp, gsupp)
            psupp = psupp[is_sym.(psupp)]
            basis[1], flag = reducebasis!(psupp, basis[1], blocks[1], cl[1], blocksize[1])
            if flag == 1
                ksupp = copy(supp[1])
                for i in 1:m
                    append!(ksupp, _sym_canon.(supp[i + 1]))
                end
                append!(ksupp, [[item[end:-1:1]; item] for item in basis[1]])
                if partition > 0
                    ksupp = _comm.(ksupp, partition)
                end
                sort!(ksupp)
                unique!(ksupp)
                blocks, cl, blocksize = get_blocks(
                    m,
                    ksupp,
                    supp[2:end],
                    basis;
                    TS=TS,
                    obj=obj,
                    QUIET=QUIET,
                    merge=merge,
                    md=md,
                )
            end
        end
    end
    if QUIET == false
        mb = maximum(maximum.(blocksize))
        println(
            "Obtained the block structure in $time seconds. The maximal size of blocks is $mb.",
        )
    end
    opt, ksupp, moment, GramMat = solvesdp(
        order,
        n,
        m,
        supp,
        coe,
        basis,
        blocks,
        cl,
        blocksize;
        numeq=numeq,
        QUIET=QUIET,
        obj=obj,
        soc=soc,
        TS=TS,
        solve=solve,
        solver=solver,
        Gram=Gram,
        partition=partition,
        constraint=constraint,
        normality=normality,
        cosmo_setting=cosmo_setting,
    )
    data = nccpop_data(
        n,
        m,
        order,
        numeq,
        supp,
        coe,
        partition,
        constraint,
        obj,
        basis,
        ksupp,
        blocks,
        cl,
        blocksize,
        moment,
        GramMat,
    )
    return opt, data
end

function nctssos_higher!(
    data::nccpop_data;
    TS="block",
    soc=false,
    merge=false,
    md=3,
    solve=true,
    solver="Mosek",
    Gram=false,
    QUIET=false,
    normality=false,
    cosmo_setting=cosmo_para(),
)
    n = data.n
    m = data.m
    numeq = data.numeq
    supp = data.supp
    coe = data.coe
    partition = data.partition
    constraint = data.constraint
    obj = data.obj
    basis = data.basis
    ksupp = data.ksupp
    blocks = data.blocks
    cl = data.cl
    blocksize = data.blocksize
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    oblocksize = deepcopy(data.blocksize)
    time = @elapsed begin
        blocks, cl, blocksize = get_blocks(
            m,
            ksupp,
            supp[2:end],
            basis;
            blocks=blocks,
            cl=cl,
            blocksize=blocksize,
            TS=TS,
            obj=obj,
            QUIET=QUIET,
            merge=merge,
            md=md,
            partition=partition,
            constraint=constraint,
        )
    end
    if blocksize == oblocksize
        println("No higher TS step of the NCTSSOS hierarchy!")
        opt = nothing
    else
        if QUIET == false
            mb = maximum(maximum.(blocksize))
            println(
                "Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.",
            )
        end
        opt, ksupp, moment, GramMat = solvesdp(
            data.order,
            n,
            m,
            supp,
            coe,
            basis,
            blocks,
            cl,
            blocksize;
            numeq=numeq,
            QUIET=QUIET,
            obj=obj,
            soc=soc,
            TS=TS,
            solve=solve,
            solver=solver,
            Gram=Gram,
            partition=partition,
            constraint=constraint,
            normality=normality,
            cosmo_setting=cosmo_setting,
        )
        data.moment = moment
        data.GramMat = GramMat
        data.ksupp = ksupp
        data.blocks = blocks
        data.cl = cl
        data.blocksize = blocksize
    end
    return opt, data
end

function get_gsupp(m, supp, gbasis, gblocks, gcl, gblocksize)
    gsupp = Vector{UInt16}[]
    for k in 1:m,
        i in 1:gcl[k],
        j in 1:gblocksize[k][i],
        r in j:gblocksize[k][i],
        s in 1:length(supp[k + 1])

        @inbounds bi = [
            gbasis[k][gblocks[k][i][j]][end:-1:1]
            supp[k + 1][s]
            gbasis[k][gblocks[k][i][r]]
        ]
        push!(gsupp, bi)
    end
    return gsupp
end

function reducebasis!(psupp, basis, blocks, cl, blocksize)
    init = 0
    flag = 0
    check = 0
    while init == 0 || check > 0
        init = 1
        check = 0
        for i in 1:cl
            if blocksize[i] > 1
                for j in 1:blocksize[i], r in 1:blocksize[i]
                    if j != r
                        @inbounds bi = [basis[blocks[i][j]][end:-1:1]; basis[blocks[i][r]]]
                        if is_sym(bi)
                            push!(psupp, bi)
                        end
                    end
                end
            end
        end
        sort!(psupp)
        unique!(psupp)
        lpsupp = length(psupp)
        for i in 1:cl
            lo = blocksize[i]
            indexb = [k for k in 1:lo]
            j = 1
            while lo >= j
                bi = [basis[blocks[i][indexb[j]]][end:-1:1]; basis[blocks[i][indexb[j]]]]
                Locb = bfind(psupp, lpsupp, bi)
                if Locb === nothing
                    check = 1
                    flag = 1
                    deleteat!(indexb, j)
                    lo = lo - 1
                else
                    j += 1
                end
            end
            blocks[i] = blocks[i][indexb]
            blocksize[i] = lo
        end
    end
    if flag == 1
        indexb = blocks[1]
        for i in 2:cl
            indexb = append!(indexb, blocks[i])
        end
        sort!(indexb)
        unique!(indexb)
        return basis[indexb], flag
    else
        return basis, flag
    end
end

function get_graph(ksupp, supp, basis; obj="eigen", partition=0, constraint=nothing)
    lb = length(basis)
    lksupp = length(ksupp)
    G = SimpleGraph(lb)
    for i in 1:lb, j in (i + 1):lb
        r = 1
        while r <= length(supp)
            bi = [basis[i][end:-1:1]; supp[r]; basis[j]]
            bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
            if bfind(ksupp, lksupp, bi) !== nothing
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

function get_blocks(
    m,
    ksupp,
    gsupp,
    basis;
    blocks=[],
    cl=[],
    blocksize=[],
    TS="block",
    obj="eigen",
    QUIET=true,
    merge=false,
    md=3,
    partition=0,
    constraint=nothing,
)
    if isempty(blocks)
        blocks = Vector{Vector{Vector{UInt16}}}(undef, m + 1)
        blocksize = Vector{Vector{UInt16}}(undef, m + 1)
        cl = Vector{UInt16}(undef, m + 1)
    end
    if TS == false
        for k in 1:(m + 1)
            blocks[k] = [[i for i in 1:length(basis[k])]]
            blocksize[k] = [length(basis[k])]
            cl[k] = 1
        end
    else
        for k in 1:(m + 1)
            if k == 1
                G = get_graph(
                    ksupp, basis[1]; obj=obj, partition=partition, constraint=constraint
                )
            else
                G = get_graph(
                    ksupp,
                    gsupp[k - 1],
                    basis[k];
                    obj=obj,
                    partition=partition,
                    constraint=constraint,
                )
            end
            if TS == "block"
                blocks[k] = connected_components(G)
                blocksize[k] = length.(blocks[k])
                cl[k] = length(blocksize[k])
            else
                blocks[k], cl[k], blocksize[k] = chordal_cliques!(
                    G; method=TS, minimize=false
                )
                if merge == true
                    blocks[k], cl[k], blocksize[k] = clique_merge!(
                        blocks[1]; d=md, QUIET=true
                    )
                end
            end
            if QUIET == false
                sb = sort(Int.(unique(blocksize[k])); rev=true)
                numb = [sum(blocksize[k] .== i) for i in sb]
                println(
                    "-----------------------------------------------------------------------------",
                )
                println(
                    "The sizes of PSD blocks for the $k-th SOHS multiplier:\n$sb\n$numb"
                )
                println(
                    "-----------------------------------------------------------------------------",
                )
            end
        end
    end
    return blocks, cl, blocksize
end

function solvesdp(
    order,
    n,
    m,
    supp,
    coe,
    basis,
    blocks,
    cl,
    blocksize;
    numeq=0,
    QUIET=true,
    obj="eigen",
    solve=true,
    solver="Mosek",
    TS="block",
    Gram=false,
    soc=false,
    partition=0,
    constraint=nothing,
    normality=false,
    cosmo_setting=cosmo_para(),
)
    ksupp = Vector{UInt16}[]
    for i in 1:cl[1], j in 1:blocksize[1][i], r in j:blocksize[1][i]
        @inbounds bi = [basis[1][blocks[1][i][j]][end:-1:1]; basis[1][blocks[1][i][r]]]
        push!(ksupp, bi)
    end
    if TS != false
        gsupp = get_gsupp(m, supp, basis[2:end], blocks[2:end], cl[2:end], blocksize[2:end])
        append!(ksupp, gsupp)
    end
    if soc == true
        df = maximum(length.(supp[1]))
        slb = basis[1][length.(basis[1]) .<= 2 * order - df]
        spb = basis[1][length.(basis[1]) .<= order - Int(ceil(df / 2))]
        lspb = length(spb)
        if TS != false
            for item in slb, s in 1:length(supp[1])
                @inbounds bi = [supp[1][s]; item]
                bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                push!(ksupp, bi)
                @inbounds bi = [item; supp[1][s]]
                bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                push!(ksupp, bi)
            end
            for i in 1:lspb, j in i:lspb, s in 1:length(supp[1])
                @inbounds bi = [spb[i][end:-1:1]; supp[1][s]; spb[j]]
                bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                push!(ksupp, bi)
                @inbounds bi = [supp[1][s]; spb[i][end:-1:1]; spb[j]]
                bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                push!(ksupp, bi)
                @inbounds bi = [spb[i][end:-1:1]; spb[j]; supp[1][s]]
                bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                push!(ksupp, bi)
            end
        end
    end
    if normality == true
        wbasis = basis[1]
        bs = length(wbasis)
        for i in 1:n, j in 1:bs, k in j:bs
            @inbounds bi = [wbasis[j][end:-1:1]; i; wbasis[k]]
            push!(ksupp, bi)
            @inbounds bi = [wbasis[j][end:-1:1]; i; i; wbasis[k]]
            push!(ksupp, bi)
        end
    end
    ksupp = reduce!.(ksupp, obj=obj, partition=partition, constraint=constraint)
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
            model = Model(
                optimizer_with_attributes(
                    COSMO.Optimizer,
                    "eps_abs" => cosmo_setting.eps_abs,
                    "eps_rel" => cosmo_setting.eps_rel,
                    "max_iter" => cosmo_setting.max_iter,
                ),
            )
        else
            @error "The solver is currently not supported!"
            return nothing, nothing, nothing, nothing
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
            cons = [AffExpr(0) for i in 1:lksupp]
            if normality == true
                for i in 1:n
                    hnom = @variable(model, [1:(2bs), 1:(2bs)], PSD)
                    for j in 1:bs, k in j:bs
                        @inbounds bi = [wbasis[j][end:-1:1]; wbasis[k]]
                        bi = reduce!(
                            bi; obj=obj, partition=partition, constraint=constraint
                        )
                        Locb = bfind(ksupp, lksupp, bi)
                        if j == k
                            @inbounds add_to_expression!(cons[Locb], hnom[j, k])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2, hnom[j, k])
                        end
                        @inbounds bi = UInt16[wbasis[j][end:-1:1]; i; i; wbasis[k]]
                        bi = reduce!(
                            bi; obj=obj, partition=partition, constraint=constraint
                        )
                        Locb = bfind(ksupp, lksupp, bi)
                        if j == k
                            @inbounds add_to_expression!(cons[Locb], hnom[j + bs, k + bs])
                        else
                            @inbounds add_to_expression!(
                                cons[Locb], 2, hnom[j + bs, k + bs]
                            )
                        end
                        @inbounds bi = UInt16[wbasis[j][end:-1:1]; i; wbasis[k]]
                        bi = reduce!(
                            bi; obj=obj, partition=partition, constraint=constraint
                        )
                        Locb = bfind(ksupp, lksupp, bi)
                        if j == k
                            @inbounds add_to_expression!(cons[Locb], 2, hnom[j, k + bs])
                        else
                            @inbounds add_to_expression!(
                                cons[Locb], 2, hnom[j, k + bs] + hnom[k, j + bs]
                            )
                        end
                    end
                end
            end
            pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1])
            for i in 1:cl[1]
                bs = blocksize[1][i]
                if bs == 1
                    @inbounds pos[i] = @variable(model, lower_bound = 0)
                    @inbounds bi = [
                        basis[1][blocks[1][i][1]][end:-1:1]
                        basis[1][blocks[1][i][1]]
                    ]
                    bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                    Locb = bfind(ksupp, lksupp, bi)
                    @inbounds add_to_expression!(cons[Locb], pos[i])
                else
                    @inbounds pos[i] = @variable(model, [1:bs, 1:bs], PSD)
                    for j in 1:bs, r in j:bs
                        @inbounds bi = [
                            basis[1][blocks[1][i][j]][end:-1:1]
                            basis[1][blocks[1][i][r]]
                        ]
                        bi = reduce!(
                            bi; obj=obj, partition=partition, constraint=constraint
                        )
                        Locb = bfind(ksupp, lksupp, bi)
                        if j == r
                            @inbounds add_to_expression!(cons[Locb], pos[i][j, r])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2, pos[i][j, r])
                        end
                    end
                end
            end
            gpos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m)
            for k in 1:m
                gpos[k] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(
                    undef, cl[k + 1]
                )
                for i in 1:cl[k + 1]
                    bs = blocksize[k + 1][i]
                    if bs == 1
                        if k <= m - numeq
                            gpos[k][i] = @variable(model, lower_bound = 0)
                        else
                            gpos[k][i] = @variable(model)
                        end
                        for s in 1:length(supp[k + 1])
                            @inbounds bi = [
                                basis[k + 1][blocks[k + 1][i][1]][end:-1:1]
                                supp[k + 1][s]
                                basis[k + 1][blocks[k + 1][i][1]]
                            ]
                            bi = reduce!(
                                bi; obj=obj, partition=partition, constraint=constraint
                            )
                            Locb = bfind(ksupp, lksupp, bi)
                            @inbounds add_to_expression!(
                                cons[Locb], coe[k + 1][s], gpos[k][i]
                            )
                        end
                    else
                        if k <= m - numeq
                            gpos[k][i] = @variable(model, [1:bs, 1:bs], PSD)
                        else
                            gpos[k][i] = @variable(model, [1:bs, 1:bs], Symmetric)
                        end
                        for j in 1:bs, r in j:bs, s in 1:length(supp[k + 1])
                            @inbounds bi = [
                                basis[k + 1][blocks[k + 1][i][j]][end:-1:1]
                                supp[k + 1][s]
                                basis[k + 1][blocks[k + 1][i][r]]
                            ]
                            bi = reduce!(
                                bi; obj=obj, partition=partition, constraint=constraint
                            )
                            Locb = bfind(ksupp, lksupp, bi)
                            if j == r
                                @inbounds add_to_expression!(
                                    cons[Locb], coe[k + 1][s], gpos[k][i][j, r]
                                )
                            else
                                @inbounds add_to_expression!(
                                    cons[Locb], 2 * coe[k + 1][s], gpos[k][i][j, r]
                                )
                            end
                        end
                    end
                end
            end
            if soc == true
                fr = @variable(model, [1:length(slb)])
                for (i, item) in enumerate(slb), s in 1:length(supp[1])
                    @inbounds bi = [supp[1][s]; item]
                    bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                    Locb = bfind(ksupp, lksupp, bi)
                    add_to_expression!(cons[Locb], coe[1][s], fr[i])
                    @inbounds bi = [item; supp[1][s]]
                    bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                    Locb = bfind(ksupp, lksupp, bi)
                    add_to_expression!(cons[Locb], -coe[1][s], fr[i])
                end
                pos = @variable(model, [1:lspb, 1:lspb], PSD)
                for i in 1:lspb, j in i:lspb, s in 1:length(supp[1])
                    @inbounds bi = [spb[i][end:-1:1]; supp[1][s]; spb[j]]
                    bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                    Locb = bfind(ksupp, lksupp, bi)
                    if i == j
                        @inbounds add_to_expression!(cons[Locb], coe[1][s], pos[i, j])
                    else
                        @inbounds add_to_expression!(cons[Locb], 2 * coe[1][s], pos[i, j])
                    end
                    @inbounds bi = [supp[1][s]; spb[i][end:-1:1]; spb[j]]
                    bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                    Locb = bfind(ksupp, lksupp, bi)
                    if i == j
                        @inbounds add_to_expression!(
                            cons[Locb], -0.5 * coe[1][s], pos[i, j]
                        )
                    else
                        @inbounds add_to_expression!(cons[Locb], -coe[1][s], pos[i, j])
                    end
                    @inbounds bi = [spb[i][end:-1:1]; spb[j]; supp[1][s]]
                    bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                    if i == j
                        @inbounds add_to_expression!(
                            cons[Locb], -0.5 * coe[1][s], pos[i, j]
                        )
                    else
                        @inbounds add_to_expression!(cons[Locb], -coe[1][s], pos[i, j])
                    end
                end
            end
            bc = zeros(lksupp)
            for i in 1:length(supp[1])
                Locb = bfind(ksupp, lksupp, supp[1][i])
                if Locb === nothing
                    @error "The monomial basis is not enough!"
                    return nothing, nothing, nothing, nothing
                else
                    bc[Locb] = coe[1][i]
                end
            end
            @variable(model, lower)
            cons[1] += lower
            @constraint(model, con, cons .== bc)
            @objective(model, Max, lower)
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
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
            GramMat = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, m + 1)
            GramMat[1] = [value.(pos[i]) for i in 1:cl[1]]
            for k in 1:m
                GramMat[k + 1] = [value.(gpos[k][i]) for i in 1:cl[k + 1]]
            end
        end
        dual_var = -dual.(con)
        moment = Vector{Matrix{Float64}}(undef, cl[1])
        for i in 1:cl[1]
            moment[i] = zeros(blocksize[1][i], blocksize[1][i])
            for j in 1:blocksize[1][i], k in j:blocksize[1][i]
                @inbounds bi = [
                    basis[1][blocks[1][i][j]][end:-1:1]
                    basis[1][blocks[1][i][k]]
                ]
                bi = reduce!(bi; obj=obj, partition=partition, constraint=constraint)
                Locb = bfind(ksupp, lksupp, bi)
                moment[i][j, k] = dual_var[Locb]
            end
            moment[i] = Symmetric(moment[i], :U)
        end
    end
    return objv, ksupp, moment, GramMat
end
