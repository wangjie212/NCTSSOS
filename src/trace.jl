mutable struct traceopt_type
    supp # support data
    coe # coefficient data
    numeq # number of equality constraints
    constraint # nothing or "projection" or "nilpotent"
    ptsupp # pure trace support
    wbasis # word basis
    tbasis # trace basis
    basis # non-trace basis
    blocks # the block structure
    cl # the number of blocks
    blocksize # the size of blocks
    ksupp # extending support at the k-th step
    sb # sizes of different blocks
    numb # numbers of different blocks
    moment # moment matrix
    GramMat # Gram matrix
end

mutable struct mixword
    ncword # non-commutative part of a mix word
    cword # commutative part of a mix word
end

function get_tbasis(n, d, ptsupp)
    basis = Vector{UInt16}[[]]
    i = 0
    temp = UInt16[]
    while i < d+1
        if sum(temp) == n*i
            temp = ones(UInt16, i+1)
            if i < d && sum(length.(ptsupp[temp])) <= d
                push!(basis, temp)
            end
            i += 1
        else
            temp2 = copy(temp)
            j = temp[1]
            ind = findfirst(x->temp[x]!=j, 1:length(temp))
            if ind === nothing
                ind = length(temp)+1
            end
            if j != 1
                temp2[1:ind-2] = ones(UInt16, ind-2)
            end
            temp2[ind-1] = j+1
            temp = temp2
            if sum(length.(ptsupp[temp])) <= d
                push!(basis, temp)
            end
        end
    end
    return basis
end

function get_tbasis4(ptsupp)
    d1 = findfirst(i->length(ptsupp[i])>1, 1:length(ptsupp)) - 1
    d2 = findfirst(i->length(ptsupp[i])>2, 1:length(ptsupp)) - 1 - d1
    d3 = findfirst(i->length(ptsupp[i])>3, 1:length(ptsupp)) - 1 - d1 - d2
    basis = get_tbasis(d1+d2, 4, ptsupp[1:d1+d2])
    for i = d1+d2+1:d1+d2+d3
        push!(basis, [i])
        for j = 1:d1
            push!(basis, [j;i])
        end
    end
    for i = d1+d2+d3+1:length(ptsupp)
        push!(basis, [i])
    end
    return basis
end

function get_tbasis5(ptsupp)
    d1 = findfirst(i->length(ptsupp[i])>1, 1:length(ptsupp)) - 1
    d2 = findfirst(i->length(ptsupp[i])>2, 1:length(ptsupp)) - 1 - d1
    d3 = findfirst(i->length(ptsupp[i])>3, 1:length(ptsupp)) - 1 - d1 - d2
    d4 = findfirst(i->length(ptsupp[i])>4, 1:length(ptsupp)) - 1 - d1 - d2 - d3
    basis = get_tbasis(d1+d2, 5, ptsupp[1:d1+d2])
    for i = d1+d2+1:d1+d2+d3
        push!(basis, [i])
        for j = 1:d1+d2
            push!(basis, [j;i])
        end
        for j = 1:d1, k = j:d1
            push!(basis, [j;k;i])
        end
    end
    for i = d1+d2+d3+1:d1+d2+d3+d4
        push!(basis, [i])
        for j = 1:d1
            push!(basis, [j;i])
        end
    end
    for i = d1+d2+d3+d4+1:length(ptsupp)
        push!(basis, [i])
    end
    return basis
end

function get_wbasis(n, d, ptsupp, bsupp; scalar=0)
    ind = [length(item) <= d for item in bsupp]
    basis = bsupp[ind]
    if scalar > 0
        ind = [i == 1 || maximum(item) <= n-scalar for item in basis]
        basis = basis[ind]
    end
    inx = findfirst(i->length(ptsupp[i])>d, 1:length(ptsupp)) - 1
    if d == 4
        tbasis = get_tbasis4(ptsupp[1:inx])
    elseif d == 5
        tbasis = get_tbasis5(ptsupp[1:inx])
    else
        tbasis = get_tbasis(inx, d, ptsupp[1:inx])
    end
    wbasis = Vector{UInt16}[]
    for i = 1:length(tbasis), j = 1:length(basis)
        if sum(length.(ptsupp[tbasis[i]])) + length(basis[j]) <= d
            push!(wbasis, [i,j])
        end
    end
    return wbasis,tbasis,basis
end

function isless_td(a, b)
    if length(a) < length(b)
        return true
    elseif length(a) > length(b)
        return false
    else
        return a < b
    end
end

function sym_cyclic(word)
    return min(_cyclic_canon(word), _cyclic_canon(reverse(word)))
end

function traceopt_first(tr_supp::Vector{Vector{Union{Vector{Vector{Int}}, mixword}}}, coe, n, d; numeq=0, TS="block", monosquare=false, QUIET=false, constraint=nothing, solve=true, Gram=false,
    solver="Mosek", cosmo_setting=cosmo_para())
    println("********************************** NCTSSOS **********************************")
    println("Version 0.2.0, developed by Jie Wang, 2020--2023")
    println("NCTSSOS is launching...")
    bsupp = get_ncbasis(n, d, binary=constraint!==nothing)
    ptsupp = get_ncbasis(n, 2d, binary=constraint!==nothing)
    if constraint !== nothing
        ind = [length(item) <= 1 || (item[1] != item[end] && sym_cyclic(item)==item) for item in ptsupp]
    else
        ind = [length(item) <= 1 || sym_cyclic(item)==item for item in ptsupp]
    end
    ptsupp = ptsupp[ind]
    ptsupp = ptsupp[2:end]
    sort!(ptsupp, lt=isless_td)
    m = length(tr_supp) - 1
    supp = Vector{Vector{Vector{Union{Vector{UInt16}, UInt16}}}}(undef, length(tr_supp))
    supp[1] = [sort(UInt16[bfind(ptsupp, length(ptsupp), tr_supp[1][i][j], lt=isless_td) for j=1:length(tr_supp[1][i])])  for i=1:length(tr_supp[1])]
    for k = 1:m
        supp[k+1] = Vector{Union{Vector{UInt16}, UInt16}}(undef, length(tr_supp[k+1]))
        for i = 1:length(tr_supp[k+1])
            supp[k+1][i] = Vector{Union{Vector{UInt16}, UInt16}}(undef, 2)
            if tr_supp[k+1][i].ncword != []
                supp[k+1][i][1] = UInt16(bfind(ptsupp, length(ptsupp), tr_supp[k+1][i].ncword, lt=isless_td))
            else
                supp[k+1][i][1] = UInt16(0)
            end
            supp[k+1][i][2] = sort(UInt16[bfind(ptsupp, length(ptsupp), tr_supp[k+1][i].cword[j], lt=isless_td) for j=1:length(tr_supp[k+1][i].cword)])
        end
    end
    dg = [maximum([length(tr_supp[i+1][j].ncword) + sum(length.(tr_supp[i+1][j].cword), init=0) for j=1:length(tr_supp[i+1])]) for i=1:m]
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    wbasis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    tbasis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    basis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    wbasis[1],tbasis[1],basis[1] = get_wbasis(n, d, ptsupp, bsupp)
    ksupp = copy(supp[1])
    for i = 1:m
        wbasis[i+1],tbasis[i+1],basis[i+1] = get_wbasis(n, d-Int(ceil(dg[i]/2)), ptsupp, bsupp)
        for j = 1:length(supp[i+1])
            bi = sort([supp[i+1][j][1]; supp[i+1][j][2]])
            if bi != []
                push!(ksupp, bi)
            end
        end
    end
    if monosquare == true
        for i = 1:length(wbasis[1])
            bi1 = sort([tbasis[1][wbasis[1][i][1]]; tbasis[1][wbasis[1][i][1]]])
            bi2 = [reverse(basis[1][wbasis[1][i][2]]); basis[1][wbasis[1][i][2]]]
            constraint_reduce!(bi2, constraint=constraint)
            bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
            push!(ksupp, bi)
        end
    end
    sort!(ksupp)
    unique!(ksupp)
    blocks,cl,blocksize,sb,numb,_ = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, supp=supp, TS=TS, QUIET=QUIET, constraint=constraint)
    end
    if QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,GramMat = ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, numeq=numeq, QUIET=QUIET, constraint=constraint,
    solve=solve, Gram=Gram, solver=solver, cosmo_setting=cosmo_setting)
    data = traceopt_type(supp, coe, numeq, constraint, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, ksupp, sb, numb, moment, GramMat)
    return opt,data
end

function traceopt_higher!(data; TS="block", QUIET=false, solve=true, solver="Mosek", Gram=false, cosmo_setting=cosmo_para())
    supp = data.supp
    coe = data.coe
    numeq = data.numeq
    constraint = data.constraint
    ptsupp = data.ptsupp
    wbasis = data.wbasis
    tbasis = data.tbasis
    basis = data.basis
    ksupp = data.ksupp
    sb = data.sb
    numb = data.numb
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,sb,numb,status = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, supp=supp, sb=sb, numb=numb, TS=TS, QUIET=QUIET, constraint=constraint)
    end
    opt = moment = nothing
    if status == 1
        if QUIET == false
            mb = maximum(maximum.(sb))
            println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,GramMat = ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, numeq=numeq, QUIET=QUIET, constraint=constraint,
        solve=solve, Gram=Gram, solver=solver, cosmo_setting=cosmo_setting)
    end
    data.ksupp = ksupp
    data.blocks = blocks
    data.cl = cl
    data.blocksize = blocksize
    data.sb = sb
    data.numb = numb
    data.moment = moment
    data.GramMat = GramMat
    return opt,data
end

function ptraceopt_first(tr_supp::Vector{Vector{Vector{Int}}}, coe, n, d; numeq=0, TS="block", monosquare=false, solver="Mosek",
    QUIET=false, constraint=nothing, solve=true, Gram=false, cosmo_setting=cosmo_para())
    return ptraceopt_first([tr_supp], [coe], n, d, numeq=numeq, TS=TS, monosquare=monosquare, solver=solver, QUIET=QUIET,
    constraint=constraint, solve=solve, Gram=Gram, cosmo_setting=cosmo_setting)
end

function ptraceopt_first(tr_supp::Vector{Vector{Vector{Vector{Int}}}}, coe, n, d; numeq=0, TS="block", monosquare=false, QUIET=false, constraint=nothing, solve=true, Gram=false,
    solver="Mosek", cosmo_setting=cosmo_para())
    println("********************************** NCTSSOS **********************************")
    println("Version 0.2.0, developed by Jie Wang, 2020--2023")
    println("NCTSSOS is launching...")
    bsupp = get_ncbasis(n, d, binary=constraint!==nothing)
    ptsupp = get_ncbasis(n, 2d, binary=constraint!==nothing)
    if constraint !== nothing
        ind = [length(item) <= 1 || (item[1] != item[end] && sym_cyclic(item)==item) for item in ptsupp]
    else
        ind = [length(item) <= 1 || sym_cyclic(item)==item for item in ptsupp]
    end
    ptsupp = ptsupp[ind]
    ptsupp = ptsupp[2:end]
    sort!(ptsupp, lt=isless_td)
    m = length(tr_supp) - 1
    supp = Vector{Vector{Vector{Union{Vector{UInt16}, UInt16}}}}(undef, length(tr_supp))
    supp[1] = [sort(UInt16[bfind(ptsupp, length(ptsupp), tr_supp[1][i][j], lt=isless_td) for j=1:length(tr_supp[1][i])])  for i=1:length(tr_supp[1])]
    for k = 1:m
        supp[k+1] = Vector{Union{Vector{UInt16}, UInt16}}(undef, length(tr_supp[k+1]))
        for i = 1:length(tr_supp[k+1])
            supp[k+1][i] = Vector{Union{Vector{UInt16}, UInt16}}(undef, 2)
            if tr_supp[k+1][i].ncword != []
                supp[k+1][i][1] = UInt16(bfind(ptsupp, length(ptsupp), tr_supp[k+1][i].ncword, lt=isless_td))
            else
                supp[k+1][i][1] = UInt16(0)
            end
            supp[k+1][i][2] = sort(UInt16[bfind(ptsupp, length(ptsupp), tr_supp[k+1][i].cword[j], lt=isless_td) for j=1:length(tr_supp[k+1][i].cword)])
        end
    end
    dg = [maximum([length(tr_supp[i+1][j].ncword) + sum(length.(tr_supp[i+1][j].cword), init=0) for j=1:length(tr_supp[i+1])]) for i=1:m]
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    wbasis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    tbasis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    basis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    wbasis[1],tbasis[1],basis[1] = get_wbasis(n, d, ptsupp, bsupp)
    ksupp = copy(supp[1])
    for i = 1:m
        wbasis[i+1],tbasis[i+1],basis[i+1] = get_wbasis(n, d-Int(ceil(dg[i]/2)), ptsupp, bsupp)
        for j = 1:length(supp[i+1])
            bi = sort([supp[i+1][j][1]; supp[i+1][j][2]])
            if bi != []
                push!(ksupp, bi)
            end
        end
    end
    if monosquare == true
        for i = 1:length(wbasis[1])
            bi1 = sort([tbasis[1][wbasis[1][i][1]]; tbasis[1][wbasis[1][i][1]]])
            bi2 = [reverse(basis[1][wbasis[1][i][2]]); basis[1][wbasis[1][i][2]]]
            constraint_reduce!(bi2, constraint=constraint)
            bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
            push!(ksupp, bi)
        end
    end
    sort!(ksupp)
    unique!(ksupp)
    blocks,cl,blocksize,sb,numb,_ = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, supp=supp, TS=TS, QUIET=QUIET, constraint=constraint)
    end
    if QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,GramMat = ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, numeq=numeq, QUIET=QUIET, constraint=constraint,
    solve=solve, Gram=Gram, solver=solver, cosmo_setting=cosmo_setting)
    data = traceopt_type(supp, coe, numeq, constraint, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, ksupp, sb, numb, moment, GramMat)
    return opt,data
end

function ptraceopt_higher!(data; TS="block", QUIET=false, solve=true, solver="Mosek", Gram=false, cosmo_setting=cosmo_para())
    supp = data.supp
    coe = data.coe
    numeq = data.numeq
    constraint = data.constraint
    ptsupp = data.ptsupp
    wbasis = data.wbasis
    tbasis = data.tbasis
    basis = data.basis
    ksupp = data.ksupp
    sb = data.sb
    numb = data.numb
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,sb,numb,status = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, supp=supp, sb=sb, numb=numb, TS=TS, QUIET=QUIET, constraint=constraint)
    end
    opt = moment = nothing
    if status == 1
        if QUIET == false
            mb = maximum(maximum.(sb))
            println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,GramMat = ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, numeq=numeq, QUIET=QUIET, constraint=constraint,
        solve=solve, Gram=Gram, solver=solver, cosmo_setting=cosmo_setting)
        data.moment = moment
        data.GramMat = GramMat
    end
    data.ksupp = ksupp
    data.blocks = blocks
    data.cl = cl
    data.blocksize = blocksize
    data.sb = sb
    data.numb = numb
    return opt,data
end

function ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize; numeq=0, QUIET=false, constraint=nothing, solve=true,
    Gram=false, solver="Mosek", cosmo_setting=cosmo_para())
    m = length(supp) - 1
    # ksupp = Vector{Vector{UInt16}}(undef, Int(sum(Int.(blocksize).^2+blocksize)/2))
    # k = 1
    ksupp = Vector{UInt32}[]
    for i = 1:cl[1], j = 1:blocksize[1][i], r = j:blocksize[1][i]
        @inbounds bi1 = sort([tbasis[1][wbasis[1][blocks[1][i][j]][1]]; tbasis[1][wbasis[1][blocks[1][i][r]][1]]])
        @inbounds bi2 = [reverse(basis[1][wbasis[1][blocks[1][i][j]][2]]); basis[1][wbasis[1][blocks[1][i][r]][2]]]
        if constraint !== nothing
            constraint_reduce!(bi2, constraint=constraint)
        end
        bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
        push!(ksupp, bi)
        # @inbounds ksupp[k] = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
        # k += 1
    end
    for k = 1:m, i = 1:cl[k+1], j = 1:blocksize[k+1][i], r = j:blocksize[k+1][i], s = 1:length(supp[k+1])
        @inbounds bi1 = sort([tbasis[k+1][wbasis[k+1][blocks[k+1][i][j]][1]]; supp[k+1][s][2]; tbasis[k+1][wbasis[k+1][blocks[k+1][i][r]][1]]])
        temp = supp[k+1][s][1] != 0 ? ptsupp[supp[k+1][s][1]] : UInt16[]
        @inbounds bi2 = [reverse(basis[k+1][wbasis[k+1][blocks[k+1][i][j]][2]]); temp; basis[k+1][wbasis[k+1][blocks[k+1][i][r]][2]]]
        if constraint !== nothing
            constraint_reduce!(bi2, constraint=constraint)
        end
        bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
        push!(ksupp, bi)
    end
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
        cons = [AffExpr(0) for i=1:lksupp]
        pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1])
        for i = 1:cl[1]
            bs = blocksize[1][i]
            if bs == 1
               @inbounds pos[i] = @variable(model, lower_bound=0)
               @inbounds bi1 = sort([tbasis[1][wbasis[1][blocks[1][i][1]][1]]; tbasis[1][wbasis[1][blocks[1][i][1]][1]]])
               @inbounds bi2 = [reverse(basis[1][wbasis[1][blocks[1][i][1]][2]]); basis[1][wbasis[1][blocks[1][i][1]][2]]]
               if constraint !== nothing
                   constraint_reduce!(bi2, constraint=constraint)
               end
               bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
               Locb = bfind(ksupp, lksupp, bi)
               @inbounds add_to_expression!(cons[Locb], pos[i])
            else
               @inbounds pos[i] = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:blocksize[1][i], r = j:blocksize[1][i]
                   @inbounds bi1 = sort([tbasis[1][wbasis[1][blocks[1][i][j]][1]]; tbasis[1][wbasis[1][blocks[1][i][r]][1]]])
                   @inbounds bi2 = [reverse(basis[1][wbasis[1][blocks[1][i][j]][2]]); basis[1][wbasis[1][blocks[1][i][r]][2]]]
                   if constraint !== nothing
                       constraint_reduce!(bi2, constraint=constraint)
                   end
                   bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
                   Locb = bfind(ksupp, lksupp, bi)
                   if j == r
                       @inbounds add_to_expression!(cons[Locb], pos[i][j,r])
                   else
                       @inbounds add_to_expression!(cons[Locb], 2, pos[i][j,r])
                   end
               end
            end
        end
        for k = 1:m, i = 1:cl[k+1]
            bs = blocksize[k+1][i]
            if bs == 1
                if k <= m-numeq
                    gpos = @variable(model, lower_bound=0)
                else
                    gpos = @variable(model)
                end
                for s = 1:length(supp[k+1])
                    @inbounds bi1 = sort([tbasis[k+1][wbasis[k+1][blocks[k+1][i][1]][1]]; supp[k+1][s][2]; tbasis[k+1][wbasis[k+1][blocks[k+1][i][1]][1]]])
                    temp = supp[k+1][s][1] != 0 ? ptsupp[supp[k+1][s][1]] : UInt16[]
                    @inbounds bi2 = [reverse(basis[k+1][wbasis[k+1][blocks[k+1][i][1]][2]]); temp; basis[k+1][wbasis[k+1][blocks[k+1][i][1]][2]]]
                    if constraint !== nothing
                        constraint_reduce!(bi2, constraint=constraint)
                    end
                    bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
                    Locb = bfind(ksupp, lksupp, bi)
                    @inbounds add_to_expression!(cons[Locb], coe[k+1][s], gpos)
                end
            else
                if k <= m-numeq
                    @inbounds gpos = @variable(model, [1:bs, 1:bs], PSD)
                else
                    @inbounds gpos = @variable(model, [1:bs, 1:bs], Symmetric)
                end
                for j = 1:blocksize[k+1][i], r = j:blocksize[k+1][i], s = 1:length(supp[k+1])
                    @inbounds bi1 = sort([tbasis[k+1][wbasis[k+1][blocks[k+1][i][j]][1]]; supp[k+1][s][2]; tbasis[k+1][wbasis[k+1][blocks[k+1][i][r]][1]]])
                    temp = supp[k+1][s][1] != 0 ? ptsupp[supp[k+1][s][1]] : UInt16[]
                    @inbounds bi2 = [reverse(basis[k+1][wbasis[k+1][blocks[k+1][i][j]][2]]); temp; basis[k+1][wbasis[k+1][blocks[k+1][i][r]][2]]]
                    if constraint !== nothing
                        constraint_reduce!(bi2, constraint=constraint)
                    end
                    bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
                    Locb = bfind(ksupp, lksupp, bi)
                    if j == r
                        @inbounds add_to_expression!(cons[Locb], coe[k+1][s], gpos[j,r])
                    else
                        @inbounds add_to_expression!(cons[Locb], 2*coe[k+1][s], gpos[j,r])
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
            GramMat = [value.(pos[i]) for i = 1:cl[1]]
        end
        dual_var = -dual.(con)
        moment = Vector{Matrix{Float64}}(undef, cl[1])
        for i = 1:cl[1]
            moment[i] = zeros(blocksize[1][i],blocksize[1][i])
            for j = 1:blocksize[1][i], k = j:blocksize[1][i]
                @inbounds bi1 = sort([tbasis[1][wbasis[1][blocks[1][i][j]][1]]; tbasis[1][wbasis[1][blocks[1][i][k]][1]]])
                @inbounds bi2 = [reverse(basis[1][wbasis[1][blocks[1][i][j]][2]]); basis[1][wbasis[1][blocks[1][i][k]][2]]]
                if constraint !== nothing
                    constraint_reduce!(bi2, constraint=constraint)
                end
                bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
                Locb = bfind(ksupp, lksupp, bi)
                moment[i][j,k] = dual_var[Locb]
            end
            moment[i] = Symmetric(moment[i],:U)
        end
    end
    return objv,ksupp,moment,GramMat
end

function trace_reduce(word1, word2, ptsupp; constraint=nothing)
    if constraint !== nothing
        while length(word2) > 2 && word2[1] == word2[end]
            if constraint == "nilpotent"
                word2 = word2[2:end-1]
            elseif constraint == "projection"
                word2 = word2[1:end-1]
            end
        end
    end
    if isempty(word2)
        ind = UInt16[]
    else
        ind = UInt16(bfind(ptsupp, length(ptsupp), sym_cyclic(word2), lt=isless_td))
    end
    return sort([word1; ind])
end

function constraint_reduce!(word; constraint="nilpotent")
    i = 1
    while i < length(word)
        if word[i] == word[i+1]
            deleteat!(word, i)
            if constraint == "nilpotent"
                deleteat!(word, i)
            end
            i = 1
        else
            i += 1
        end
    end
    return word
end

function get_ncgraph(ksupp, ptsupp, wbasis, tbasis, basis; vargroup=nothing, constraint=nothing, type="trace", bilocal=false)
    lb = length(wbasis)
    G = SimpleGraph(lb)
    lksupp = length(ksupp)
    for i = 1:lb, j = i+1:lb
        @inbounds bi1 = sort([tbasis[wbasis[i][1]]; tbasis[wbasis[j][1]]])
        @inbounds bi2 = [reverse(basis[wbasis[i][2]]); basis[wbasis[j][2]]]
        if vargroup !== nothing
            res_comm!(bi2, vargroup)
        end
        if constraint !== nothing
            constraint_reduce!(bi2, constraint=constraint)
        end
        if type == "trace"
            bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
        end
        if bilocal == false || bilocal_zeros(bi2) == false
            if bilocal == true
                wx,wz,flag = bilocal_reduce(bi2)
            end
            if bilocal == false || flag == false || (bilocal_zeros(wx) == false && bilocal_zeros(wz) == false)
                if type == "state"
                    bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
                end
                if bfind(ksupp, lksupp, bi) !== nothing
                    add_edge!(G, i, j)
                end
            end
        end
    end
    return G
end

function get_nccgraph(ksupp, ptsupp, supp, wbasis, tbasis, basis; vargroup=nothing, constraint=nothing, type="trace", bilocal=false)
    lb = length(wbasis)
    G = SimpleGraph(lb)
    lksupp = length(ksupp)
    for i = 1:lb, j = i+1:lb
        r = 1
        while r <= length(supp)
            if type == "trace"
                @inbounds bi1 = sort([tbasis[wbasis[i][1]]; supp[r][2]; tbasis[wbasis[j][1]]])
                temp = supp[r][1] != 0 ? ptsupp[supp[r][1]] : UInt16[]
                @inbounds bi2 = [reverse(basis[wbasis[i][2]]); temp; basis[wbasis[j][2]]]
                if constraint !== nothing
                    constraint_reduce!(bi2, constraint=constraint)
                end
                bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
            else
                @inbounds bi1 = sort([tbasis[wbasis[i][1]]; supp[r]; tbasis[wbasis[j][1]]])
                @inbounds bi2 = [reverse(basis[wbasis[i][2]]); basis[wbasis[j][2]]]
                if vargroup !== nothing
                    res_comm!(bi2, vargroup)
                end
                if constraint !== nothing
                    constraint_reduce!(bi2, constraint=constraint)
                end
                bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
            end
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

function get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis; supp=[], vargroup=nothing, sb=[], numb=[], TS="block", QUIET=false, constraint=nothing, type="trace", bilocal=false)
    m = length(wbasis) - 1
    blocks = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    blocksize = Vector{Vector{Int}}(undef, m+1)
    cl = Vector{Int}(undef, m+1)
    if TS == false
        blocksize = [[length(wbasis[k])] for k = 1:m+1]
        blocks = [[Vector(1:length(wbasis[k]))] for k = 1:m+1]
        cl = ones(Int, m+1)
        status = 1
        nsb = Int.(blocksize[1])
        nnumb = [1]
    else
        G = get_ncgraph(ksupp, ptsupp, wbasis[1], tbasis[1], basis[1], vargroup=vargroup, constraint=constraint, type=type, bilocal=bilocal)
        if TS == "block"
            blocks[1] = connected_components(G)
            blocksize[1] = length.(blocks[1])
            cl[1] = length(blocksize[1])
        else
            blocks[1],cl[1],blocksize[1] = chordal_cliques!(G, method=TS)
        end
        nsb = sort(unique(blocksize[1]), rev=true)
        nnumb = [sum(blocksize[1].== i) for i in nsb]
        if isempty(sb) || nsb != sb || nnumb != numb
            status = 1
            if QUIET == false
                println("-----------------------------------------------------------------------------")
                println("The sizes of PSD blocks:\n$nsb\n$nnumb")
                println("-----------------------------------------------------------------------------")
            end
            for k = 1:m
                G = get_nccgraph(ksupp, ptsupp, supp[k+1], wbasis[k+1], tbasis[k+1], basis[k+1], vargroup=vargroup, constraint=constraint, type=type, bilocal=bilocal)
                if TS == "block"
                    blocks[k+1] = connected_components(G)
                    blocksize[k+1] = length.(blocks[k+1])
                    cl[k+1] = length(blocksize[k+1])
                else
                    blocks[k+1],cl[k+1],blocksize[k+1] = chordal_cliques!(G, method=TS, minimize=false)
                end
            end
        else
            status = 0
            println("No higher TS step of the NCTSSOS hierarchy!")
        end
    end
    return blocks,cl,blocksize,nsb,nnumb,status
end

function Werner_witness_first(dY, sigma, n, d; TS="block", monosquare=false, QUIET=false, solve=true, solver="Mosek", cosmo_setting=cosmo_para())
    println("********************************** NCTSSOS **********************************")
    println("Version 0.2.0, developed by Jie Wang, 2020--2023")
    println("NCTSSOS is launching...")
    bsupp = get_ncbasis(n, 2d)
    ind = [findfirst(j -> bsupp[i][j] == bsupp[i][j+1], 1:length(bsupp[i])-1) === nothing for i = 1:length(bsupp)]
    bsupp = bsupp[ind]
    ind = [length(bsupp[i]) <= 1 || (bsupp[i][1] != bsupp[i][end] && sym_cyclic(bsupp[i])==bsupp[i]) for i=1:length(bsupp)]
    ptsupp = bsupp[ind]
    ptsupp = ptsupp[2:end]
    sort!(ptsupp, lt=isless_td)
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    wbasis,tbasis,basis = get_wbasis(n, d, ptsupp, bsupp)
    htrace = generate_htrace(n)
    supp = Vector{Vector{UInt16}}(undef, length(htrace))
    for i = 1:length(htrace)
        supp[i] = sort([bfind(ptsupp, length(ptsupp), sym_cyclic(htrace[i][j]), lt=isless_td) for j=1:length(htrace[i])])
    end
    if monosquare == true
        for i = 1:length(wbasis)
            bi1 = sort([tbasis[wbasis[i][1]]; tbasis[wbasis[i][1]]])
            bi2 = [reverse(basis[wbasis[i][2]]); basis[wbasis[i][2]]]
            constraint_reduce!(bi2, constraint=constraint)
            bi = trace_reduce(bi1, bi2, ptsupp, constraint=constraint)
            push!(supp, bi)
        end
    end
    sort!(supp)
    unique!(supp)
    blocks,cl,blocksize,sb,numb,_ = get_ncblocks(supp, ptsupp, wbasis, tbasis, basis, TS=TS, QUIET=QUIET, constraint="projection")
    end
    if QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
    end
    opt,ksupp = Werner_SDP(dY, sigma, htrace, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, QUIET=QUIET, solve=solve, solver=solver)
    data = traceopt_type(htrace, dY, sigma, ptsupp, wbasis, tbasis, basis, ksupp, sb, numb, nothing)
    return opt,data
end

function Werner_witness_higher!(data; TS="block", QUIET=false, solve=true, solver="Mosek", cosmo_setting=cosmo_para())
    htrace = data.supp
    dY = data.coe
    sigma = data.constraint
    ptsupp = data.ptsupp
    wbasis = data.wbasis
    tbasis = data.tbasis
    basis = data.basis
    ksupp = data.ksupp
    sb = data.sb
    numb = data.numb
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    blocks,cl,blocksize,sb,numb,status = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, sb=sb, numb=numb, TS=TS, QUIET=QUIET, constraint="projection")
    opt = nothing
    if status == 1
        if QUIET == false
            mb = maximum(maximum.(sb))
            println("Obtained the block structure. The maximal size of blocks is $mb.")
        end
        opt,ksupp = Werner_SDP(dY, sigma, htrace, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, QUIET=QUIET, solve=solve, solver=solver)
    end
    data.ksupp = ksupp
    data.sb = sb
    data.numb = numb
    return opt,data
end

function Werner_SDP(dY, sigma, htrace, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize; QUIET=false, solve=true, solver="Mosek", cosmo_setting=cosmo_para())
    ksupp = Vector{Vector{UInt16}}(undef, Int(sum(Int.(blocksize).^2+blocksize)/2))
    k = 1
    for i = 1:cl, j = 1:blocksize[i], r = j:blocksize[i]
        @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][r]][1]]])
        @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][r]][2]]]
        constraint_reduce!(bi2, constraint="projection")
        @inbounds ksupp[k] = trace_reduce(bi1, bi2, ptsupp, constraint="projection")
        k += 1
    end
    sort!(ksupp)
    unique!(ksupp)
    lksupp = length(ksupp)
    if QUIET == false
        println("There are $lksupp affine constraints.")
    end
    objv = nothing
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
        cons = [AffExpr(0) for i=1:lksupp]
        for i = 1:cl
            bs = blocksize[i]
            if bs == 1
               @inbounds pos = @variable(model, lower_bound=0)
               @inbounds bi1 = sort([tbasis[wbasis[blocks[i][1]][1]]; tbasis[wbasis[blocks[i][1]][1]]])
               @inbounds bi2 = [reverse(basis[wbasis[blocks[i][1]][2]]); basis[wbasis[blocks[i][1]][2]]]
               constraint_reduce!(bi2, constraint="projection")
               bi = trace_reduce(bi1, bi2, ptsupp, constraint="projection")
               Locb = bfind(ksupp, lksupp, bi)
               @inbounds add_to_expression!(cons[Locb], pos)
            else
               @inbounds pos = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:blocksize[i], r = j:blocksize[i]
                   @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][r]][1]]])
                   @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][r]][2]]]
                   constraint_reduce!(bi2, constraint="projection")
                   bi = trace_reduce(bi1, bi2, ptsupp, constraint="projection")
                   Locb = bfind(ksupp, lksupp, bi)
                   if j == r
                       @inbounds add_to_expression!(cons[Locb], pos[j,r])
                   else
                       @inbounds add_to_expression!(cons[Locb], 2, pos[j,r])
                   end
               end
            end
        end
        bc = [AffExpr(0) for i=1:lksupp]
        coe = @variable(model, [1:length(htrace)])
        tcons = AffExpr(0)
        for i = 1:length(htrace)
            for k = 1:length(dY)
                tcons += sigma[k]*prod(tr(prod(dY[k][htrace[i][j]])) for j=1:length(htrace[i]))*coe[i]
            end
            temp = sort([bfind(ptsupp, length(ptsupp), sym_cyclic(htrace[i][j]), lt=isless_td) for j=1:length(htrace[i])])
            Locb = bfind(ksupp, lksupp, temp)
            if Locb === nothing
               @error "The monomial basis is not enough!"
               return nothing,nothing
            else
               bc[Locb] += coe[i]
            end
        end
        @constraint(model, tcons==-16)
        @variable(model, epsilon)
        cons[1] -= epsilon
        @constraint(model, cons.==bc)
        @objective(model, Min, epsilon)
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
    end
    return objv,ksupp
end

function generate_htrace(n)
    return _generate_htrace(UInt16[i for i=1:n])
end

function _generate_htrace(var)
    if isempty(var)
         htrace = [Vector{UInt16}[]]
    else
        htrace = Vector{Vector{UInt16}}[]
        for i = 1:length(var)
            sset = _selete(var[2:end], i-1)
            pushfirst!.(sset, var[1])
            for mem in sset
                cbasis = _cyclic_basis(mem)
                for emem in cbasis
                    sub_htrace = _generate_htrace(UInt16.(setdiff(var, mem)))
                    for j = 1:length(sub_htrace)
                        pushfirst!(sub_htrace[j], emem)
                    end
                    append!(htrace, sub_htrace)
                end
            end
        end
    end
    return htrace
end

function _selete(var, d)
    if d > 0
        set = Vector{UInt16}[]
        for i = 1:length(var)-d+1
            iset = _selete(var[i+1:end], d-1)
            pushfirst!.(iset, var[i])
            append!(set, iset)
        end
    else
        set = [UInt16[]]
    end
    return set
end

function _cyclic_basis(var)
    basis = _permutation(var, ones(length(var)))
    return unique(_cyclic_canon.(basis))
end
