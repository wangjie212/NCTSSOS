mutable struct traceopt_type
    supp # support data
    coe # coefficient data
    constraint # "projection" or "unipotent"
    ptsupp # pure trace support
    wbasis # word basis
    tbasis # trace basis
    basis # non-trace basis
    ksupp # extending support at the k-th step
    sb # sizes of different blocks
    numb # numbers of different blocks
    moment # moment matrix
    GramMat # Gram matrix
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

function ptraceopt_first(tr_supp, coe, n, d; TS="block", monosquare=false, QUIET=false, constraint="unipotent", solve=true, Gram=false, solver="Mosek")
    println("********************************** NCTSSOS **********************************")
    println("Version 0.2.0, developed by Jie Wang, 2020--2022")
    println("NCTSSOS is launching...")
    bsupp = get_ncbasis(n, d, binary=true)
    ptsupp = get_ncbasis(n, 2d, binary=true)
    ind = [length(item) <= 1 || (item[1] != item[end] && sym_cyclic(item)==item) for item in ptsupp]
    ptsupp = ptsupp[ind]
    ptsupp = ptsupp[2:end]
    sort!(ptsupp, lt=isless_td)
    supp = Vector{Vector{UInt16}}(undef, length(tr_supp))
    for i = 1:length(tr_supp)
        supp[i] = sort([ncbfind(ptsupp, length(ptsupp), tr_supp[i][j], lt=isless_td) for j=1:length(tr_supp[i])])
    end
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    wbasis,tbasis,basis = get_wbasis(n, d, ptsupp, bsupp)
    ksupp = copy(supp)
    if monosquare == true
        for i = 1:length(wbasis)
            bi1 = sort([tbasis[wbasis[i][1]]; tbasis[wbasis[i][1]]])
            bi2 = [reverse(basis[wbasis[i][2]]); basis[wbasis[i][2]]]
            constraint_reduce!(bi2, constraint=constraint)
            while length(bi2) > 2 && bi2[1] == bi2[end]
                if constraint == "unipotent"
                    bi2 = bi2[2:end-1]
                elseif constraint == "projection"
                    bi2 = bi2[1:end-1]
                end
            end
            bi = trace_reduce(bi1, bi2, ptsupp)
            push!(ksupp, bi)
        end
    end
    sort!(ksupp)
    unique!(ksupp)
    blocks,cl,blocksize,sb,numb,_ = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, TS=TS, QUIET=QUIET, constraint=constraint)
    end
    if QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,GramMat = ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, QUIET=QUIET, constraint=constraint, solve=solve, Gram=Gram, solver=solver)
    data = traceopt_type(supp, coe, constraint, ptsupp, wbasis, tbasis, basis, ksupp, sb, numb, moment, GramMat)
    return opt,data
end

function ptraceopt_higher!(data; TS="block", QUIET=false, solve=true, solver="Mosek", Gram=false)
    supp = data.supp
    coe = data.coe
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
    blocks,cl,blocksize,sb,numb,status = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, sb=sb, numb=numb, TS=TS, QUIET=QUIET, constraint=constraint)
    end
    opt = moment = nothing
    if status == 1
        if QUIET == false
            mb = maximum(maximum.(sb))
            println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,GramMat = ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, QUIET=QUIET, constraint=constraint, solve=solve, Gram=Gram, solver=solver)
    end
    data.ksupp = ksupp
    data.sb = sb
    data.numb = numb
    data.moment = moment
    data.GramMat = GramMat
    return opt,data
end

function ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize; QUIET=false, constraint="unipotent", solve=true, Gram=false, solver="Mosek")
    ksupp = Vector{Vector{UInt16}}(undef, Int(sum(Int.(blocksize).^2+blocksize)/2))
    k = 1
    for i = 1:cl, j = 1:blocksize[i], r = j:blocksize[i]
        @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][r]][1]]])
        @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][r]][2]]]
        constraint_reduce!(bi2, constraint=constraint)
        while length(bi2) > 2 && bi2[1] == bi2[end]
            if constraint == "unipotent"
                bi2 = bi2[2:end-1]
            elseif constraint == "projection"
                bi2 = bi2[1:end-1]
            end
        end
        @inbounds ksupp[k] = trace_reduce(bi1, bi2, ptsupp)
        k += 1
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
            model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => 1e-4, "eps_rel" => 1e-4, "max_iter" => 10000))
        else
            @error "The solver is currently not supported!"
            return nothing,nothing,nothing,nothing
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons = [AffExpr(0) for i=1:lksupp]
        pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl)
        for i = 1:cl
            bs = blocksize[i]
            if bs == 1
               @inbounds pos[i] = @variable(model, lower_bound=0)
               @inbounds bi1 = sort([tbasis[wbasis[blocks[i][1]][1]]; tbasis[wbasis[blocks[i][1]][1]]])
               @inbounds bi2 = [reverse(basis[wbasis[blocks[i][1]][2]]); basis[wbasis[blocks[i][1]][2]]]
               constraint_reduce!(bi2, constraint=constraint)
               while length(bi2) > 2 && bi2[1] == bi2[end]
                   if constraint == "unipotent"
                       bi2 = bi2[2:end-1]
                   elseif constraint == "projection"
                       bi2 = bi2[1:end-1]
                   end
               end
               bi = trace_reduce(bi1, bi2, ptsupp)
               Locb = ncbfind(ksupp, lksupp, bi)
               @inbounds add_to_expression!(cons[Locb], pos[i])
            else
               @inbounds pos[i] = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:blocksize[i], r = j:blocksize[i]
                   @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][r]][1]]])
                   @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][r]][2]]]
                   constraint_reduce!(bi2, constraint=constraint)
                   while length(bi2) > 2 && bi2[1] == bi2[end]
                       if constraint == "unipotent"
                            bi2 = bi2[2:end-1]
                       elseif constraint == "projection"
                            bi2 = bi2[1:end-1]
                       end
                   end
                   bi = trace_reduce(bi1, bi2, ptsupp)
                   Locb = ncbfind(ksupp, lksupp, bi)
                   if Locb == 0
                       @error "The word does not exist!"
                       return nothing,nothing
                   end
                   if j == r
                       @inbounds add_to_expression!(cons[Locb], pos[i][j,r])
                   else
                       @inbounds add_to_expression!(cons[Locb], 2, pos[i][j,r])
                   end
               end
            end
        end
        bc = zeros(lksupp)
        for i = 1:length(supp)
            Locb = ncbfind(ksupp, lksupp, supp[i])
            if Locb == 0
               @error "The monomial basis is not enough!"
               return nothing,nothing
            else
               bc[Locb] = coe[i]
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
            GramMat = [value.(pos[i]) for i = 1:cl]
        end
        dual_var = -dual.(con)
        moment = Vector{Matrix{Float64}}(undef, cl)
        for i = 1:cl
            moment[i] = zeros(blocksize[i],blocksize[i])
            for j = 1:blocksize[i], k = j:blocksize[i]
                @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][k]][1]]])
                @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][k]][2]]]
                constraint_reduce!(bi2, constraint=constraint)
                while length(bi2) > 2 && bi2[1] == bi2[end]
                    if constraint == "unipotent"
                        bi2 = bi2[2:end-1]
                    elseif constraint == "projection"
                        bi2 = bi2[1:end-1]
                    end
                end
                bi = trace_reduce(bi1, bi2, ptsupp)
                Locb = ncbfind(ksupp, lksupp, bi)
                moment[i][j,k] = dual_var[Locb]
            end
            moment[i] = Symmetric(moment[i],:U)
        end
    end
    return objv,ksupp,moment,GramMat
end

function trace_reduce(word1, word2, ptsupp)
    if isempty(word2)
        ind = UInt16[]
    else
        ind = UInt16(ncbfind(ptsupp, length(ptsupp), sym_cyclic(word2), lt=isless_td))
    end
    return sort([word1; ind])
end

function constraint_reduce!(word; constraint="unipotent")
    i = 1
    while i < length(word)
        if word[i] == word[i+1]
            deleteat!(word, i)
            if constraint == "unipotent"
                deleteat!(word, i)
            end
            i = 1
        else
            i += 1
        end
    end
    return word
end

function get_ncgraph(ksupp, ptsupp, wbasis, tbasis, basis; vargroup=nothing, constraint="unipotent", type="trace", bilocal=false)
    lb = length(wbasis)
    G = SimpleGraph(lb)
    lksupp = length(ksupp)
    for i = 1:lb, j = i+1:lb
        @inbounds bi1 = sort([tbasis[wbasis[i][1]]; tbasis[wbasis[j][1]]])
        @inbounds bi2 = [reverse(basis[wbasis[i][2]]); basis[wbasis[j][2]]]
        if vargroup !== nothing
            res_comm!(bi2, vargroup)
        end
        constraint_reduce!(bi2, constraint=constraint)
        if type == "trace"
            while length(bi2) > 2 && bi2[1] == bi2[end]
                if constraint == "unipotent"
                    bi2 = bi2[2:end-1]
                elseif constraint == "projection"
                    bi2 = bi2[1:end-1]
                end
            end
            bi = trace_reduce(bi1, bi2, ptsupp)
        else
            bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
        end
        if bilocal == false || bilocal_zeros(bi2) == false
            if bilocal == true
                wx,wz,flag = bilocal_reduce(bi2)
            end
            if bilocal == false || flag == false || (bilocal_zeros(wx) == false && bilocal_zeros(wz) == false)
                if ncbfind(ksupp, lksupp, bi) != 0
                    add_edge!(G, i, j)
                end
            end
        end
    end
    return G
end

function get_nccgraph(ksupp, ptsupp, supp, wbasis, tbasis, basis; vargroup=nothing, constraint="unipotent", type="trace", bilocal=false)
    lb = length(wbasis)
    G = SimpleGraph(lb)
    lksupp = length(ksupp)
    for i = 1:lb, j = i+1:lb
        r = 1
        while r <= length(supp)
            @inbounds bi1 = sort([tbasis[wbasis[i][1]]; supp[r]; tbasis[wbasis[j][1]]])
            @inbounds bi2 = [reverse(basis[wbasis[i][2]]); basis[wbasis[j][2]]]
            if vargroup !== nothing
                res_comm!(bi2, vargroup)
            end
            constraint_reduce!(bi2, constraint=constraint)
            bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
            if ncbfind(ksupp, lksupp, bi) != 0
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

function get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis; supp=[], vargroup=nothing, sb=[], numb=[], TS="block", QUIET=false, constraint="unipotent", type="trace", bilocal=false)
    if type == "state"
        m = length(wbasis) - 1
        blocks = Vector{Vector{Vector{UInt16}}}(undef, m+1)
        blocksize = Vector{Vector{Int}}(undef, m+1)
        cl = Vector{Int}(undef, m+1)
    end
    if TS == false
        if type == "trace"
            blocksize = [length(wbasis)]
            blocks = [[i for i=1:length(wbasis)]]
            cl = 1
        else
            blocksize = [[length(wbasis[k])] for k = 1:m+1]
            blocks = [[Vector(1:length(wbasis[k]))] for k = 1:m+1]
            cl = ones(Int, m+1)
        end
    else
        if type == "trace"
            G = get_ncgraph(ksupp, ptsupp, wbasis, tbasis, basis, vargroup=vargroup, constraint=constraint, type=type)
        else
            G = get_ncgraph(ksupp, ptsupp, wbasis[1], tbasis[1], basis[1], vargroup=vargroup, constraint=constraint, type=type, bilocal=bilocal)
        end
        if TS == "block"
            if type == "trace"
                blocks = connected_components(G)
                blocksize = length.(blocks)
                cl = length(blocksize)
            else
                blocks[1] = connected_components(G)
                blocksize[1] = length.(blocks[1])
                cl[1] = length(blocksize[1])
            end
        else
            if type == "trace"
                blocks,cl,blocksize = chordal_cliques!(G, method=TS)
            else
                blocks[1],cl[1],blocksize[1] = chordal_cliques!(G, method=TS)
            end
        end
    end
    if type == "trace"
        nsb = sort(unique(blocksize), rev=true)
        nnumb = [sum(blocksize.== i) for i in nsb]
    else
        nsb = sort(unique(blocksize[1]), rev=true)
        nnumb = [sum(blocksize[1].== i) for i in nsb]
    end
    if isempty(sb) || nsb!=sb || nnumb!=numb
        status = 1
        if QUIET == false
            println("-----------------------------------------------------------------------------")
            println("The sizes of PSD blocks:\n$nsb\n$nnumb")
            println("-----------------------------------------------------------------------------")
        end
    else
        status = 0
        println("No higher TS step of the NCTSSOS hierarchy!")
    end
    if type == "state"
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
    end
    return blocks,cl,blocksize,nsb,nnumb,status
end

function Werner_witness_first(dY, sigma, n, d; TS="block", monosquare=false, QUIET=false, solve=true, solver="Mosek")
    println("********************************** NCTSSOS **********************************")
    println("Version 0.2.0, developed by Jie Wang, 2020--2022")
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
        supp[i] = sort([ncbfind(ptsupp, length(ptsupp), sym_cyclic(htrace[i][j]), lt=isless_td) for j=1:length(htrace[i])])
    end
    if monosquare == true
        for i = 1:length(wbasis)
            bi1 = sort([tbasis[wbasis[i][1]]; tbasis[wbasis[i][1]]])
            bi2 = [reverse(basis[wbasis[i][2]]); basis[wbasis[i][2]]]
            constraint_reduce!(bi2, constraint="projection")
            while length(bi2) > 2 && bi2[1] == bi2[end]
                bi2 = bi2[1:end-1]
            end
            bi = trace_reduce(bi1, bi2, ptsupp)
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

function Werner_witness_higher!(data; TS="block", QUIET=false, solve=true, solver="Mosek")
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

function Werner_SDP(dY, sigma, htrace, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize; QUIET=false, solve=true, solver="Mosek")
    ksupp = Vector{Vector{UInt16}}(undef, Int(sum(Int.(blocksize).^2+blocksize)/2))
    k = 1
    for i = 1:cl, j = 1:blocksize[i], r = j:blocksize[i]
        @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][r]][1]]])
        @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][r]][2]]]
        constraint_reduce!(bi2, constraint="projection")
        while length(bi2) > 2 && bi2[1] == bi2[end]
            bi2 = bi2[1:end-1]
        end
        @inbounds ksupp[k] = trace_reduce(bi1, bi2, ptsupp)
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
            model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => 1e-4, "eps_rel" => 1e-4, "max_iter" => 10000))
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
               while length(bi2) > 2 && bi2[1] == bi2[end]
                     bi2 = bi2[1:end-1]
               end
               bi = trace_reduce(bi1, bi2, ptsupp)
               Locb = ncbfind(ksupp, lksupp, bi)
               @inbounds add_to_expression!(cons[Locb], pos)
            else
               @inbounds pos = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:blocksize[i], r = j:blocksize[i]
                   @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][r]][1]]])
                   @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][r]][2]]]
                   constraint_reduce!(bi2, constraint="projection")
                   while length(bi2) > 2 && bi2[1] == bi2[end]
                         bi2 = bi2[1:end-1]
                   end
                   bi = trace_reduce(bi1, bi2, ptsupp)
                   Locb = ncbfind(ksupp, lksupp, bi)
                   if Locb == 0
                       @error "The word does not exist!"
                       return nothing,nothing
                   end
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
            temp = sort([ncbfind(ptsupp, length(ptsupp), sym_cyclic(htrace[i][j]), lt=isless_td) for j=1:length(htrace[i])])
            Locb = ncbfind(ksupp, lksupp, temp)
            if Locb == 0
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
