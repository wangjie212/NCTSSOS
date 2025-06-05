mutable struct cstateopt_type
    supp # support data
    coe # coefficient data
    scalar # number of scalar variables
    vargroup # variables commute across groups
    constraint # nothing or "projection" or "unipotent"
    ptsupp # real pure state support
    iptsupp # imaginary pure state support
    wbasis # word basis
    tbasis # real state basis
    itbasis # imaginary state basis
    basis # non-state basis
    blocks # block structure
    cl # number of blocks
    blocksize # size of blocks
    ksupp # extended support at the k-th step
    moment # moment matrix
    GramMat # Gram matrix
end

function get_cwbasis(n, d, ptsupp, iptsupp, bsupp)
    ind = [length(item) <= d for item in bsupp]
    basis = bsupp[ind]
    inx = findfirst(i->length(ptsupp[i])>d, 1:length(ptsupp)) - 1
    iinx = findfirst(i->length(iptsupp[i])>d, 1:length(iptsupp)) - 1
    if d == 4
        tbasis = get_tbasis4(ptsupp[1:inx])
        itbasis = get_tbasis4(iptsupp[1:iinx])
    elseif d == 5
        tbasis = get_tbasis5(ptsupp[1:inx])
        itbasis = get_tbasis5(iptsupp[1:iinx])
    else
        tbasis = get_tbasis(inx, d, ptsupp[1:inx])
        itbasis = get_tbasis(iinx, d, iptsupp[1:iinx])
    end
    wbasis = Vector{UInt16}[]
    for i = 1:length(tbasis), j = 1:length(itbasis), k = 1:length(basis)
        if sum(length.(ptsupp[tbasis[i]])) + sum(length.(iptsupp[itbasis[j]])) + length(basis[k]) <= d
            push!(wbasis, [i,j,k])
        end
    end
    return wbasis,tbasis,itbasis,basis
end

function cpstateopt_first(st_supp::Vector{Vector{Vector{Vector{Int}}}}, coe, n, d; vargroup=[n], TS="block", solver="Mosek", writetofile=false,
    QUIET=false, constraint=nothing, solve=true, Gram=false, bilocal=false, cosmo_setting=cosmo_para(), zero_moments=false)
    return cpstateopt_first([st_supp], [coe], n, d, vargroup=vargroup, TS=TS, solver=solver, writetofile=writetofile, QUIET=QUIET,
    constraint=constraint, solve=solve, Gram=Gram, bilocal=bilocal, cosmo_setting=cosmo_setting, zero_moments=zero_moments)
end

function cpstateopt_first(st_supp::Vector{Vector{Vector{Vector{Vector{Int}}}}}, coe, n, d; vargroup=[n], TS="block", solver="Mosek", writetofile=false,
    QUIET=false, constraint=nothing, solve=true, Gram=false, bilocal=false, cosmo_setting=cosmo_para(), zero_moments=false)
    println("********************************** NCTSSOS **********************************")
    println("NCTSSOS is launching...")
    bsupp = get_ncbasis(n, d, binary=constraint!==nothing)
    ind = [iscomm(item, vargroup) for item in bsupp]
    bsupp = bsupp[ind]
    ptsupp = get_ncbasis(vargroup[1], 2d, ind=Vector{UInt16}(1:vargroup[1]), binary=constraint!==nothing)
    l = vargroup[1]
    for i = 2:length(vargroup)
        nptsupp = Vector{UInt16}[]
        temp = get_ncbasis(vargroup[i], 2d, ind=Vector{UInt16}(l+1:l+vargroup[i]), binary=constraint!==nothing)
        for item1 in ptsupp, item2 in temp
            if length(item1) + length(item2) <= 2d
                push!(nptsupp, [item1;item2])
            end
        end
        ptsupp = nptsupp
        l += vargroup[i]
    end
    if bilocal == false
        ind = [sym(item, vargroup)==item for item in ptsupp]
    else
        ind = [isbilocal(item, bilocal) && sym(item, vargroup)==item for item in ptsupp]
    end
    ptsupp = ptsupp[ind]
    ptsupp = ptsupp[2:end]
    if zero_moments == true
        others = [[1], [2], [3], [4], [5], [6], [7], [8], [9],
         [1;5], [1;6], [2;4], [2;6], [3;4],
         [3;5], [4;8], [4;9], [5;7], [5;9], [6;7], [6;8],
         [1;4;7], [1;4;8], [1;4;9], [1;5;7], [1;5;8], [1;6;7], [1;6;9],
         [2;4;7], [2;4;8], [2;5;7], [2;5;8], [2;5;9], [2;6;8], [2;6;9],
         [3;4;7], [3;4;9], [3;5;8], [3;5;9], [3;6;7], [3;6;8], [3;6;9]]
         sort!(others)
         ind = [bfind(others, length(others), item) === nothing for item in ptsupp]
         ptsupp = ptsupp[ind]
    end
    sort!(ptsupp, lt=isless_td)
    ind = [!issym(item, vargroup) for item in ptsupp]
    iptsupp = ptsupp[ind]
    supp = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(st_supp))
    for i = 1:length(st_supp) 
        supp[i] = Vector{Vector{Vector{UInt16}}}(undef, length(st_supp[i]))
        for k = 1:length(st_supp[i])
            supp[i][k] = Vector{Vector{UInt16}}(undef, 2)
            temp1 = UInt16[]
            temp2 = UInt16[]
            for j = 1:length(st_supp[i][k][1])
                ind = bfind(ptsupp, length(ptsupp), st_supp[i][k][1][j], lt=isless_td)
                if ind !== nothing
                    push!(temp1, ind)
                end
            end
            for j = 1:length(st_supp[i][k][2])
                ind = bfind(iptsupp, length(iptsupp), st_supp[i][k][2][j], lt=isless_td)
                if ind !== nothing
                    push!(temp2, ind)
                end
            end
            supp[i][k][1] = sort(temp1)
            supp[i][k][2] = sort(temp2)
        end
    end
    m = length(st_supp) - 1
    dg = [maximum([sum(length.(st_supp[i+1][j][1])+length.(st_supp[i+1][j][2])) for j=1:length(st_supp[i+1])]) for i=1:m]
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    wbasis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    tbasis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    itbasis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    basis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    wbasis[1],tbasis[1],itbasis[1],basis[1] = get_cwbasis(n, d, ptsupp, iptsupp, bsupp)
    ksupp = copy(supp[1])
    for i = 1:m
        wbasis[i+1],tbasis[i+1],itbasis[i+1],basis[i+1] = get_cwbasis(n, d-Int(ceil(dg[i]/2)), ptsupp, iptsupp, bsupp)
        append!(ksupp, supp[i+1])
    end
    sort!(ksupp)
    unique!(ksupp)
    blocks,cl,blocksize = get_blocks(ksupp, ptsupp, iptsupp, wbasis, tbasis, itbasis, basis, supp=supp, vargroup=vargroup, TS=TS, 
    QUIET=QUIET, constraint=constraint, bilocal=bilocal, zero_moments=zero_moments)
    end
    if QUIET == false
        mb = maximum(maximum.(blocksize))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,GramMat = cpstate_SDP(supp, coe, ptsupp, iptsupp, wbasis, tbasis, itbasis, basis, blocks, cl, blocksize, vargroup, solver=solver, QUIET=QUIET,
    constraint=constraint, solve=solve, writetofile=writetofile, Gram=Gram, bilocal=bilocal, cosmo_setting=cosmo_setting, zero_moments=zero_moments)
    data = cstateopt_type(supp, coe, 0, vargroup, constraint, ptsupp, iptsupp, wbasis, tbasis, itbasis, basis, blocks, cl, blocksize, ksupp, moment, GramMat)
    return opt,data
end

function cpstateopt_higher!(data; TS="block", solver="Mosek", writetofile=false, QUIET=false, solve=true, Gram=false, bilocal=false, cosmo_setting=cosmo_para(), zero_moments=false)
    supp = data.supp
    coe = data.coe
    constraint = data.constraint
    vargroup = data.vargroup
    ptsupp = data.ptsupp
    iptsupp = data.iptsupp
    wbasis = data.wbasis
    tbasis = data.tbasis
    itbasis = data.itbasis
    basis = data.basis
    ksupp = data.ksupp
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    oblocksize = deepcopy(data.blocksize)
    time = @elapsed begin
    blocks,cl,blocksize = get_blocks(ksupp, ptsupp, iptsupp, wbasis, tbasis, itbasis, basis, supp=supp, vargroup=vargroup, TS=TS, QUIET=QUIET, 
    constraint=constraint, bilocal=bilocal, zero_moments=zero_moments)
    end
    if blocksize == oblocksize
        println("No higher TS step of the NCTSSOS hierarchy!")
        opt = nothing
    else
        if QUIET == false
            mb = maximum(maximum.(blocksize))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,GramMat = cpstate_SDP(supp, coe, ptsupp, iptsupp, wbasis, tbasis, itbasis, basis, blocks, cl, blocksize, vargroup, solver=solver, QUIET=QUIET,
        constraint=constraint, solve=solve, writetofile=writetofile, Gram=Gram, bilocal=bilocal, cosmo_setting=cosmo_setting, zero_moments=zero_moments)
        data.moment = moment
        data.GramMat = GramMat
        data.ksupp = ksupp
        data.blocks = blocks
        data.cl = cl
        data.blocksize = blocksize
    end
    return opt,data
end

function get_graph(ksupp, ptsupp, iptsupp, wbasis, tbasis, itbasis, basis; vargroup=nothing, constraint=nothing, bilocal=false, zero_moments=false)
    lb = length(wbasis)
    G = SimpleGraph(lb)
    lksupp = length(ksupp)
    for i = 1:lb, j = i+1:lb
        @inbounds bi1 = [[tbasis[wbasis[i][1]]; tbasis[wbasis[j][1]]], [itbasis[wbasis[i][2]]; itbasis[wbasis[j][2]]]]
        @inbounds bi2 = [reverse(basis[wbasis[i][3]]); basis[wbasis[j][3]]]
        if vargroup !== nothing
            res_comm!(bi2, vargroup)
        end
        if constraint !== nothing
            constraint_reduce!(bi2, constraint=constraint)
        end
        if zero_moments == false || mom_zeros(bi2) == false
            if bilocal != false
                wx,wz,flag = bilocal_reduce(bi2, bilocal)
            end
            if zero_moments == false || flag == false || (mom_zeros(wx) == false && mom_zeros(wz) == false)
                temp1,temp2 = cstate_reduce(bi1, bi2, ptsupp, iptsupp, vargroup, bilocal=bilocal)
                if isempty(temp2)
                    bi = temp1[1]
                else
                    bi = [temp1[1]; temp2[1]]
                end
                r = 1
                while r <= length(bi)
                    if bfind(ksupp, lksupp, bi[r]) !== nothing
                        break
                    else
                        r += 1
                    end
                end
                if r <= length(bi)
                    add_edge!(G, i, j)
                end
            end
        end
    end
    return G
end

function get_graph(ksupp, ptsupp, iptsupp, supp, wbasis, tbasis, itbasis, basis; vargroup=nothing, constraint=nothing, bilocal=false)
    lb = length(wbasis)
    G = SimpleGraph(lb)
    lksupp = length(ksupp)
    for i = 1:lb, j = i+1:lb
        r = 1
        while r <= length(supp)
            @inbounds bi1 = [[tbasis[wbasis[i][1]]; tbasis[wbasis[j][1]]; supp[r][1]], [itbasis[wbasis[i][2]]; itbasis[wbasis[j][2]]; supp[r][2]]]
            @inbounds bi2 = [reverse(basis[wbasis[i][3]]); basis[wbasis[j][3]]]
            if vargroup !== nothing
                res_comm!(bi2, vargroup)
            end
            if constraint !== nothing
                constraint_reduce!(bi2, constraint=constraint)
            end
            temp1,temp2 = cstate_reduce(bi1, bi2, ptsupp, iptsupp, vargroup, bilocal=bilocal)
            if isempty(temp2)
                bi = temp1[1]
            else
                bi = [temp1[1]; temp2[1]]
            end
            s = 1
            while s <= length(bi)
                if bfind(ksupp, lksupp, bi[s]) !== nothing
                    break
                else
                    s += 1
                end
            end
            if s <= length(bi)
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

function get_blocks(ksupp, ptsupp, iptsupp, wbasis, tbasis, itbasis, basis; supp=[], vargroup=nothing, TS="block", QUIET=false, 
    constraint=nothing, bilocal=false, zero_moments=false)
    m = length(wbasis) - 1
    blocks = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    blocksize = Vector{Vector{Int}}(undef, m+1)
    cl = Vector{Int}(undef, m+1)
    if TS == false
        blocksize = [[length(wbasis[k])] for k = 1:m+1]
        blocks = [[Vector(1:length(wbasis[k]))] for k = 1:m+1]
        cl = ones(Int, m+1)
    else
        for k = 1:m+1
            if k == 1
                G = get_graph(ksupp, ptsupp, iptsupp, wbasis[1], tbasis[1], itbasis[1], basis[1], vargroup=vargroup, constraint=constraint, bilocal=bilocal, zero_moments=zero_moments)
            else
                G = get_graph(ksupp, ptsupp, iptsupp, supp[k], wbasis[k], tbasis[k], itbasis[k], basis[k], vargroup=vargroup, constraint=constraint, bilocal=bilocal)
            end
            if TS == "block"
                blocks[k] = connected_components(G)
                blocksize[k] = length.(blocks[k])
                cl[k] = length(blocksize[k])
            else
                blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS)
            end
            if QUIET == false
                sb = sort(Int.(unique(blocksize[k])), rev=true)
                numb = [sum(blocksize[k].== i) for i in sb]
                println("-----------------------------------------------------------------------------")
                println("The sizes of PSD blocks for the $k-th SOHS multiplier:\n$sb\n$numb")
                println("-----------------------------------------------------------------------------")
            end
        end
    end
    return blocks,cl,blocksize
end

function cstate_reduce(word1, word2, ptsupp, iptsupp, vargroup; bilocal=false)
    if isempty(word2)
        return [[[sort(word1[1]), sort(word1[2])]], [1]],[]
    elseif bilocal == false
        sw = sym(word2, vargroup)
        ind = UInt32(bfind(ptsupp, length(ptsupp), sw, lt=isless_td))
        if issym(word2, vargroup)
            return [[[sort([word1[1]; ind]), sort(word1[2])]], [1]],[]
        else
            iind = UInt32(bfind(iptsupp, length(iptsupp), sw, lt=isless_td))
            c = 1
            if sw != word2
                c = -1
            end
            return [[[sort([word1[1]; ind]), sort(word1[2])]], [1]], [[[sort(word1[1]), sort([word1[2]; iind])]], [c]]
        end
    else
        wx,wz,flag = bilocal_reduce(word2, bilocal)
        if flag == true
            swx = sym(wx, vargroup)
            swz = sym(wz, vargroup)
            temp1 = bfind(ptsupp, length(ptsupp), swx, lt=isless_td)
            temp2 = bfind(ptsupp, length(ptsupp), swz, lt=isless_td)
            c1 = c2 = 1
            ssx = ssz = false
            if issym(swx, vargroup)
                ssx = true
            else
                temp3 = bfind(iptsupp, length(iptsupp), swx, lt=isless_td)
                if swx != wx
                    c1 = -1
                end
            end
            if issym(swz, vargroup)
                ssz = true
            else
                temp4 = bfind(iptsupp, length(iptsupp), swz, lt=isless_td)
                if swz != wz
                    c2 = -1
                end
            end
            ind = UInt32[temp1; temp2]
            if ssx == true && ssz == true
                return [[[sort([word1[1]; temp1; temp2]), sort(word1[2])]], [1]],[]
            elseif ssx == false && ssz == true
                return [[[sort([word1[1]; temp1; temp2]), sort(word1[2])]], [1]], 
            [[[sort([word1[1]; temp2]), sort([word1[2]; temp3])]], [c1]]
            elseif ssx == true  && ssz == false
                return [[[sort([word1[1]; temp1; temp2]), sort(word1[2])]], [1]], 
            [[[sort([word1[1]; temp1]), sort([word1[2]; temp4])]], [c2]]
            else
                return [[[sort([word1[1]; temp1; temp2]), sort(word1[2])], [sort(word1[1]), sort([word1[2]; temp3; temp4])]], [1; -1*c1*c2]], 
            [[[sort([word1[1]; temp1]), sort([word1[2]; temp4])], [sort([word1[1]; temp2]), sort([word1[2]; temp3])]], [c2; c1]]
            end
        else
            sw = sym(word2, vargroup)
            ind = UInt32(bfind(ptsupp, length(ptsupp), sw, lt=isless_td))
            if issym(word2, vargroup)
                return [[[sort([word1[1]; ind]), sort(word1[2])]], [1]],[]
            else
                iind = UInt32(bfind(iptsupp, length(iptsupp), sw, lt=isless_td))
                c = 1
                if sw != word2
                    c = -1
                end
                return [[[sort([word1[1]; ind]), sort(word1[2])]], [1]], [[[sort(word1[1]), sort([word1[2]; iind])]], [c]]
            end
        end
    end
end

function cpstate_SDP(supp, coe, ptsupp, iptsupp, wbasis, tbasis, itbasis, basis, blocks, cl, blocksize, vargroup; solver="Mosek", writetofile=false, QUIET=false, constraint=nothing,
    solve=true, Gram=false, bilocal=false, cosmo_setting=cosmo_para(), zero_moments=false)
    m = length(supp) - 1
    ksupp = Vector{Vector{UInt32}}[]
    for i = 1:cl[1], j = 1:blocksize[1][i], r = j:blocksize[1][i]
        @inbounds bi1 = [[tbasis[1][wbasis[1][blocks[1][i][j]][1]]; tbasis[1][wbasis[1][blocks[1][i][r]][1]]], [itbasis[1][wbasis[1][blocks[1][i][j]][2]]; itbasis[1][wbasis[1][blocks[1][i][r]][2]]]]
        @inbounds bi2 = [reverse(basis[1][wbasis[1][blocks[1][i][j]][3]]); basis[1][wbasis[1][blocks[1][i][r]][3]]]
        res_comm!(bi2, vargroup)
        if constraint !== nothing
            constraint_reduce!(bi2, constraint=constraint)
        end
        if zero_moments == false || mom_zeros(bi2) == false
            if bilocal != false
                wx,wz,flag = bilocal_reduce(bi2, bilocal)
            end
            if zero_moments == false || flag == false || (mom_zeros(wx) == false && mom_zeros(wz) == false)
                temp1,temp2 = cstate_reduce(bi1, bi2, ptsupp, iptsupp, vargroup, bilocal=bilocal)
                if isempty(temp2)
                    bi = temp1[1]
                else
                    bi = [temp1[1]; temp2[1]]
                end
                append!(ksupp, bi)
            end
        end
    end
    for k = 1:m, i = 1:cl[k+1], j = 1:blocksize[k+1][i], r = j:blocksize[k+1][i], s = 1:length(supp[k+1])
        @inbounds bi1 = [[tbasis[k+1][wbasis[k+1][blocks[k+1][i][j]][1]]; tbasis[k+1][wbasis[k+1][blocks[k+1][i][r]][1]]; supp[k+1][s][1]], [itbasis[k+1][wbasis[k+1][blocks[k+1][i][j]][2]]; itbasis[k+1][wbasis[k+1][blocks[k+1][i][r]][2]]; supp[k+1][s][2]]]
        @inbounds bi2 = [reverse(basis[k+1][wbasis[k+1][blocks[k+1][i][j]][3]]); basis[k+1][wbasis[k+1][blocks[k+1][i][r]][3]]]
        res_comm!(bi2, vargroup)
        if constraint !== nothing
            constraint_reduce!(bi2, constraint=constraint)
        end
        temp1,temp2 = cstate_reduce(bi1, bi2, ptsupp, iptsupp, vargroup, bilocal=bilocal)
        if isempty(temp2)
            bi = temp1[1]
        else
            bi = [temp1[1]; temp2[1]]
        end
        append!(ksupp, bi)
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
        if solver == "COSMO"
            model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter))
        else
            model = Model(optimizer_with_attributes(Mosek.Optimizer))
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = [AffExpr(0) for i=1:lksupp]
        pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1])
        for i = 1:cl[1]
            bs = blocksize[1][i]
            if bs == 1
               @inbounds pos[i] = @variable(model, lower_bound=0)
               @inbounds bi1 = [[tbasis[1][wbasis[1][blocks[1][i][1]][1]]; tbasis[1][wbasis[1][blocks[1][i][1]][1]]], [itbasis[1][wbasis[1][blocks[1][i][1]][2]]; itbasis[1][wbasis[1][blocks[1][i][1]][2]]]]
               @inbounds bi2 = [reverse(basis[1][wbasis[1][blocks[1][i][1]][3]]); basis[1][wbasis[1][blocks[1][i][1]][3]]]
               res_comm!(bi2, vargroup)
               if constraint !== nothing
                   constraint_reduce!(bi2, constraint=constraint)
               end
               if zero_moments == false || mom_zeros(bi2) == false
                   if bilocal != false
                       wx,wz,flag = bilocal_reduce(bi2, bilocal)
                   end
                   if zero_moments == false || flag == false || (mom_zeros(wx) == false && mom_zeros(wz) == false)
                       bi = cstate_reduce(bi1, bi2, ptsupp, iptsupp, vargroup, bilocal=bilocal)[1][1][1]
                       Locb = bfind(ksupp, lksupp, bi)
                       @inbounds add_to_expression!(cons[Locb], pos[i])
                   end
               end
            else
               @inbounds pos[i] = @variable(model, [1:2bs, 1:2bs], PSD)
               for j = 1:blocksize[1][i], r = j:blocksize[1][i]
                   @inbounds bi1 = [[tbasis[1][wbasis[1][blocks[1][i][j]][1]]; tbasis[1][wbasis[1][blocks[1][i][r]][1]]], [itbasis[1][wbasis[1][blocks[1][i][j]][2]]; itbasis[1][wbasis[1][blocks[1][i][r]][2]]]]
                   @inbounds bi2 = [reverse(basis[1][wbasis[1][blocks[1][i][j]][3]]); basis[1][wbasis[1][blocks[1][i][r]][3]]]
                   res_comm!(bi2, vargroup)
                   if constraint !== nothing
                       constraint_reduce!(bi2, constraint=constraint)
                   end
                   if zero_moments == false || mom_zeros(bi2) == false
                       if bilocal != false
                           wx,wz,flag = bilocal_reduce(bi2, bilocal)
                       end
                       if zero_moments == false || flag == false || (mom_zeros(wx) == false && mom_zeros(wz) == false)
                           temp1,temp2 = cstate_reduce(bi1, bi2, ptsupp, iptsupp, vargroup, bilocal=bilocal)
                           for t = 1:length(temp1[1])
                               Locb = bfind(ksupp, lksupp, temp1[1][t])
                               if j != r
                                   @inbounds add_to_expression!(cons[Locb], 2*temp1[2][t], pos[i][j,r]+pos[i][j+bs,r+bs])
                               else
                                   @inbounds add_to_expression!(cons[Locb], temp1[2][t], pos[i][j,r]+pos[i][j+bs,r+bs])
                               end
                           end
                           if j != r && !isempty(temp2)
                               for t = 1:length(temp2[1])
                                   Locb = bfind(ksupp, lksupp, temp2[1][t])
                                   @inbounds add_to_expression!(cons[Locb], -2*temp2[2][t], pos[i][j,r+bs]-pos[i][r,j+bs])
                               end
                           end
                       end
                   end
               end
            end
        end
        for k = 1:m, i = 1:cl[k+1]
            bs = blocksize[k+1][i]
            if bs == 1
                @inbounds gpos = @variable(model, lower_bound=0)
                for s = 1:length(supp[k+1])
                    @inbounds bi1 = [[tbasis[k+1][wbasis[k+1][blocks[k+1][i][1]][1]]; tbasis[k+1][wbasis[k+1][blocks[k+1][i][1]][1]]; supp[k+1][s][1]], [itbasis[k+1][wbasis[k+1][blocks[k+1][i][1]][2]]; itbasis[k+1][wbasis[k+1][blocks[k+1][i][1]][2]]; supp[k+1][s][2]]]
                    @inbounds bi2 = [reverse(basis[k+1][wbasis[k+1][blocks[k+1][i][1]][3]]); basis[k+1][wbasis[k+1][blocks[k+1][i][1]][3]]]
                    res_comm!(bi2, vargroup)
                    if constraint !== nothing
                        constraint_reduce!(bi2, constraint=constraint)
                    end
                    bi = cstate_reduce(bi1, bi2, ptsupp, iptsupp, vargroup, bilocal=bilocal)[1][1][1]
                    Locb = bfind(ksupp, lksupp, bi)
                    @inbounds add_to_expression!(cons[Locb], coe[k+1][s], gpos)
                end
            else
                @inbounds gpos = @variable(model, [1:2bs, 1:2bs], PSD)
                for j = 1:blocksize[k+1][i], r = j:blocksize[k+1][i], s = 1:length(supp[k+1])
                    @inbounds bi1 = [[tbasis[k+1][wbasis[k+1][blocks[k+1][i][j]][1]]; tbasis[k+1][wbasis[k+1][blocks[k+1][i][r]][1]]; supp[k+1][s]], [itbasis[k+1][wbasis[k+1][blocks[k+1][i][j]][2]]; itbasis[k+1][wbasis[k+1][blocks[k+1][i][r]][2]]; supp[k+1][s]]]
                    @inbounds bi2 = [reverse(basis[k+1][wbasis[k+1][blocks[k+1][i][j]][3]]); basis[k+1][wbasis[k+1][blocks[k+1][i][r]][3]]]
                    res_comm!(bi2, vargroup)
                    if constraint !== nothing
                        constraint_reduce!(bi2, constraint=constraint)
                    end
                    temp1,temp2 = cstate_reduce(bi1, bi2, ptsupp, iptsupp, vargroup, bilocal=bilocal)
                    for t = 1:length(temp1[1])
                        Locb = bfind(ksupp, lksupp, temp1[1][t])
                        @inbounds add_to_expression!(cons[Locb], temp1[2][t]*coe[k+1][s], gpos[j,r]+gpos[j+bs,r+bs])
                    end
                    if j != r && !isempty(temp2)
                        for t = 1:length(temp2[1])
                            Locb = bfind(ksupp, lksupp, temp2[1][t])
                            @inbounds add_to_expression!(icons[Locb], temp2[2][t]*coe[k+1][s], gpos[j,r+bs]-gpos[r,j+bs])
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
        @constraint(model, con[i=1:lksupp], cons[i]==bc[i])
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
        if Gram == true
            GramMat[i] = Vector{ComplexF64}(undef, cl[1])
            for i = 1:cl[1]
                bs = blocksize[1][i]
                temp = value.(pos[1][i][1:bs,bs+1:2bs])
                GramMat[i] = value.(pos[1][i][1:bs,1:bs]+pos[1][i][bs+1:2bs,bs+1:2bs]) + im*(temp-temp')
            end
        end
        var = -dual.(con)
        moment = Vector{Matrix{ComplexF64}}(undef, cl[1])
        for i = 1:cl[1]
            rtemp = zeros(Float64, blocksize[1][i], blocksize[1][i])
            itemp = zeros(Float64, blocksize[1][i], blocksize[1][i])
            for j = 1:blocksize[1][i], k = j:blocksize[1][i]
                @inbounds bi1 = [[tbasis[1][wbasis[1][blocks[1][i][j]][1]]; tbasis[1][wbasis[1][blocks[1][i][k]][1]]], [itbasis[1][wbasis[1][blocks[1][i][j]][2]]; itbasis[1][wbasis[1][blocks[1][i][k]][2]]]]
                @inbounds bi2 = [reverse(basis[1][wbasis[1][blocks[1][i][j]][3]]); basis[1][wbasis[1][blocks[1][i][k]][3]]]
                res_comm!(bi2, vargroup)
                if constraint !== nothing
                    constraint_reduce!(bi2, constraint=constraint)
                end
                if zero_moments == false || mom_zeros(bi2) == false
                    if bilocal != false
                        wx,wz,flag = bilocal_reduce(bi2, bilocal)
                    end
                    if zero_moments == false || flag == false || (mom_zeros(wx) == false && mom_zeros(wz) == false)
                        temp1,temp2 = cstate_reduce(bi1, bi2, ptsupp, iptsupp, vargroup, bilocal=bilocal)
                        for t = 1:length(temp1[1])
                            Locb = bfind(ksupp, lksupp, temp1[1][t])
                            rtemp[j,k] += temp1[2][t]*var[Locb]
                        end
                        if !isempty(temp2)
                            for t = 1:length(temp2[1])
                                Locb = bfind(ksupp, lksupp, temp2[1][t])
                                itemp[j,k] += temp2[2][t]*var[Locb]
                            end
                        else
                            itemp[j,k] = 0
                        end
                    else
                        rtemp[j,k] = 0
                        itemp[j,k] = 0
                    end
                else
                    rtemp[j,k] = 0
                    itemp[j,k] = 0
                end
            end
            rtemp = (rtemp + rtemp')/2
            itemp = (itemp - itemp')/2
            moment[i] = rtemp - itemp*im
        end
    end
    return objv,ksupp,moment,GramMat
end