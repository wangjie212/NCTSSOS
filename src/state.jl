mutable struct stateopt_type
    supp # support data
    coe # coefficient data
    scalar # number of scalar variables
    vargroup # variables commute across groups
    constraint # nothing or "projection" or "unipotent"
    ptsupp # pure state support
    wbasis # word basis
    tbasis # state basis
    basis # non-state basis
    blocks # the block structure
    cl # the number of blocks
    blocksize # the size of blocks
    ksupp # extending support at the k-th step
    sb # sizes of different blocks
    numb # numbers of different blocks
    moment # moment matrix
    GramMat # Gram matrix
end

function sym(word, vargroup)
    cword = copy(word)
    ind = [gind(word[i], vargroup) for i = 1:length(word)]
    uind = unique(ind)
    nind = [count(ind .== uind[i]) for i = 1:length(uind)]
    k = 0
    for i = 1:length(uind)
        cword[k+1:k+nind[i]] = reverse(cword[k+1:k+nind[i]])
        k += nind[i]
    end
    return min(word, cword)
end

function iscomm(a, vargroup)
    for i = 1:length(a)-1
        if a[i] > a[i+1] && gind(a[i], vargroup) != gind(a[i+1], vargroup)
            return false
        end
    end
    return true
end

function gind(k, vargroup)
    return findfirst(i -> k <= sum(vargroup[1:i]), 1:length(vargroup))
end

function res_comm!(a, vargroup)
    i = 1
    while i < length(a)
        if a[i] > a[i+1] && gind(a[i], vargroup) != gind(a[i+1], vargroup)
            temp = a[i]
            a[i] = a[i+1]
            a[i+1] = temp
            if i > 1
                i -= 1
            end
        else
            i += 1
        end
    end
    return a
end

function pstateopt_first(st_supp::Vector{Vector{Vector{Int}}}, coe, n, d; scalar=0, vargroup=[n], TS="block", monosquare=false, solver="Mosek",
    QUIET=false, constraint=nothing, solve=true, Gram=false, bilocal=false, cosmo_setting=cosmo_para())
    return pstateopt_first([st_supp], [coe], n, d, scalar=scalar, vargroup=vargroup, TS=TS, monosquare=monosquare, solver=solver, QUIET=QUIET,
    constraint=constraint, solve=solve, Gram=Gram, bilocal=bilocal, cosmo_setting=cosmo_setting)
end

function pstateopt_first(st_supp::Vector{Vector{Vector{Vector{Int}}}}, coe, n, d; scalar=0, vargroup=[n], TS="block", monosquare=false, solver="Mosek",
    QUIET=false, constraint=nothing, solve=true, Gram=false, bilocal=false, cosmo_setting=cosmo_para())
    println("********************************** NCTSSOS **********************************")
    println("Version 0.2.0, developed by Jie Wang, 2020--2023")
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
        ind = [isbilocal(item) && sym(item, vargroup)==item for item in ptsupp]
    end
    ptsupp = ptsupp[ind]
    ptsupp = ptsupp[2:end]
    if bilocal == true
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
    supp = Vector{Vector{Vector{UInt16}}}(undef, length(st_supp))
    for i = 1:length(st_supp)
        supp[i] = [sort([bfind(ptsupp, length(ptsupp), st_supp[i][k][j], lt=isless_td) for j=1:length(st_supp[i][k])]) for k = 1:length(st_supp[i])]
    end
    m = length(st_supp) - 1
    dg = [maximum([sum(length.(st_supp[i+1][j])) for j=1:length(st_supp[i+1])]) for i=1:m]
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    wbasis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    tbasis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    basis = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    wbasis[1],tbasis[1],basis[1] = get_wbasis(n, d, ptsupp, bsupp, scalar=scalar)
    ksupp = copy(supp[1])
    for i = 1:m
        wbasis[i+1],tbasis[i+1],basis[i+1] = get_wbasis(n, d-Int(ceil(dg[i]/2)), ptsupp, bsupp, scalar=scalar)
        append!(ksupp, supp[i+1])
    end
    if monosquare == true
        for i = 1:length(wbasis[1])
            bi1 = sort([tbasis[1][wbasis[1][i][1]]; tbasis[1][wbasis[1][i][1]]])
            bi2 = [reverse(basis[1][wbasis[1][i][2]]); basis[1][wbasis[1][i][2]]]
            res_comm!(bi2, vargroup)
            if constraint !== nothing
                constraint_reduce!(bi2, constraint=constraint)
            end
            bi = state_reduce(bi1, bi2, ptsupp, vargroup)
            push!(ksupp, bi)
        end
    end
    sort!(ksupp)
    unique!(ksupp)
    blocks,cl,blocksize,sb,numb,_ = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, supp=supp, vargroup=vargroup, TS=TS, QUIET=QUIET, constraint=constraint, type="state", bilocal=bilocal)
    end
    if QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,GramMat = pstate_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, vargroup, solver=solver, QUIET=QUIET,
    constraint=constraint, solve=solve, Gram=Gram, bilocal=bilocal, cosmo_setting=cosmo_setting)
    data = stateopt_type(supp, coe, scalar, vargroup, constraint, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, ksupp, sb, numb, moment, GramMat)
    return opt,data
end

function pstateopt_higher!(data; TS="block", solver="Mosek", QUIET=false, solve=true, Gram=false, bilocal=false, cosmo_setting=cosmo_para())
    supp = data.supp
    coe = data.coe
    constraint = data.constraint
    vargroup = data.vargroup
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
    blocks,cl,blocksize,sb,numb,status = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, supp=supp, vargroup=vargroup, sb=sb, numb=numb, TS=TS, QUIET=QUIET, constraint=constraint, type="state", bilocal=bilocal)
    end
    opt = moment = nothing
    if status == 1
        if QUIET == false
            mb = maximum(maximum.(sb))
            println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,GramMat = pstate_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, vargroup, solver=solver, QUIET=QUIET,
        constraint=constraint, solve=solve, Gram=Gram, bilocal=bilocal, cosmo_setting=cosmo_setting)
        data.moment = moment
        data.GramMat = GramMat
    end
    data.ksupp = ksupp
    data.blocks = blocks
    data.cl = cl
    data.blocksize = blocksize
    data.moment = moment
    data.sb = sb
    data.numb = numb
    return opt,data
end

function pstate_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, vargroup; solver="Mosek", QUIET=false, constraint=nothing,
    solve=true, Gram=false, bilocal=false, cosmo_setting=cosmo_para())
    m = length(supp) - 1
    # ksupp = Vector{Vector{UInt32}}(undef, Int(sum(Int.(blocksize[1]).^2+blocksize[1])/2))
    # k = 1
    ksupp = Vector{UInt32}[]
    for i = 1:cl[1], j = 1:blocksize[1][i], r = j:blocksize[1][i]
        @inbounds bi1 = sort([tbasis[1][wbasis[1][blocks[1][i][j]][1]]; tbasis[1][wbasis[1][blocks[1][i][r]][1]]])
        @inbounds bi2 = [reverse(basis[1][wbasis[1][blocks[1][i][j]][2]]); basis[1][wbasis[1][blocks[1][i][r]][2]]]
        res_comm!(bi2, vargroup)
        if constraint !== nothing
            constraint_reduce!(bi2, constraint=constraint)
        end
        if bilocal == false || bilocal_zeros(bi2) == false
            # @inbounds ksupp[k] = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
            # k += 1
            if bilocal == true
                wx,wz,flag = bilocal_reduce(bi2)
            end
            if bilocal == false || flag == false || (bilocal_zeros(wx) == false && bilocal_zeros(wz) == false)
                bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
                push!(ksupp, bi)
            end
        end
    end
    for k = 1:m, i = 1:cl[k+1], j = 1:blocksize[k+1][i], r = j:blocksize[k+1][i], s = 1:length(supp[k+1])
        @inbounds bi1 = sort([tbasis[k+1][wbasis[k+1][blocks[k+1][i][j]][1]]; supp[k+1][s]; tbasis[k+1][wbasis[k+1][blocks[k+1][i][r]][1]]])
        @inbounds bi2 = [reverse(basis[k+1][wbasis[k+1][blocks[k+1][i][j]][2]]); basis[k+1][wbasis[k+1][blocks[k+1][i][r]][2]]]
        res_comm!(bi2, vargroup)
        if constraint !== nothing
            constraint_reduce!(bi2, constraint=constraint)
        end
        bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
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
        if solver == "COSMO"
            model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter))
        else
            model = Model(optimizer_with_attributes(Mosek.Optimizer))
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
               res_comm!(bi2, vargroup)
               if constraint !== nothing
                   constraint_reduce!(bi2, constraint=constraint)
               end
               if bilocal == false || bilocal_zeros(bi2) == false
                   if bilocal == true
                       wx,wz,flag = bilocal_reduce(bi2)
                   end
                   if bilocal == false || flag == false || (bilocal_zeros(wx) == false && bilocal_zeros(wz) == false)
                       bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
                       Locb = bfind(ksupp, lksupp, bi)
                       @inbounds add_to_expression!(cons[Locb], pos[i])
                   end
               end
            else
               @inbounds pos[i] = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:blocksize[1][i], r = j:blocksize[1][i]
                   @inbounds bi1 = sort([tbasis[1][wbasis[1][blocks[1][i][j]][1]]; tbasis[1][wbasis[1][blocks[1][i][r]][1]]])
                   @inbounds bi2 = [reverse(basis[1][wbasis[1][blocks[1][i][j]][2]]); basis[1][wbasis[1][blocks[1][i][r]][2]]]
                   res_comm!(bi2, vargroup)
                   if constraint !== nothing
                       constraint_reduce!(bi2, constraint=constraint)
                   end
                   if bilocal == false || bilocal_zeros(bi2) == false
                    if bilocal == true
                        wx,wz,flag = bilocal_reduce(bi2)
                    end
                    if bilocal == false || flag == false || (bilocal_zeros(wx) == false && bilocal_zeros(wz) == false)
                       bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
                       Locb = bfind(ksupp, lksupp, bi)
                       if j == r
                           @inbounds add_to_expression!(cons[Locb], pos[i][j,r])
                       else
                           @inbounds add_to_expression!(cons[Locb], 2, pos[i][j,r])
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
                    @inbounds bi1 = sort([tbasis[k+1][wbasis[k+1][blocks[k+1][i][1]][1]]; supp[k+1][s]; tbasis[k+1][wbasis[k+1][blocks[k+1][i][1]][1]]])
                    @inbounds bi2 = [reverse(basis[k+1][wbasis[k+1][blocks[k+1][i][1]][2]]); basis[k+1][wbasis[k+1][blocks[k+1][i][1]][2]]]
                    res_comm!(bi2, vargroup)
                    if constraint !== nothing
                        constraint_reduce!(bi2, constraint=constraint)
                    end
                    bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
                    Locb = bfind(ksupp, lksupp, bi)
                    @inbounds add_to_expression!(cons[Locb], coe[k+1][s], gpos)
                end
            else
                @inbounds gpos = @variable(model, [1:bs, 1:bs], PSD)
                for j = 1:blocksize[k+1][i], r = j:blocksize[k+1][i], s = 1:length(supp[k+1])
                    @inbounds bi1 = sort([tbasis[k+1][wbasis[k+1][blocks[k+1][i][j]][1]]; supp[k+1][s]; tbasis[k+1][wbasis[k+1][blocks[k+1][i][r]][1]]])
                    @inbounds bi2 = [reverse(basis[k+1][wbasis[k+1][blocks[k+1][i][j]][2]]); basis[k+1][wbasis[k+1][blocks[k+1][i][r]][2]]]
                    res_comm!(bi2, vargroup)
                    if constraint !== nothing
                        constraint_reduce!(bi2, constraint=constraint)
                    end
                    bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
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
                res_comm!(bi2, vargroup)
                if constraint !== nothing
                    constraint_reduce!(bi2, constraint=constraint)
                end
                if bilocal == false || bilocal_zeros(bi2) == false
                    if bilocal == true
                        wx,wz,flag = bilocal_reduce(bi2)
                    end
                    if bilocal == false || flag == false || (bilocal_zeros(wx) == false && bilocal_zeros(wz) == false)
                        bi = state_reduce(bi1, bi2, ptsupp, vargroup, bilocal=bilocal)
                        Locb = bfind(ksupp, lksupp, bi)
                        moment[i][j,k] = dual_var[Locb]
                    else
                        moment[i][j,k] = 0
                    end
                else
                    moment[i][j,k] = 0
                end
            end
            moment[i] = Symmetric(moment[i],:U)
        end
    end
    return objv,ksupp,moment,GramMat
end

function state_reduce(word1, word2, ptsupp, vargroup; bilocal=false)
    if isempty(word2)
        ind = UInt32[]
    elseif bilocal == false
        ind = UInt32(bfind(ptsupp, length(ptsupp), sym(word2, vargroup), lt=isless_td))
    else
        wx,wz,flag = bilocal_reduce(word2)
        if flag == true
            temp1 = bfind(ptsupp, length(ptsupp), sym(wx, vargroup), lt=isless_td)
            temp2 = bfind(ptsupp, length(ptsupp), sym(wz, vargroup), lt=isless_td)
            ind = UInt32[temp1; temp2]
        else
            ind = UInt32(bfind(ptsupp, length(ptsupp), sym(word2, vargroup), lt=isless_td))
        end
    end
    return sort([word1; ind])
end

function bilocal_zeros(word)
    others = [[1], [2], [3], [4], [5], [6], [7], [8], [9],
         [1;5], [1;6], [2;4], [2;6], [3;4],
         [3;5], [4;8], [4;9], [5;7], [5;9], [6;7], [6;8],
         [1;4;7], [1;4;8], [1;4;9], [1;5;7], [1;5;8], [1;6;7], [1;6;9],
         [2;4;7], [2;4;8], [2;5;7], [2;5;8], [2;5;9], [2;6;8], [2;6;9],
         [3;4;7], [3;4;9], [3;5;8], [3;5;9], [3;6;7], [3;6;8], [3;6;9]]
    sort!(others)
    return bfind(others, length(others), word) !== nothing
end

function bilocal_reduce(word)
    wx = word[word .<= 3]
    wz = word[word .>= 7]
    # wz = wz[wz .<= 9]
    if length(wx) > 0 && length(wz) > 0 && length(wx) + length(wz) == length(word)
        return wx,wz,true
    else
        return wx,wz,false
    end
end

function isbilocal(word)
    wx = word[word .<= 3]
    wz = word[word .>= 7]
    # wz = wz[wz .<= 9]
    if length(wx) > 0 && length(wz) > 0 && length(wx) + length(wz) == length(word)
        return false
    else
        return true
    end
end
