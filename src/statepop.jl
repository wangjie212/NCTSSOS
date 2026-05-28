mutable struct statepop_data
    pop # polynomial optimiztion problem
    numeq # number of equality constraints
    scalar # number of scalar variables
    vargroup # variables commute across groups
    constraint # nothing or "projection" or "unipotent"
    basis # monomial bases
    ksupp # extended support at the k-th TS step
    blocksize # sizes of blocks
    blocks # block structure
    moment # Moment matrix
    GramMat # Gram matrices
end

function statepop(pop::Vector{T}, n, d; numeq=0, scalar=0, vargroup=nothing, TS="block", type="state", monosquare=false, 
    QUIET=false, constraint=nothing, solve=true, Gram=false, writetofile=false, model=nothing, mosek_setting=mosek_para(), 
    bilocal=false, zero_moments=[]) where {T<:poly}
    println("********************************** NCTSSOS **********************************")
    println("NCTSSOS is launching...")
    cost = pop[1]
    pop_cons = [ncstatepoly{Float64}([tuple(Int[], Vector{Int}[])], [1]); pop[2:end]]
    wbasis = get_ncbasis(n-scalar, d, ind=Vector(scalar+1:n), binary=constraint!==nothing)
    wbasis = wbasis[[!iscomm(item, vargroup=vargroup) for item in wbasis]]
    psupp = wbasis[2:end]
    if type == "state"
        psupp = psupp[[sym_inv(item, vargroup=vargroup) == item for item in psupp]]
    else
        if constraint !== nothing
            ind = [length(item) <= 1 || (item[1] != item[end] && sym_cyclic(item, vargroup=vargroup, constraint=constraint) == item) for item in psupp]
        else
            ind = [length(item) <= 1 || sym_cyclic(item, vargroup=vargroup, constraint=constraint)==item for item in psupp]
        end
        psupp = psupp[ind]
    end
    if bilocal != false
        psupp = psupp[[!bilocal_reduce(item, bilocal)[3] for item in psupp]]
    end
    if !isempty(zero_moments)
        psupp = psupp[[bfind(zero_moments, item) === nothing for item in psupp]]
    end
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    basis = Vector{Vector{Union{Vector{Vector{Int}}, Tuple{Vector{Int}, Vector{Vector{Int}}}}}}(undef, length(pop))
    ksupp = copy(cost.supp)
    for i = 1:length(pop)
        if typeof(pop_cons[i]) <: statepoly
            basis[i] = get_statebasis(psupp, d-Int(maxdeg(pop_cons[i])/2), scalar=scalar)
            append!(ksupp, pop_cons[i].supp)
        else
            basis[i] = get_ncstatebasis(wbasis, psupp, d-Int(maxdeg(pop_cons[i])/2), scalar=scalar)
            if type == "state"
                append!(ksupp, [sort([[sym_inv(item[1], vargroup=vargroup)]; item[2]]) for item in pop_cons[i].supp])
            else
                append!(ksupp, [sort([[sym_cyclic(item[1], vargroup=vargroup, constraint=constraint)]; item[2]]) for item in pop_cons[i].supp])
            end
        end
    end
    if monosquare == true
        for item in basis[1]
            bi,c = add(star(item), item, vargroup=vargroup, type=type, scalar=scalar, constraint=constraint, bilocal=bilocal, zero_moments=zero_moments)
            if c != 0
                push!(ksupp, bi)
            end
        end
    end
    sort!(ksupp)
    unique!(ksupp)
    blocks,cl,blocksize = get_blocks(pop_cons, ksupp, basis, vargroup=vargroup, TS=TS, type=type, scalar=scalar, QUIET=QUIET, constraint=constraint, 
    bilocal=bilocal, zero_moments=zero_moments)
    end
    if QUIET == false
        mb = maximum(maximum.(blocksize))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,GramMat = solvesdp(cost, pop_cons, basis, blocks, cl, blocksize, type=type, TS=TS, vargroup=vargroup, numeq=numeq, QUIET=QUIET, constraint=constraint,
    solve=solve, scalar=scalar, model=model, mosek_setting=mosek_setting, Gram=Gram, writetofile=writetofile, bilocal=bilocal, zero_moments=zero_moments)
    data = statepop_data(pop, numeq, scalar, vargroup, constraint, basis, ksupp, blocksize, blocks, moment, GramMat)
    return opt,data
end

function statepop(data; TS="block", type="state", writetofile=false, model=nothing, mosek_setting=mosek_para(), QUIET=false, solve=true, Gram=false, bilocal=false, zero_moments=[])
    cost = data.pop[1]
    pop_cons = [ncstatepoly{Float64}([tuple(Int[], Vector{Int}[])], [1]); data.pop[2:end]]
    numeq = data.numeq
    scalar = data.scalar
    constraint = data.constraint
    vargroup = data.vargroup
    basis = data.basis
    ksupp = data.ksupp
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize = get_blocks(pop_cons, ksupp, basis, vargroup=vargroup, TS=TS, QUIET=QUIET, constraint=constraint, type=type, scalar=scalar, bilocal=bilocal, zero_moments=zero_moments)
    end
    if blocksize == data.blocksize
        println("No higher TS step!")
        opt = nothing
    else
        if QUIET == false
            mb = maximum(maximum.(blocksize))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,GramMat = solvesdp(cost, pop_cons, basis, blocks, cl, blocksize, type=type, TS=TS, vargroup=vargroup, numeq=numeq, writetofile=writetofile, QUIET=QUIET,
        constraint=constraint, scalar=scalar, solve=solve, model=model, mosek_setting=mosek_setting, Gram=Gram, bilocal=bilocal, zero_moments=zero_moments)
        data.moment = moment
        data.GramMat = GramMat
        data.ksupp = ksupp
        data.blocks = blocks
        data.blocksize = blocksize
    end
    return opt,data
end

function solvesdp(cost, pop_cons, basis, blocks, cl, blocksize; type="state", TS="block", vargroup=nothing, numeq=0, scalar=0, writetofile=false, QUIET=false, constraint=nothing,
    solve=true, model=nothing, mosek_setting=mosek_para(), Gram=false, bilocal=false, zero_moments=[])
    ksupp = Vector{Vector{Int}}[]
    for i = 1:cl[1], j = 1:blocksize[1][i], r = j:blocksize[1][i]
        bi,c = add(star(basis[1][blocks[1][i][j]]), basis[1][blocks[1][i][r]], vargroup=vargroup, type=type, scalar=scalar, constraint=constraint, bilocal=bilocal, zero_moments=zero_moments)
        if c != 0
            push!(ksupp, bi)
        end
    end
    if TS != false
        for k = 2:length(pop_cons), i = 1:cl[k], j = 1:blocksize[k][i], r = j:blocksize[k][i], item in pop_cons[k].supp
            bi,c = add(star(basis[k][blocks[k][i][j]]), item, basis[k][blocks[k][i][r]], vargroup=vargroup, type=type, scalar=scalar, constraint=constraint, bilocal=bilocal, zero_moments=zero_moments)
            if c != 0
                push!(ksupp, bi)
            end
        end
    end
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
        cons = [AffExpr(0) for i=1:length(ksupp)]
        pos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(pop_cons))
        for (k, p) in enumerate(pop_cons)
            pos[k] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k])
            for i = 1:cl[k]
                if k <= length(pop_cons) - numeq
                    @inbounds pos[k][i] = @variable(model, [1:blocksize[k][i], 1:blocksize[k][i]], PSD)
                else
                    @inbounds pos[k][i] = @variable(model, [1:blocksize[k][i], 1:blocksize[k][i]], Symmetric)
                end
                for j = 1:blocksize[k][i], r = j:blocksize[k][i], (s, item) in enumerate(p.supp)
                    bi,c = add(star(basis[k][blocks[k][i][j]]), item, basis[k][blocks[k][i][r]], vargroup=vargroup, type=type, scalar=scalar, constraint=constraint, bilocal=bilocal, zero_moments=zero_moments)
                    if c != 0
                        Locb = bfind(ksupp, bi)
                        if j == r
                            @inbounds add_to_expression!(cons[Locb], p.coe[s], pos[k][i][j,r])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2*p.coe[s], pos[k][i][j,r])
                        end
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
                GramMat = [value.(pos[i]) for i = 1:cl[1]]
            end
            dual_var = -dual(con)
            moment = Vector{Matrix{Float64}}(undef, cl[1])
            for i = 1:cl[1]
                moment[i] = zeros(blocksize[1][i], blocksize[1][i])
                for j = 1:blocksize[1][i], k = j:blocksize[1][i]
                    bi,c = add(star(basis[1][blocks[1][i][j]]), basis[1][blocks[1][i][k]], vargroup=vargroup, type=type, scalar=scalar, constraint=constraint, bilocal=bilocal, zero_moments=zero_moments)
                    if c != 0
                        Locb = bfind(ksupp, bi)
                        moment[i][j,k] = dual_var[Locb]   
                    else
                        moment[i][j,k] = 0
                    end
                end
                moment[i] = Symmetric(moment[i], :U)
            end
        end
    end
    return objv,ksupp,moment,GramMat
end

function get_graph(tsupp::Vector{Vector{Vector{Int}}}, supp, basis; vargroup=nothing, constraint=nothing, type="state", scalar=0, bilocal=false, zero_moments=[])
    lb = length(basis)
    G = SimpleGraph(lb)
    for i = 1:lb, j = i+1:lb
        r = 1
        while r <= length(supp)
            bi,c = add(star(basis[i]), supp[r], basis[j], vargroup=vargroup, type=type, scalar=scalar, constraint=constraint, bilocal=bilocal, zero_moments=zero_moments)
            if c == 0 || bfind(tsupp, bi) !== nothing
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

function get_blocks(pop_cons, tsupp::Vector{Vector{Vector{Int}}}, basis; vargroup=nothing, TS="block", QUIET=false, constraint=nothing, type="state", scalar=0, bilocal=false, zero_moments=[])
    blocks = Vector{Vector{Vector{Int}}}(undef, length(pop_cons))
    blocksize = Vector{Vector{Int}}(undef, length(pop_cons))
    cl = Vector{Int}(undef, length(pop_cons))
    if TS == false
        blocksize = [[length(basis[k])] for k = 1:length(pop_cons)]
        blocks = [[Vector(1:length(basis[k]))] for k = 1:length(pop_cons)]
        cl = ones(Int, length(pop_cons))
    else
        for (k, p) in enumerate(pop_cons)
            G = get_graph(tsupp, p.supp, basis[k], vargroup=vargroup, constraint=constraint, type=type, scalar=scalar, bilocal=bilocal, zero_moments=zero_moments)
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
                println("The sizes of PSD blocks for the $k-th multiplier:\n$sb\n$numb")
                println("-----------------------------------------------------------------------------")
            end
        end
    end
    return blocks,cl,blocksize
end

function bilocal_reduce(word, loc)
    wx = word[word .<= loc[1]]
    wz = word[word .>= loc[2]]
    if length(wx) > 0 && length(wz) > 0 && length(wx) + length(wz) == length(word)
        return wx,wz,true
    else
        return wx,wz,false
    end
end

function add(a::Vector{Vector{Int}}...; vargroup=nothing, type="state", scalar=0, constraint=nothing, bilocal=false, zero_moments=[])
    bi = sort(vcat(a...))
    if scalar > 0
        if length(bi) > 2 && bi[1][end] <= scalar && bi[2][end] <= scalar && bi[3][end] <= scalar
            bi = [[sort([bi[1]; bi[2]; bi[3]])]; bi[4:end]]
        elseif length(bi) > 1 && bi[1][end] <= scalar && bi[2][end] <= scalar
            bi = [[sort([bi[1]; bi[2]])]; bi[3:end]]
        end
    end
    return bi,1
end

function add(a::Tuple{Vector{Int}, Vector{Vector{Int}}}...; vargroup=nothing, type="state", scalar=0, constraint=nothing, bilocal=false, zero_moments=[])
    bi1 = vcat([b[1] for b in a]...)
    bi2 = sort(vcat([b[2] for b in a]...))
    if scalar > 0
        if length(bi2) > 2 && bi2[1][end] <= scalar && bi2[2][end] <= scalar && bi2[3][end] <= scalar
            bi2 = [[sort([bi2[1]; bi2[2]; bi2[3]])]; bi2[4:end]]
        elseif length(bi2) > 1 && bi2[1][end] <= scalar && bi2[2][end] <= scalar
            bi2 = [[sort([bi2[1]; bi2[2]])]; bi2[3:end]]
        end
    end
    if vargroup !== nothing
        res_comm!(bi1, vargroup)
    end
    constraint_reduce!(bi1, constraint=constraint)
    if isempty(zero_moments) || bfind(zero_moments, bi1) === nothing
        if !isempty(bi1)
            if bilocal != false
                wx,wz,flag = bilocal_reduce(bi1, bilocal)
                if flag == true
                    if isempty(zero_moments) || (bfind(zero_moments, wx) === nothing && bfind(zero_moments, wz) === nothing)
                        if type == "state"
                            return sort([[sym_inv(wx, vargroup=vargroup), sym_inv(wz, vargroup=vargroup)]; bi2]),1
                        else
                            return sort([[sym_cyclic(wx, vargroup=vargroup, constraint=constraint), sym_cyclic(wz, vargroup=vargroup, constraint=constraint)]; bi2]),1
                        end
                    else
                        return [],0
                    end
                else
                    if type == "state"
                        return sort([[sym_inv(bi1, vargroup=vargroup)]; bi2]),1
                    else
                        return sort([[sym_cyclic(bi1, vargroup=vargroup, constraint=constraint)]; bi2]),1
                    end
                end
            else
                if type == "state"
                    return sort([[sym_inv(bi1, vargroup=vargroup)]; bi2]),1
                else
                    return sort([[sym_cyclic(bi1, vargroup=vargroup, constraint=constraint)]; bi2]),1
                end
            end
        else
            return bi2,1
        end
    else
        return [],0
    end
end
