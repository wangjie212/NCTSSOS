mutable struct ncupop_type
    supp # support data
    coe # coefficient data
    partition # the first 'partition' variables commute with the remaining variables
    comm_var # the first 'comm_var' variables commute each other
    constraint # nothing or "projection" or "unipotent"
    basis # basis
    obj # "eigen" or "trace"
    ksupp # extended support at the k-th step
    moment # moment matrix
    GramMat # Gram matrix
end

mutable struct cosmo_para
    eps_abs::Float64
    eps_rel::Float64
    max_iter::Int64
end

cosmo_para() = cosmo_para(1e-8, 1e-8, 1e4)

"""
    opt,data = nctssos_first(f, x; newton=true, reducebasis=true, TS="block", solver="Mosek", writetofile=false, 
    obj="eigen", partition=0, comm_var=0, constraint=nothing, merge=false, md=3, 
    solve=true, Gram=false, QUIET=false)

Compute the first step of the NCTSSOS hierarchy for unconstrained noncommutative polynomial optimization.
If `newton=true`, then compute a monomial basis by the Newton chip method.
If `reducebasis=true`, then remove monomials from the monomial basis by diagonal inconsistency.
If `TS="block"`, use maximal chordal extensions; if `TS="MD"`, use approximately smallest chordal
extensions. If obj="eigen", minimize the eigenvalue; if obj="trace", then minimize the trace.
If `merge=true`, perform the PSD block merging. Return the optimum and other auxiliary data.

# Arguments
- `f`: the objective function for unconstrained noncommutative polynomial optimization
- `x`: the set of noncommuting variables
- `partition`: the first 'partition' variables commutes with the remaining variables
- `comm_var`: the first 'comm_var' variables commutes each other
- `constraint`: nothing or "projection" or "unipotent"
- `md`: the tunable parameter for merging blocks
"""
function nctssos_first(f::Poly{T}, x; order=0, newton=true, reducebasis=true, monosquare=true,
    TS="block", obj="eigen", partition=0, comm_var=0, constraint=nothing, merge=false, md=3, solve=true, Gram=false, solver="Mosek", writetofile=false, 
    QUIET=false, cosmo_setting=cosmo_para()) where T<:Number
    n,supp,coe = poly_info(f, x)
    opt,data = nctssos_first(supp, coe, n, order=order, newton=newton, reducebasis=reducebasis, monosquare=monosquare, TS=TS, obj=obj,
    partition=partition, comm_var=comm_var, constraint=constraint, merge=merge, md=md, solve=solve, Gram=Gram, solver=solver, writetofile=writetofile, 
    QUIET=QUIET, cosmo_setting=cosmo_setting)
    return opt,data
end

function nctssos_first(supp::Vector{Vector{UInt16}}, coe, n::Int; order=0, newton=true, reducebasis=true, monosquare=true, TS="block",
    obj="eigen", partition=0, comm_var=0, constraint=nothing, merge=false, md=3, solve=true, Gram=false, solver="Mosek", writetofile=false, 
    QUIET=false, cosmo_setting=cosmo_para())
    println("********************************** NCTSSOS **********************************")
    println("NCTSSOS is launching...")
    if order == 0
        order = Int(maximum(length.(supp))/2)
    end
    if obj == "trace"
        supp,coe = cyclic_canon(supp, coe)
    else
        supp,coe = sym_canon(supp, coe)
    end
    if newton == true && constraint === nothing
        if obj == "trace"
            basis = newton_cyclic(supp, n, order)
        else
            basis = newton_ncbasis(supp)
        end
    else
        basis = get_ncbasis(n, order, binary=constraint!==nothing)
    end
    ind = [_comm(item, partition, comm_var) == item for item in basis]
    basis = basis[ind]
    ksupp = copy(supp)
    if monosquare == true && TS != false
        if obj == "trace"
            append!(ksupp, [_cyclic_canon([item[end:-1:1]; item]) for item in basis])
        else
            append!(ksupp, [[item[end:-1:1]; item] for item in basis])
        end
        ksupp = _comm.(ksupp, partition, comm_var)
        if constraint !== nothing
            constraint_reduce!.(ksupp, constraint = constraint)
        end
        sort!(ksupp)
        unique!(ksupp)
    end
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    blocks,cl,blocksize = get_blocks(ksupp, basis, TS=TS, QUIET=QUIET, merge=merge, md=md, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    if reducebasis == true && obj == "eigen" && constraint === nothing
        psupp = copy(supp)
        psupp = psupp[is_sym.(psupp)]
        push!(psupp, UInt16[])
        basis,flag = reducebasis!(psupp, basis, blocks, cl, blocksize)
        if flag == 1
            ksupp = copy(supp)
            if monosquare == true
                if obj == "trace"
                    append!(ksupp, [_cyclic_canon([item[end:-1:1]; item]) for item in basis])
                else
                    append!(ksupp, [[item[end:-1:1]; item] for item in basis])
                end
                ksupp = _comm.(ksupp, partition, comm_var)
                sort!(ksupp)
                unique!(ksupp)
            end
            blocks,cl,blocksize = get_blocks(ksupp, basis, TS=TS, QUIET=QUIET, merge=merge, md=md, obj=obj)
        end
    end
    if QUIET == false
        mb = maximum(blocksize)
        println("Obtained the block structure.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,GramMat = solvesdp(supp, coe, basis, blocks, cl, blocksize, QUIET=QUIET, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint, solve=solve,
    Gram=Gram, solver=solver, writetofile=writetofile, cosmo_setting=cosmo_setting)
    data = ncupop_type(supp, coe, partition, comm_var, constraint, basis, obj, ksupp, moment, GramMat)
    return opt,data
end

"""
    opt,data = nctssos_higher!(data, TS="block", merge=false, md=3, solve=true, Gram=false, QUIET=false)

Compute higher steps of the NCTSSOS hierarchy.
Return the optimum and other auxiliary data.
"""
function nctssos_higher!(data::ncupop_type; TS="block", merge=false, md=3, solve=true, solver="Mosek", writetofile=false, Gram=false, QUIET=false)
    supp = data.supp
    basis = data.basis
    coe = data.coe
    partition = data.partition
    comm_var = data.comm_var
    constraint = data.constraint
    obj = data.obj
    ksupp = data.ksupp
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    oblocksize = deepcopy(data.blocksize)
    blocks,cl,blocksize = get_blocks(ksupp, basis, TS=TS, QUIET=QUIET, merge=merge, md=md, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    if data.blocksize == oblocksize
        println("No higher TS step of the NCTSSOS hierarchy!")
        opt = nothing
    else   
        if QUIET == false
            mb = maximum(blocksize)
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,GramMat = solvesdp(supp, coe, basis, blocks, cl, blocksize, QUIET=QUIET, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint, solve=solve,
        Gram=Gram, solver=solver, writetofile=writetofile, cosmo_setting=cosmo_setting)
        data.moment=moment
        data.GramMat = GramMat
        data.ksupp = ksupp
    end
    return opt,data
end

function cc(a::Vector{UInt16}, n::Int)
    ua = unique(a)
    ca = zeros(UInt8, n)
    for item in ua
        ca[item] = count(x->isequal(item, x), a)
    end
    return ca
end

function remove(csupp, dw, n)
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    t = @variable(model)
    alpha = @variable(model, [1:n])
    @constraint(model, [i=1:length(csupp)], alpha'*(csupp[i].-2*dw)<=t)
    @objective(model, Min, t)
    optimize!(model)
    if objective_value(model) >= 0
        return true
    else
        return false
    end
end

function newton_cyclic(supp, n, d)
    pbasis = get_basis(n,d)
    basis = [UInt16[]]
    csupp = cc.(supp, n)
    pushfirst!(csupp, zeros(UInt8, n))
    sort!(csupp)
    unique!(csupp)
    for i = 2:size(pbasis,2)
        if remove(csupp, pbasis[:,i], n)
            append!(basis, permutation(pbasis[:,i]))
        end
    end
    sort!(basis)
    return basis
end

function newton_ncbasis(supp)
    nbasis = [UInt16[]]
    for bi in supp
        if iseven(length(bi))
            k = Int(length(bi)/2)
            w = bi[end-k+1:end]
            if isequal(w, bi[k:-1:1])
                for j = 1:k
                    push!(nbasis, w[end-j+1:end])
                end
            end
        end
    end
    sort!(nbasis)
    unique!(nbasis)
    return nbasis
end

function get_graph(ksupp, basis; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    lb = length(basis)
    G = SimpleGraph(lb)
    lksupp = length(ksupp)
    for i = 1:lb, j = i+1:lb
        bi = [basis[i][end:-1:1]; basis[j]]
        bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint, comm_var=comm_var)
        if bfind(ksupp, lksupp, bi) !== nothing
           add_edge!(G, i, j)
        end
    end
    return G
end

function get_blocks(ksupp, basis; TS="block", obj="eigen", partition=0, comm_var=0, constraint=nothing, QUIET=true, merge=false, md=3)
    if TS == false
        blocksize = [length(basis)]
        blocks = [[i for i=1:length(basis)]]
        cl = 1
    else
        G = get_graph(ksupp, basis, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
        if TS == "block"
            blocks = connected_components(G)
            blocksize = length.(blocks)
            cl = length(blocksize)
        else
            blocks,cl,blocksize = chordal_cliques!(G, method=TS)
            if merge == true
                blocks,cl,blocksize = clique_merge!(blocks, d=md, QUIET=true)
            end
        end
    end
    if QUIET == false
        sb = sort(unique(blocksize), rev=true)
        numb = [sum(blocksize.== i) for i in sb]
        println("-----------------------------------------------------------------------------")
        println("The sizes of PSD blocks:\n$sb\n$numb")
        println("-----------------------------------------------------------------------------")
    end
    return blocks,cl,blocksize
end

function solvesdp(supp, coe, basis, blocks, cl, blocksize; QUIET=true, obj="eigen", partition=0, comm_var=0, constraint=nothing, solve=true, Gram=false,
    solver="Mosek", writetofile=false, cosmo_setting=cosmo_para())
    ksupp = Vector{Vector{UInt16}}(undef, Int(sum(Int.(blocksize).^2+blocksize)/2))
    k = 1
    for i = 1:cl, j = 1:blocksize[i], r = j:blocksize[i]
        @inbounds bi = [basis[blocks[i][j]][end:-1:1]; basis[blocks[i][r]]]
        @inbounds ksupp[k] = bi
        k += 1
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
        cons = [AffExpr(0) for i=1:lksupp]
        pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl)
        for i = 1:cl
            bs = blocksize[i]
            if bs == 1
               @inbounds pos[i] = @variable(model, lower_bound=0)
               @inbounds bi = [basis[blocks[i][1]][end:-1:1]; basis[blocks[i][1]]]
               bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
               Locb = bfind(ksupp,lksupp,bi)
               @inbounds add_to_expression!(cons[Locb], pos[i])
            else
               @inbounds pos[i] = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:blocksize[i], r = j:blocksize[i]
                   @inbounds bi = [basis[blocks[i][j]][end:-1:1]; basis[blocks[i][r]]]
                   bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                   Locb = bfind(ksupp, lksupp, bi)
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
            Locb = bfind(ksupp, lksupp, supp[i])
            if Locb === nothing
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing,nothing
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
            GramMat = [value.(pos[i]) for i = 1:cl]
        end
        dual_var = -dual.(con)
        moment = Vector{Matrix{Float64}}(undef, cl)
        for i = 1:cl
            moment[i] = zeros(blocksize[i],blocksize[i])
            for j = 1:blocksize[i], k = j:blocksize[i]
                @inbounds bi = [basis[blocks[i][j]][end:-1:1]; basis[blocks[i][k]]]
                bi = reduce!(bi, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
                Locb = bfind(ksupp, lksupp, bi)
                moment[i][j,k] = dual_var[Locb]
            end
            moment[i] = Symmetric(moment[i],:U)
        end
    end
    return objv,ksupp,moment,GramMat
end
