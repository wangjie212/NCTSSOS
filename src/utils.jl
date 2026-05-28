function get_basis(n::Int, d::Int)
    lb = binomial(n+d, d)
    basis = Vector{Vector{Int}}(undef, lb)
    basis[1] = Int[]
    i = 0
    t = 1
    while i <= d
        t += 1
        if sum(basis[t-1]) == n*i
            if i < d
                basis[t] = ones(Int, i+1)
            end
            i += 1
        else
            basis[t] = copy(basis[t-1])
            ind = findfirst(x -> basis[t][x] != j, 1:length(basis[t]))
            if ind === nothing
                ind = length(basis[t]) + 1
            end
            if j != 1
                basis[t][1:ind-2] = ones(Int, ind-2)
            end
            basis[t][ind-1] = j + 1
        end
    end
    return basis
end

function get_statebasis(psupp, d; scalar=0)
    basis = Vector{Vector{Int}}[[]]
    if d <= 0
        return basis
    end
    i = 0
    temp = Int[]
    while i <= d
        if sum(temp) == length(psupp)*i
            temp = ones(Int, i+1)
            if i < d && sum(length.(psupp[temp])) <= d
                push!(basis, psupp[temp])
            end
            i += 1
        else
            ntemp = copy(temp)
            j = temp[1]
            ind = findfirst(x -> temp[x] != j, 1:length(temp))
            if ind === nothing
                ind = length(temp) + 1
            end
            if j != 1
                ntemp[1:ind-2] = ones(Int, ind-2)
            end
            ntemp[ind-1] = j + 1
            temp = ntemp
            if sum(length.(psupp[temp])) <= d
                push!(basis, psupp[temp])
            end
        end
    end
    if scalar > 0
        nbasis = Vector{Vector{Int}}[]
        sb = get_basis(scalar, d)
        for a in sb, b in basis
            if length(a) + sum(length.(b)) <= d
                if isempty(a)
                    push!(nbasis, b)
                else
                    push!(nbasis, [[a]; b])
                end
            end
        end
        basis = nbasis
    end
    return basis
end

function get_ncstatebasis(wbasis, psupp, d; scalar=0)
    if d == 0
        return [tuple(Int[], Vector{Int}[])]
    end
    basis = Tuple{Vector{Int}, Vector{Vector{Int}}}[]
    sbasis = get_statebasis(psupp, d, scalar=scalar)
    for a in wbasis, b in sbasis
        if length(a) + sum(length.(b)) <= d
            push!(basis, tuple(a, b))
        end
    end
    return basis
end

function _cyclic_canon(a::Vector{Int})
    if isempty(a)
        return a
    else
        return minimum([[a[i+1:length(a)]; a[1:i]] for i = 0:length(a)-1])
    end
end

function cyclic_canon(p::ncpoly{T}) where {T<:Union{Number,AffExpr,GenericAffExpr}}
    nsupp = [min(_cyclic_canon(word), _cyclic_canon(reverse(word))) for word in p.supp]
    sort!(nsupp)
    unique!(nsupp)
    ncoe = zeros(T, length(nsupp))
    for (i, item) in enumerate(p.supp)
        bi = min(_cyclic_canon(item), _cyclic_canon(reverse(item)))
        Locb = bfind(nsupp, bi)
        ncoe[Locb] += p.coe[i]
    end
    return ncpoly(nsupp, ncoe)
end

function _sym_canon(a::Vector{Int})
    i = 1
    while i <= Int(ceil((length(a)-1)/2))
        if a[i] < a[end+1-i]
            return a
        elseif a[i] > a[end+1-i]
            return reverse(a)
        else
            i += 1
        end
    end
    return a
end

function sym_canon(p::ncpoly{T}) where {T<:Union{Number,AffExpr,GenericAffExpr}}
    nsupp = _sym_canon.(p.supp)
    sort!(nsupp)
    unique!(nsupp)
    ncoe = zeros(T, length(nsupp))
    for (i, item) in enumerate(p.supp)
        Locb = bfind(nsupp, _sym_canon(item))
        ncoe[Locb] += p.coe[i]
    end
    return ncpoly(nsupp, ncoe)
end

function cc(a::Vector{Int}, n::Int)
    ua = unique(a)
    ca = zeros(Int, n)
    ca[ua] = [count(a .== x) for x in ua]
    return ca
end

function remove(supp, dw)
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    t = @variable(model)
    alpha = @variable(model, [1:length(dw)])
    @constraint(model, [i=1:length(supp)], alpha'*(supp[i] .- 2dw) <= t)
    @objective(model, Min, t)
    optimize!(model)
    return objective_value(model) >= 0
end

function get_fbasis(n, d)
    lb = binomial(n+d, d)
    basis = zeros(Int, n, lb)
    i = 0
    t = 1
    while i <= d
        t += 1
        if basis[n, t-1] == i
           if i < d
              basis[1, t] = i + 1
           end
           i += 1
        else
            j = findfirst(x -> basis[x, t-1] != 0, 1:n)
            basis[:, t] = basis[:, t-1]
            if j == 1
               basis[1, t] -= 1
               basis[2, t] += 1
            else
               basis[1, t] = basis[j, t] - 1
               basis[j, t] = 0
               basis[j+1, t] += 1
            end
        end
    end
    return basis
end

function newton_cyclic(p, n, d)
    pbasis = get_fbasis(n, d)[:, 2:end]
    basis = [Int[]]
    csupp = cc.(p.supp, n)
    pushfirst!(csupp, zeros(Int, n))
    sort!(csupp)
    unique!(csupp)
    for a in eachcol(pbasis)
        if remove(csupp, a)
            append!(basis, permutation(a))
        end
    end
    return sort!(basis)
end

function newton_ncbasis(p)
    nbasis = [Int[]]
    for bi in p.supp
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

function get_ncbasis(n, d; ind=Vector(1:n), binary=false)
    basis = [Int[]]
    for i = 1:d
        append!(basis, _get_ncbasis_deg(n, i, ind=ind, binary=binary))
    end
    return basis
end

function _get_ncbasis_deg(n, d; ind=Vector(1:n), binary=false)
    if d > 0
        basis = Vector{Int}[]
        for i = 1:n
            temp = _get_ncbasis_deg(n, d-1, ind=ind, binary=binary)
            if binary == false || d == 1
                push!.(temp, ind[i])
                append!(basis, temp)
            else
                for item in temp
                    if item[end] != ind[i]
                        push!(basis, [item; ind[i]])
                    end
                end
            end
        end
        return basis
    else
        return [Int[]]
    end
end

function constraint_reduce!(word::Vector{Int}; constraint="unipotent")
    if constraint !== nothing
        i = 1
        while i < length(word)
            if word[i] == word[i+1]
                deleteat!(word, i)
                if constraint == "unipotent"
                    deleteat!(word, i)
                    if i > 1
                        i -= 1
                    end
                end
            else
                i += 1
            end
        end
    end
    return word
end

function reduce!(word::Vector{Int}; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    if obj == "trace"
        word = min(_cyclic_canon(word), _cyclic_canon(reverse(word)))
    else
        if constraint === nothing
            word = min(_comm(word, partition, comm_var), _comm(reverse(word), partition, comm_var))
        else
            word = min(constraint_reduce!(_comm(word, partition, comm_var), constraint = constraint), constraint_reduce!(_comm(reverse(word), partition, comm_var), constraint = constraint))
        end
    end
    return word
end

function wreduce(w::NCMono, x; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    ind = w.z .> 0
    vars = w.vars[ind]
    exp = w.z[ind]
    word = Int[]
    for j = 1:length(vars)
        k = bfind(x, vars[j], rev=true)
        append!(word, k*ones(Int, exp[j]))
    end
    word = reduce!(word, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    return prod(x[word])
end

function _comm(word::Vector{Int}, partition, comm_var)
    if partition > 0
        w1 = copy(word[word .<= partition])
        w2 = word[word .> partition]
    else
        w1 = copy(word)
        w2 = Int[]
    end
    if comm_var > 0
        i = 1
        while i < length(w1)
            if w1[i] <= comm_var && w1[i+1] <= comm_var && w1[i] > w1[i+1]
                w1[i],w1[i+1] = w1[i+1],w1[i]
                if i > 1
                    i -= 1
                else
                    i = 2
                end
            else
                i += 1
            end
        end
    end
    return [w1; w2]
end

function bfind(A, a; rev=false)
    low = 1
    high = length(A)
    while low <= high
        mid = Int(ceil((low + high)/2))
        if A[mid] == a
           return mid
        elseif A[mid] < a
            if rev == false
                low = mid + 1
            else
                high = mid - 1
            end
        else
            if rev == false
                high = mid - 1
            else
                low = mid + 1
            end
        end
    end
    return nothing
end

function permutation(a)
    b = sparse(a)
    ua = convert(Vector{Int}, b.nzind)
    na = convert(Vector{Int}, b.nzval)
    return _permutation(ua, na)
end

function _permutation(ua, na)
    if !isempty(ua)
        perm = Vector{Int}[]
        for i = 1:length(ua)
            nua = copy(ua)
            nna = copy(na)
            if na[i] == 1
                deleteat!(nua, i)
                deleteat!(nna, i)
            else
                nna[i] -= 1
            end
            temp = _permutation(nua, nna)
            push!.(temp, ua[i])
            append!(perm, temp)
        end
        return perm
    else
        return [Int[]]
    end
end

function sym_inv(word; vargroup=nothing)
    if vargroup !== nothing
        cword = copy(word)
        ind = [gind(a, vargroup) for a in word]
        uind = unique(ind)
        nind = [count(ind .== a) for a in uind]
        k = 0
        for i = 1:length(uind)
            cword[k+1:k+nind[i]] = reverse(cword[k+1:k+nind[i]])
            k += nind[i]
        end
        return min(word, cword)
    else
        return min(word, reverse(word))
    end
end

function sym_cyclic(word; vargroup=nothing, constraint=nothing)
    if constraint !== nothing
        while length(word) > 2 && word[1] == word[end]
            if constraint == "unipotent"
                word = word[2:end-1]
            elseif constraint == "projection"
                word = word[1:end-1]
            end
        end
    end
    return min(_cyclic_canon(word), _cyclic_canon(reverse(word)))
end

function iscomm(a; vargroup=nothing)
    if vargroup !== nothing
        for i = 1:length(a)-1
            if a[i] > a[i+1] && gind(a[i], vargroup) != gind(a[i+1], vargroup)
                return true
            end
        end
    end
    return false
end

function gind(k, vargroup)
    return findfirst(i -> k <= sum(vargroup[1:i]), 1:length(vargroup))
end

function res_comm!(a, vargroup)
    i = 1
    while i < length(a)
        if a[i] > a[i+1] && gind(a[i], vargroup) != gind(a[i+1], vargroup)
            a[i], a[i+1] = a[i+1], a[i]            
            if i > 1
                i -= 1
            end
        else
            i += 1
        end
    end
    return a
end

# generate an SOHS polynomial with variables vars and degree 2d
function add_SOHS!(model, vars, d; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    basis = vcat([MultivariatePolynomials.monomials(vars, i) for i = 0:d]...)
    if constraint !== nothing
        basis = basis[[all(item.z .< 2) for item in basis]]
    end
    if partition > 0
        ind = Int[]
        for (i, item) in enumerate(basis)
            vs = item.vars[item.z .> 0]
            if findfirst(j -> vs[j] < vars[partition] && vs[j+1] >= vars[partition], 1:length(vs)-1) === nothing
                push!(ind, i)
            end
        end
        basis = basis[ind]
    end
    sohs = 0
    pos = @variable(model, [1:length(basis), 1:length(basis)], PSD)
    for j = 1:length(basis), k = j:length(basis)
        word = wreduce(star(basis[j])*basis[k], vars, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
        if j == k
            @inbounds sohs += pos[j,k]*word
        else
            @inbounds sohs += 2*pos[j,k]*word
        end
    end
    return sohs
end

# generate a polynomial with variables vars and degree d
function add_poly!(model, vars, d; partition=0, constraint=nothing)
    basis = vcat([MultivariatePolynomials.monomials(vars, i) for i = 0:d]...)
    if constraint !== nothing
        basis = basis[[all(item.z .< 2) for item in basis]]
    end
    if partition > 0
        ind = Int[]
        for (i, item) in enumerate(basis)
            vs = item.vars[item.z .> 0]
            if findfirst(j -> vs[j] < vars[partition] && vs[j+1] >= vars[partition], 1:length(vs)-1) === nothing
                push!(ind, i)
            end
        end
        basis = basis[ind]
    end
    free = @variable(model, [1:length(basis)])
    poly = free'*basis
    return poly
end

function arrange(p, vars; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    mons = monomials(p)
    coe = coefficients(p)
    mons = [wreduce(mon, vars, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint) for mon in mons]
    nmons = unique(sort(mons))
    ncoe = zeros(typeof(coe[1]), length(nmons))
    for (i, item) in enumerate(coe)
        Locb = bfind(nmons, mons[i])
        ncoe[Locb] += item
    end
    return nmons,ncoe
end
