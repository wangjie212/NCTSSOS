# TODO: rest of the code has not been refactored yet

function reduce!(word::Vector{UInt16}; obj="eigen", partition=0, constraint=nothing)
    if obj == "trace"
        word = min(_cyclic_canon(word), _cyclic_canon(reverse(word)))
    else
        if partition > 0 && constraint === nothing
            word = min(_comm(word, partition), _comm(reverse(word), partition))
        elseif partition == 0 && constraint !== nothing
            cword = copy(word)
            word = min(
                reduce_cons!(word; constraint=constraint),
                reduce_cons!(reverse(cword); constraint=constraint),
            )
        elseif partition > 0 && constraint !== nothing
            word = min(
                reduce_cons!(_comm(word, partition); constraint=constraint),
                reduce_cons!(_comm(reverse(word), partition); constraint=constraint),
            )
        else
            word = _sym_canon(word)
        end
    end
    return word
end

function reduce(word::Monomial{false}, x; obj="eigen", partition=0, constraint=nothing)
    if obj == "trace"
        word = min(_cyclic_canon(word), _cyclic_canon(star(word)))
    else
        if partition > 0 && constraint === nothing
            word = min(_comm(word, x, partition), _comm(star(word), x, partition))
        elseif partition == 0 && constraint !== nothing
            word = min(
                reduce_cons(word; constraint=constraint),
                reduce_cons(star(word); constraint=constraint),
            )
        elseif partition > 0 && constraint !== nothing
            word = min(
                reduce_cons(_comm(word, x, partition); constraint=constraint),
                reduce_cons(_comm(star(word), x, partition); constraint=constraint),
            )
        else
            word = _sym_canon(word)
        end
    end
    return word
end

function _comm(word::Vector{UInt16}, partition)
    ind1 = word .<= partition
    ind2 = word .> partition
    return [word[ind1]; word[ind2]]
end

function _comm(w::Monomial{false}, x, partition)
    ind1 = w.vars .>= x[partition]
    ind2 = w.vars .< x[partition]
    return prod([w.vars[ind1]; w.vars[ind2]] .^ [w.z[ind1]; w.z[ind2]])
end

function permutation(a)
    b = sparse(a)
    ua = convert(Vector{UInt16}, b.nzind)
    na = convert(Vector{UInt16}, b.nzval)
    return _permutation(ua, na)
end

function _permutation(ua, na)
    if !isempty(ua)
        perm = Vector{UInt16}[]
        for i in 1:length(ua)
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
        return [UInt16[]]
    end
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

function _cyclic_basis(var)
    basis = _permutation(var, ones(length(var)))
    return unique(_cyclic_canon.(basis))
end

function sym(word, vargroup)
    cword = copy(word)
    ind = [gind(word[i], vargroup) for i in 1:length(word)]
    uind = unique(ind)
    nind = [count(ind .== uind[i]) for i in 1:length(uind)]
    k = 0
    for i in 1:length(uind)
        cword[(k + 1):(k + nind[i])] = reverse(cword[(k + 1):(k + nind[i])])
        k += nind[i]
    end
    return min(word, cword)
end

function iscomm(a, vargroup)
    for i in 1:(length(a) - 1)
        if a[i] > a[i + 1] && gind(a[i], vargroup) != gind(a[i + 1], vargroup)
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
        if a[i] > a[i + 1] && gind(a[i], vargroup) != gind(a[i + 1], vargroup)
            temp = a[i]
            a[i] = a[i + 1]
            a[i + 1] = temp
            if i > 1
                i -= 1
            end
        else
            i += 1
        end
    end
    return a
end

function issym(word, vargroup)
    ind = [gind(word[i], vargroup) for i in 1:length(word)]
    uind = unique(ind)
    nind = [count(ind .== uind[i]) for i in 1:length(uind)]
    k = 0
    for i in 1:length(uind)
        temp = word[(k + 1):(k + nind[i])]
        if reverse(temp) != temp
            return false
        end
        k += nind[i]
    end
    return true
end

function star(w::Monomial{false})
    return prod(reverse(w.vars) .^ reverse(w.z))
end

function star(p::Polynomial{false})
    return coefficients(p)' * star.(monomials(p))
end

# generate an SOHS polynomial with variables vars and degree 2d
function add_SOHS!(model, vars, d; obj="eigen", partition=0, constraint=nothing)
    basis = vcat([MultivariatePolynomials.monomials(vars, i) for i in 0:d]...)
    if constraint !== nothing
        basis = basis[[all(item.z .< 2) for item in basis]]
    end
    if partition > 0
        ind = Int[]
        for (i, item) in enumerate(basis)
            vs = item.vars[item.z .> 0]
            if findfirst(
                j -> vs[j] < vars[partition] && vs[j + 1] >= vars[partition],
                1:(length(vs) - 1),
            ) === nothing
                push!(ind, i)
            end
        end
        basis = basis[ind]
    end
    sohs = 0
    pos = @variable(model, [1:length(basis), 1:length(basis)], PSD)
    for j in 1:length(basis), k in j:length(basis)
        word = reduce(
            star(basis[j]) * basis[k],
            vars;
            obj=obj,
            partition=partition,
            constraint=constraint,
        )
        if j == k
            @inbounds sohs += pos[j, k] * word
        else
            @inbounds sohs += 2 * pos[j, k] * word
        end
    end
    return sohs
end

# generate a polynomial with variables vars and degree d
function add_poly!(model, vars, d; obj="eigen", partition=0, constraint=nothing)
    basis = vcat([MultivariatePolynomials.monomials(vars, i) for i in 0:d]...)
    if constraint !== nothing
        basis = basis[[all(item.z .< 2) for item in basis]]
    end
    if partition > 0
        ind = Int[]
        for (i, item) in enumerate(basis)
            vs = item.vars[item.z .> 0]
            if findfirst(
                j -> vs[j] < vars[partition] && vs[j + 1] >= vars[partition],
                1:(length(vs) - 1),
            ) === nothing
                push!(ind, i)
            end
        end
        basis = basis[ind]
    end
    free = @variable(model, [1:length(basis)])
    poly = free' * basis
    return poly
end

function arrange(p, vars; obj="eigen", partition=0, constraint=nothing)
    mons = monomials(p)
    coe = coefficients(p)
    mons = [
        reduce(mon, vars; obj=obj, partition=partition, constraint=constraint) for
        mon in mons
    ]
    nmons = unique(sort(mons))
    ncoe = zeros(typeof(coe[1]), length(nmons))
    for (i, item) in enumerate(coe)
        Locb = bfind(nmons, length(nmons), mons[i])
        ncoe[Locb] += coe[i]
    end
    return nmons, ncoe
end
