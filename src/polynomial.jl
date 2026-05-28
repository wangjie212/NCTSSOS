const NCMono = DP.Monomial{DP.NonCommutative{DP.CreationOrder}, Graded{LexOrder}}
const NCTerm{T} = Term{T, NCMono}
const NCPoly{T} = DP.Polynomial{DP.NonCommutative{DP.CreationOrder}, Graded{LexOrder}, T}
const NCPolyLike = Union{NCMono,NCTerm,NCPoly}
export NCMono, NCPoly
export ncpoly, statepoly, ncstatepoly

abstract type poly end

mutable struct ncpoly{T} <: poly
    supp::Vector{Vector{Int}}
    coe::Vector{T}
end

mutable struct statepoly{T} <: poly
    supp::Vector{Vector{Vector{Int}}}
    coe::Vector{T}
end

mutable struct ncstatepoly{T} <: poly
    supp::Vector{Tuple{Vector{Int}, Vector{Vector{Int}}}}
    coe::Vector{T}
end

function maxdeg(p::ncpoly)
    return maximum(length.(p.supp))
end

function maxdeg(p::statepoly{T}) where {T<:Number}
    return maximum([sum(length.(item)) for item in p.supp])
end

function maxdeg(p::ncstatepoly{T}) where {T<:Number}
    return maximum([length(item[1]) + sum(length.(item[2])) for item in p.supp])
end

function ncpoly(p::NCPoly{T}, x) where {T<:Union{Number,AffExpr,GenericAffExpr}}
    coe = MP.coefficients(p)
    mons = MP.monomials(p)
    supp = [Int[] for i=1:length(mons)]
    for (i, mon) in enumerate(mons)
        ind = mon.z .> 0
        vars = mon.vars[ind]
        exp = mon.z[ind]
        for j in eachindex(vars)
            k = bfind(x, vars[j], rev=true)
            append!(supp[i], k*ones(Int, exp[j]))
        end
    end
    return ncpoly{T}(supp, coe)
end

function star(w::NCMono)
    return prod(reverse(w.vars).^reverse(w.z))
end

function star(p::NCPoly{T}) where T<:Union{Number, AffExpr}
    return coefficients(p)'*star.(monomials(p))
end

function star(a::Tuple{Vector{Int}, Vector{Vector{Int}}})
    return tuple(reverse(a[1]), a[2])
end

function star(a::Vector{Vector{Int}})
    return a
end
