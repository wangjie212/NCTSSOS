abstract type AbstratMethod end

struct MomentSOS <: AbstratMethod 
	order::Int
	clique_func #function to get clique
end

get_order(method::MomentSOS) = method.order

get_clique_func(method::MomentSOS) = method.clique_func

function solve(method::MomentSOS, pop::PolynomialOptimizationProblem)
	order = method.order
	clique_func = method.clique_func
	
	

end


