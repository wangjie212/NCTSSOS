using NCTSSOS, DynamicPolynomials, Clarabel, MosekTools, CliqueTrees
N = 6
τ_tuple2idx = Dict(vec([(i, j) => (((i - 1) * (i - 2)) ÷ 2) + j for i in 1:N for j in 1:N if j < i]))
@polyvar τ[1:N*(N+1)÷2-N]
@polyvar x[1:N]
@polyvar y[1:N]
@polyvar z[1:N]

objective = sum(1.0 * τ .^ 6 .- 1.0 * τ .^ 3)

cons = [(τ[((i-1)*(i-2))÷2+j] * ((x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2) - 1) for i in 1:N for j in 1:N if j < i]
cons = [cons; -1 .* cons]

pop = PolynomialOptimizationProblem(objective, cons)

solver_config = SolverConfig(optimizer=Mosek.Optimizer, mom_order=3, cs_algo=MF(), ts_algo=MMD())

result = cs_nctssos(pop, solver_config)

result.objective

result_higher = cs_nctssos_higher(pop, result, solver_config)

result_higher.objective
