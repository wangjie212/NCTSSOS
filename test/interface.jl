using Test, NCTSSOS
using Clarabel
using MosekTools
using Yao
using LinearAlgebra


@testset "Majumdar Gosh Model" begin
    num_sites = 16 
    J1_interactions = unique!([tuple(sort([i, mod1(i + 1, num_sites)])...) for i in 1:num_sites])
    J2_interactions = unique!([tuple(sort([i, mod1(i + 2, num_sites)])...) for i in 1:num_sites])

    J1 = 2.0
    J2 = 1.0

    # make_ham(num_sites) = sum([sum([J1 / 4 * ft_op[1] * kron(num_sites, i => ft_op[2], j => ft_op[2]) for ft_op in zip([1.0, 1.0, 1.0, -1.0], [X, Y, Z, I2])]) for (i, j) in J1_interactions]) +
    #                       sum([sum([J2 / 4 * ft_op[1] * kron(num_sites, i => ft_op[2], j => ft_op[2]) for ft_op in zip([1.0, 1.0, 1.0, -1.0], [X, Y, Z, I2])]) for (i, j) in J2_interactions])

    # myham = make_ham(num_sites)
    true_ans = -num_sites / 4 * 6 # I guess?
    # eigvals(Matrix(myham))[1]

    ij2idx_dict = Dict(zip([(i,j) for i in 1:num_sites, j in 1:num_sites if j > i], 1:(num_sites*(num_sites-1)÷2)))
    @ncpolyvar hij[1:(num_sites*(num_sites-1)÷2)]

    objective = (sum([J1 * hij[ij2idx_dict[(i,j)]] for (i,j) in J1_interactions]) + sum([J2 * hij[ij2idx_dict[(i,j)]] for (i,j) in J2_interactions]))

    gs = 
    [
         unique([
            (
                hij[ij2idx_dict[tuple(sort([i, j])...)]] * hij[ij2idx_dict[tuple(sort([k, l])...)]] - 
                hij[ij2idx_dict[tuple(sort([k, l])...)]] * hij[ij2idx_dict[tuple(sort([i, j])...)]]             
            ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites, l in 1:num_sites if
            (length(unique([i, j, k, l])) == 4)

         ]);
         unique([
            (
                hij[ij2idx_dict[tuple(sort([i, j])...)]] * hij[ij2idx_dict[tuple(sort([j, k])...)]] +
                hij[ij2idx_dict[tuple(sort([j, k])...)]] * hij[ij2idx_dict[tuple(sort([i, j])...)]] - 0.5 *(hij[ij2idx_dict[tuple(sort([i, j])...)]] + hij[ij2idx_dict[tuple(sort([j, k])...)]] - hij[ij2idx_dict[tuple(sort([i, k])...)]])
            ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
            (i != j && j != k && i != k)
        ])
    ]

    

    pop = PolyOpt(-objective; constraints=gs, is_equality=[true for _ in gs], is_projective=true)

    solver_config = SolverConfig(optimizer=Mosek.Optimizer; mom_order=1)

    result = cs_nctssos(pop, solver_config)
    @test isapprox(result.objective, true_ans; atol=1e-4)
end

@testset "Problem Creation Interface" begin
    n = 2
    @ncpolyvar x[1:n]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    # change struct name
    pop = PolyOpt(f; constraints=[g, h1], is_equality=[false, true])

    solver_config = SolverConfig(optimizer=Clarabel.Optimizer; mom_order=2, cs_algo=MF(), ts_algo=MMD())

    result = cs_nctssos(pop, solver_config)
    @test isapprox(result.objective, -1.0; atol=1e-4)

    result_higher = cs_nctssos_higher(pop, result, solver_config)
    @test isapprox(result.objective, result_higher.objective; atol=1e-4)
end


@testset "README Example Unconstrained" begin
    @ncpolyvar x[1:3]
    f = 1.0 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2] + x[2]*x[1] + x[2]*x[3] + x[3]*x[2]

    pop =  PolyOpt(f)

    solver_config_dense = SolverConfig(optimizer=Clarabel.Optimizer)

    result_dense = cs_nctssos(pop, solver_config_dense)

    result_cs = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF()))

    @test isapprox(result_dense.objective, result_cs.objective, atol=1e-4)

    result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))

    @test isapprox(result_cs.objective, result_cs_ts.objective, atol=1e-4)
end

@testset "README Example Constrained" begin
    @ncpolyvar x[1:2]
    f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1]*x[2] + x[2]*x[1] - 2.0

    pop =  PolyOpt(f; constraints=[g, h1], is_equality=[false, true])

    result_dense = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer))

    result_cs = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF()))

    @test isapprox(result_dense.objective, result_cs.objective, atol=1e-4)

    result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))

    @test isapprox(result_cs.objective, result_cs_ts.objective, atol=1e-4)

    result_cs_ts_higher = cs_nctssos_higher(pop, result_cs_ts, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))

    @test isapprox(result_dense.objective, result_cs_ts_higher.objective, atol=1e-4)
end
