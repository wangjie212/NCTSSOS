import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using DynamicPolynomials
using NCTSSOS

# CHSH, dense
println("********** CHSH, dense **********")

@ncpolyvar X[1:2]
@ncpolyvar Y[1:2]

CHSH = -symmetrize(X[1]*Y[1] +  X[1]*Y[2]  + X[2]*Y[1]  - X[2]*Y[2])

rational_certificate(CHSH, [], [], [X;Y], 2; partition=2, constraint="unipotent", QUIET=false, QUIETTS=true, tol=10e-30)

# CHSH, sparse

println("********** CHSH, sparse **********")

newbound, oldbound, shift = rational_certificate_sparse(
  CHSH,
  DynamicPolynomials.Polynomial[],    # now Vector{<:Polynomial} rather than Vector{Any}
  DynamicPolynomials.Polynomial[], 
  [X;Y],
  1;
  partition  = 2,
  constraint = "unipotent",
  QUIET      = false,
  tol        = 1e-20
)
newbound = -newbound
oldbound = -oldbound

# I3322, dense

println("********** I3322, dense **********")

@ncpolyvar x[1:6]
I3322 = -symmetrize(x[1]*(x[4] + x[5] + x[6]) + x[2]*(x[4] + x[5] - x[6]) + x[3]*(x[4] - x[5]) - x[1] - 2*x[4] - x[5])

newbound, oldbound, shift = rational_certificate(I3322, [], [], x, 3; partition=3, constraint="projector", QUIET=false, QUIETTS=true, tol=10e-20)

# I3322, sparse

println("********** I3322, sparse **********")

newbound_sparse, oldbound_sparse, diff_sparse = rational_certificate_sparse(
  I3322,
  DynamicPolynomials.Polynomial[],    # now Vector{<:Polynomial} rather than Vector{Any}
  DynamicPolynomials.Polynomial[], 
  x,
  2;
  partition  = 3,
  constraint = "projector",
  QUIET      = false,
  tol        = 1e-20
)

# I4322

println("********** I4322, dense **********")

@ncpolyvar x[1:7]

I4322 = -symmetrize(-2*x[1] - x[2] - x[3] - x[5] + x[1]*(x[5] + x[6] + x[7]) + x[2]*(x[5] - x[6] + x[7]) + x[3]*(x[5] + x[6] - x[7]) + x[4]*(x[5] - x[6] - x[7]))

rational_certificate(I4322, [], [], x, 2; partition=4, constraint="projector", QUIET=false, QUIETTS=true, tol=10e-20)
