using DynamicPolynomials
using NCTSSOS


# CHSH inequality
@ncpolyvar x[1:2]
@ncpolyvar y[1:2]
f = x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2]

opt,data = nctssos_first([-f], [x;y], 1, TS=false, partition=2, constraint="unipotent", Gram=true)

# I_3322 inequality
@ncpolyvar x[1:3]
@ncpolyvar y[1:3]
f = x[1]*(y[1] + y[2] + y[3]) + x[2]*(y[1] + y[2] - y[3]) + x[3]*(y[1] - y[2]) - x[1] - 2*y[1] - y[2]

opt,data = nctssos_first([-f], [x;y], 2, TS=false, partition=3, normality=1, constraint="projection", Gram=true)
