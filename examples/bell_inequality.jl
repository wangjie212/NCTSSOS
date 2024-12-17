using DynamicPolynomials
using NCTSSOS


# I_3322 inequality
@ncpolyvar x[1:3]
@ncpolyvar y[1:3]
f = x[1]*(y[1] + y[2] + y[3]) + x[2]*(y[1] + y[2] - y[3]) + x[3]*(y[1] - y[2]) - x[1] - 2*y[1] - y[2]

opt,data = nctssos_first([-f], [x;y], 3, TS=false, normality=true, partition=3, constraint="projection")
