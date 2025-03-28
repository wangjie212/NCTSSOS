using DynamicPolynomials
using NCTSSOS


# CHSH inequality
@ncpolyvar x[1:2]
@ncpolyvar y[1:2]
f = x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2]

opt,data = nctssos_first([-f], [x;y], 1, TS=false, partition=2, constraint="unipotent", Gram=true)

vars = [x;y]
b = 1*prod.([vars[item] for item in data.basis[1]])
ib = 1*prod.([vars[reverse(item)] for item in data.basis[1]])
mon,coe = arrange(f + opt + ib'*data.GramMat[1][1]*b, [x;y], partition=2, constraint="unipotent")

# I_3322 inequality
@ncpolyvar x[1:3]
@ncpolyvar y[1:3]
f = x[1]*(y[1] + y[2] + y[3]) + x[2]*(y[1] + y[2] - y[3]) + x[3]*(y[1] - y[2]) - x[1] - 2*y[1] - y[2]

opt,data = nctssos_first([-f], [x;y], 2, TS=false, normality=false, partition=3, constraint="projection", Gram=true)

vars = [x;y]
b = 1*prod.([vars[item] for item in data.basis[1]])
ib = 1*prod.([vars[reverse(item)] for item in data.basis[1]])
mon,coe = arrange(f + opt + ib'*data.GramMat[1][1]*b, [x;y], partition=3, constraint="projection")
