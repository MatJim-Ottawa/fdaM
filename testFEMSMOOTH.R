# FEM Smooth Test file

# library 
library("fda")

# source files

# basisfd seems to be internal to the FDA package, including modified version here...
source("basisfd.R")

source("create.FEM.basis.R")
source("eval.FEM.basis.R")
source("insideIndex.R")
source("makenodes.R")
source("plot.FEM.R")
source("stiff1.R")
source("smooth.FEM.basis.R")
source("tricoefCal.R")
source("mass.R")
# Import DATA from Matlab

precdata = read.csv("data_prec.csv", header = FALSE)
p = read.csv("mesh_P.csv", header = FALSE)
e = read.csv("mesh_E.csv", header = FALSE)
t = read.csv("mesh_T.csv", header = FALSE)

t = cbind(t[,1],t[,2],t[,3])
precdata = cbind(precdata[,1],precdata[,2],precdata[,3],precdata[,4],precdata[,5],precdata[,6],precdata[,7],precdata[,8],precdata[,9],precdata[,10],precdata[,11],precdata[,12])
lat = read.csv("mesh_Lat.csv", header = FALSE)
lng = read.csv("mesh_Lng.csv", header = FALSE)
lat = matrix(lat[,1])
lng = matrix(lng[,1])
# Create basis

order = 1

basisobj = create.FEM.basis(p,e,t,order)

# Precipitation

lambda = 1e-2

res = smooth.FEM.basis(lng,lat,precdata,basisobj,lambda)

print(res$gcv)
print(res$df)
print(sqrt(res$SSE/35))

# Plot the results

AZ = 15;
EL = 25;

