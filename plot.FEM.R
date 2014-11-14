plot.FEM <- function(FEMfd, AZ = 37.5, EL = 30.0 , clims, label)
{
  coefmat = FEMfd$coefs
  nsurf   = dim(coefmat)[2]
  
  basisobj  = FEMfd$basis
  params    = basisobj$params
  
  p = params
  
}