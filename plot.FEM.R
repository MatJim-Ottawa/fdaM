plot.FEM <- function(FEMfd, AZ = 37.5, EL = 30.0 , clims, label)
{
  coefmat = FEMfd$coefs
  nsurf   = dim(coefmat)[2]
  
  basisobj  = FEMfd$basis
  params    = basisobj$params
  
  p = params$p
  t = params$t[,1:3]
  
  
  for(isurf in 1:nsurf)
  {
    #geometry package: surf.tri
    print('implement plotting in plot.FEM')
  }
  
}