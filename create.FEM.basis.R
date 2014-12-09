create.FEM.basis <- function(p, e, t, order = 2, dl = NULL)
{
  if(is.null(t))
  {
    stop('Less than three input arguments.')
  }
  
  if(dim(t)[2] == 3)
  {
    t = cbind(t,matrix(1,dim(t)[1],1))
  }
  
  if(dim(t)[1] == 3)
  {
    t = rbind(t,matrix(1,dim(t)[1],1))
  }
  
  if(is.null(e))
  {
    if(dim(p)[1] != 2 & dim(p)[2] != 2 | dim(t)[1] != 4 & dim(t)[2] != 4)
    {
      stop('Dimensions of at least one of P,E and T are not correct')
    }
  }
  else
  {
    if(dim(p)[1] != 2 & dim(p)[2] != 2 | dim(e)[1] != 2 & dim(e)[2] != 2 | dim(t)[1] != 4 & dim(t)[2] != 4)
    {
      stop('Dimensions of at least one of P,E and T are not correct')
    }
  }
  
  if(dim(p)[2] !=2 & dim(p)[1] == 2)
  {
    p = t(p)
  }
  
  if(dim(e)[2] != 2 & dim(e)[1] == 2)
  {
    e = t(e)
  }
  
  if(dim(t)[2] != 4 & dim(t)[1] == 4)
  {
    t = t(t)
  }
  
  type = 'FEM'
  
  rangeval = NULL
  
  nodeStruct = makenodes(p, t[,1:3], order)
  
  petstr = list()
  
  petstr$p = p
  petstr$e = e
  petstr$t = t
  petstr$order = order
  petstr$nodes = nodeStruct$nodes
  petstr$nodeindex = nodeStruct$nodeindex
  petstr$J = nodeStruct$J
  petstr$metric = nodeStruct$metric
  petstr$dl = dl
  
  params = petstr
  
  nbasis = dim(nodeStruct$nodes)[1]
  
  basisobj = basisfd(type, rangeval, nbasis, params)
  
  return(basisobj)
  
}