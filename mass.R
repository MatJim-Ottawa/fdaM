mass <- function(nodeStruct)
{
  order     = nodeStruct$order
  nodes     = nodeStruct$nodes
  nodeindex = nodeStruct$nodeindex
  Jvec      = nodeStruct$J
  
  nele = dim(nodeindex)[1]
  nnod = dim(nodes)[1]
  
  if (order == 2)
  {
    K0M = rbind(cbind( 6, 0, -1, -4, -1,  0),cbind( 0, 32,  0, 16, -4, 16),
                cbind(-1,  0,  6,  0, -1, -4),
                cbind(-4, 16,  0, 32,  0, 16),
                cbind(-1, -4, -1,  0,  6,  0),
                cbind( 0, 16, -4, 16,  0, 32))/360
    
    K0 = matrix(0,nnod, nnod)
    
    for(el in 1:nele )
    {
      ind = nodeindex[el,]
      K0[ind,ind]  = K0[ind,ind] + K0M*Jvec[el]
    }
  }
  else if(order == 1)
  {
    K0M = rbind(cbind(2,1,1),
                cbind(1,2,1),
                cbind(1,1,2))/24
    
    K0 = matrix(0,nnod, nnod)
    
    for(el in 1:nele)
    {
      ind = nodeindex[el,]
      K0[ind,ind] = K0[ind,ind] + K0M*Jvec[el]    }
  }
  else
  {
    stop('ORDER not 1 or 2.')
  }
}