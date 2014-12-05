stiff1 <- function(nodeStruct)
{
  order = nodeStruct$order
  nodes = nodeStruct$nodes
  nodeindex = nodeStruct$nodeindex
  Jvec = nodeStruct$J
  metric = nodeStruct$metric
  
  nele = dim(nodeindex)[1]
  nnod = dim(nodes)[1]
  
  K1 = matrix(0,nnod,nnod)
  
  
  if(order == 2)
  {
    KXX = rbind(cbind( 3, -4, 1, 0, 0, 0),
                cbind(-4, 8, -4, 0, 0, 0),
                cbind( 1, -4, 3, 0, 0, 0),
                cbind( 0, 0, 0, 8, 0 -8),
                cbind( 0, 0, 0, 0, 0, 0),
                cbind( 0, 0, 0, -8, 0, 8))/6
    
    KXY = rbind(cbind( 3, 0, 0, 0, 1, -4),
                cbind(-4, 4, 0, -4, 0, 4),
                cbind( 1, -4, 0, 4, -1, 0),
                cbind( 0, -4, 0, 4, 4, -4),
                cbind( 0, 0, 0, 0, 0, 0),
                cbind( 0, 4, 0, -4, -4, 4))/6
    
    KYY = rbind(cbind( 3, 0, 0, 0, 1, -4),
                cbind( 0, 8, 0, -8, 0, 0),
                cbind( 0, 0, 0, 0, 0, 0),
                cbind( 0, -8, 0, 8, 0, 0),
                cbind( 1, 0, 0, 0, 3, -4),
                cbind(-4, 0, 0, 0, -4, 8))/6
    
    # Assemble the stiffness matrix
    
    for (el in 1:nele)
    {
      ind = nodeindex[el,]
      K1M = ( metric(el,1,1) * KXX + metric(el,1,2) * KXY + metric(el,2,1) * t(KXY) +  metric(el,2,2) * KYY )
      
      K1[ind,ind] = K1[ind,ind] + K1M * Jvec(el)
    }

  }
  else if (order == 1)
  {
    KXX = rbind(cbind( 1, -1, 0),
                cbind(-1, 1, 0),
                cbind( 0, 0, 0))/2
    
    KXY = rbind(cbind( 1, 0, -1),
                cbind(-1, 0, 1),
                cbind( 0, 0, 0))/2
    
    KYY = rbind(cbind( 1, 0, -1),
                cbind( 0, 0, 0),
                cbind(-1, 0, 1))/2
    
    for (el in 1: nele)
    {
      ind = nodeindex[el,]
      K1M = (metric(el,1,1)*KXX + metric(el,1,2)*KXY + metric(el,2,1)*t(KXY) + metric(el,2,2)*KYY)
      K1[ind,ind] = K[ind,ind] + K1M*Jvec(el)
    }
    
  }
  else
  {
    stop('ORDER not 1 or 2.')
  }
}