makenodes <- function(p, t, order = 2)
{
  if(dim(p)[2] > 2)
  {
    nodes = t(p)
    t = t(t[1:3,])
  }
  else
  {
    nodes = p
  }
  
  nele = dim(t)[1]
  nver = dim(p)[1]
  
  Jvec = matrix(0,nele,1)
  metric = array(rep(0,nele*2*2),dim  = c(nele,2,2))
  
  if(order == 2)
  {
    rec = matrix(0,nver,nver)
    ind = rbind(c(1,2), c(2,3), c(3,1))
    nodeindex = matrix(0, nele, 6)
    nodeindex[,c(1,3,5)] = t
    
    
    for(i in 1:nele)
    {
      for(j in 1:3)
      {
        if(rec[t, ind[j,1], t[i, ind[j,2]]] == 0)
        {
          nodes = rbind(nodes,cbind(0.5, 0.5) %*% nodes[t[i,ind[j,]],])
          rec[t[i,ind[j,1]], t[i,ind[j,2]]] = dim(nodes)[1]
          rec[t[i,ind[j,2]], t[i,ind[j,1]]] = dim(nodes)[1]
          
          nodeindex[i,2*j] = dim(nodes)[1]
        }
        else
        {
          nodeindex[i,2*j] = rec[t[i,ind[j,1]], t[i,ind[j,2]]]
        }
      }
      
      diff1x = nodes[nodeindex[i,3],1] - nodes[nodeindex[i,1],1]
      diff1y = nodes[nodeindex[i,3],2] - nodes[nodeindex[i,1],2]
      diff2x = nodes[nodeindex[i,5],1] - nodes[nodeindex[i,1],1]
      diff2y = nodes[nodeindex[i,5],2] - nodes[nodeindex[i,1],2]
      
      
      Jvec[i] = (diff1x * diff2y - diff2x * diff1y) / 2
      
      Ae1 = rbind(c(diff2y, -diff1y), c(-diff2x, diff1x) )/Jvec[i]
      
      metric[i,,] = t(Ae1) %*% Ae1
      
    }
    
  }
  else if(order == 1)
  {
    nodeindex = t[,1:3]
    
    for(i in 1:nele)
    {
      diff1x = nodes[nodeindex[i,2],1] - nodes[nodeindex[i,1],1]
      diff1y = nodes[nodeindex[i,2],2] - nodes[nodeindex[i,1],2]
      diff2x = nodes[nodeindex[i,3],2] - nodes[nodeindex[i,1],2]
      diff2y = nodes[nodeindex[i,3],1] - nodes[nodeindex[i,1],1]
      
      Jvec[i] = (diff1x * diff2y - diff2x * diff1y) / 2
      
      Ae1 = rbind(c(diff2y, -diff1y), c(-diff2x, diff1x) )/Jvec[i]
      
      metric[i,,] = t(Ae1) %*% Ae1
    }
  }
  else
  {
    stop('ORDER not 1 or 2.')
  }
  
  nodeStruct$order = order
  nodeStruct$nodes = nodes
  nodeStruct$nodeindex = nodeindex
  nodeStruct$J = Jvec
  nodeStruct$metric = metric
  
  return(nodeStruct)
}