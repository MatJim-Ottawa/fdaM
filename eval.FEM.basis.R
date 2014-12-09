eval.FEM.basis <- function(Xvec, Yvec, basisobj, nderivs = matrix(0,1,2))
{
  if(dim(Xvec)[2] == 2)
  {
    if(is.null(Yvec))
    {
      stop('First argument is a coordinate matrix and the second argument is not supplied')
    }
    
    basisobj = Yvec
    Yvec = Xvec[,2]
    Xvec = Xvec[,1]
  }
  else
  {
    if(is.null(basisobj))
    {
      stop('First and second arguments are coordinate vecotors and the third argument is not supplied.')
    }
  }
    
  if(basisobj$type != 'FEM')
  {
    stop('The basis object for BASISOBJ is not of type FEM.')
  }
  nbasis = basisobj$nbasis
  
  if(!is.double(Xvec))
  {
    stop('Xvec is not a numerical array')
  }
  else
  {
    Xvec = Xvec
  }
  
  #Check Yvec
  
  if(!is.double(Yvec))
  {
    stop('Yvec is not a numerical array.')
  }
  else if(length(Yvec) != length(Xvec))
  {
    stop('Yvec is not the same length as Xvec.')
  }
  else
  {
    Yvec = Yvec
  }
  
  if(length(nderivs) != 2)
  {
    stop('NDERIVS not of length 2.')
  }
  
  if(sum(nderivs)>2)
  {
    stop('Maximum derivative order is greater than two.')
  }
  
  N = length(Xvec)
  
  # Augment Xvec and Yvec by ones for computing barycentric coordinates
  Pgpts = cbind(matrix(1,N,1),Xvec,Yvec)
  
  # Get nodes and index
  
  params = basisobj$params
  p = params$p
  t = params$t
  t = t[,1:3]
  
  order = params$order
  nodes = params$nodes
  nodeindex = params$nodeindex
  Jvec = params$J
  
  # 1st, 2nd, 3rd vertices of triangles
  
  if(order == 2)
  {
    v1 = nodes[nodeindex[,1],]
    v2 = nodes[nodeindex[,3],]
    v3 = nodes[nodeindex[,5],]
  }
  else if(order == 1)
  {
    v1 = nodes[nodeindex[,1],]
    v2 = nodes[nodeindex[,2],]
    v3 = nodes[nodeindex[,3],]
  }
  else
  {
    stop('ORDER is neither 1 or 2.')
  }
  
  # Denominator of change of coordinates change matrix
  
  modJac = Jvec
  ones3 = matrix(1,3,1)
  modJacMat = modJac %*% t(ones3)
  
  M1 = cbind(v2[,1]*v3[,2] - v3[,1]*v2[,2], v2[,2] - v3[,2], v3[,1] - v2[,1])/modJacMat/2
  M2 = cbind(v3[,1]*v1[,2] - v1[,1]*v3[,2], v3[,2] - v1[,2], v1[,1] - v3[,1])/modJacMat/2
  M3 = cbind(v1[,1]*v2[,2] - v2[,1]*v1[,2], v1[,2] - v2[,2], v2[,1] - v1[,1])/modJacMat/2
  
  tricoef = tricoefCal(p,t)
  
  ind = matrix(0,N,1)
  for(i in 1:N)
  {
    ind[i] = insideIndex(Xvec[i], Yvec[i], p, t, tricoef)
  }
  
  # Sparse?
  evalmat = matrix(0,N,nbasis)
  
  for(i in 1:N)
  {
    indi = ind[i]
    
    if(!is.nan(indi))
    {
      baryc1 = (M1[indi,]*Pgpts[i,]) %*% ones3
      baryc2 = (M2[indi,]*Pgpts[i,]) %*% ones3
      baryc3 = (M3[indi,]*Pgpts[i,]) %*% ones3
      
      if(order == 2)
      {
        if(sum(nderivs) == 0)
        {
          evalmat[i,nodeindex[indi,1]] = 2* baryc1^2 - baryc1
          evalmat[i,nodeindex[indi,2]] = 2* baryc2^2 - baryc2
          evalmat[i,nodeindex[indi,3]] = 2* baryc3^2 - baryc3
          evalmat[i,nodeindex[indi,4]] = 4* baryc1 * baryc2
          evalmat[i,nodeindex[indi,5]] = 4* baryc2 * baryc3
          evalmat[i,nodeindex[indi,6]] = 4* baryc3 * baryc1
        }
        else if(nderivs[1] == 1 && nderivs[2] == 0)
        {
          evalmat[i,nodeindex[indi,1]] = (4* baryc1 - 1) * M1[indi,2]
          evalmat[i,nodeindex[indi,2]] = (4* baryc2 - 1) * M2[indi,2]
          evalmat[i,nodeindex[indi,3]] = (4* baryc3 - 1) * M3[indi,2]
          evalmat[i,nodeindex[indi,4]] = (4* baryc2 ) * M1[indi,2] + 4*baryc1 * M2[indi,2]
          evalmat[i,nodeindex[indi,5]] = (4* baryc3 ) * M2[indi,2] + 4*baryc2 * M3[indi,2]
          evalmat[i,nodeindex[indi,6]] = (4* baryc1 ) * M3[indi,2] + 4*baryc3 * M1[indi,2]
        }
        else if(nderivs[1] == 0 && nderivs[2] == 1)
        {
          evalmat[i,nodeindex[indi,1]] = (4*baryc1 - 1)*M1[indi,3]
          evalmat[i,nodeindex[indi,2]] = (4*baryc2 - 1)*M2[indi,3]
          evalmat[i,nodeindex[indi,3]] = (4*baryc3 - 1)*M3[indi,3]
          evalmat[i,nodeindex[indi,4]] = 4*baryc2*M1[indi,3] + 4*baryc1*M2[indi,3]
          evalmat[i,nodeindex[indi,5]] = 4*baryc3*M2[indi,3] + 4*baryc2*M3[indi,3]
          evalmat[i,nodeindex[indi,6]] = 4*baryc1*M3[indi,3] + 4*baryc3*M1[indi,3]
        }
        else if(nderivs[1] == 1 && nderivs[2] == 1)
        {
          evalmat[i,nodeindex[indi,1]] = 4*M1[indi,2]%*%M1[indi,3];
          evalmat[i,nodeindex[indi,2]] = 4*M2[indi,2]%*%M2[indi,3];
          evalmat[i,nodeindex[indi,3]] = 4*M3[indi,2]%*%M3[indi,3];
          evalmat[i,nodeindex[indi,4]] = 4*M2[indi,2]%*%M1[indi,3] + 4*M2[indi,3]%*%M1[indi,2];
          evalmat[i,nodeindex[indi,5]] = 4*M3[indi,2]%*%M2[indi,3] + 4*M3[indi,3]%*%M2[indi,2];
          evalmat[i,nodeindex[indi,6]] = 4*M1[indi,2]%*%M3[indi,3] + 4*M1[indi,3]%*%M3[indi,2];
        }
        else if(nderivs[1] == 2 && nderivs[2] == 0)
        {
          evalmat[i,nodeindex[indi,1]] = 4*M1[indi,2]%*%M1[indi,2];
          evalmat[i,nodeindex[indi,2]] = 4*M2[indi,2]%*%M2[indi,2];
          evalmat[i,nodeindex[indi,3]] = 4*M3[indi,2]%*%M3[indi,2];
          evalmat[i,nodeindex[indi,4]] = 8*M2[indi,2]%*%M1[indi,2];
          evalmat[i,nodeindex[indi,5]] = 8*M3[indi,2]%*%M2[indi,2];
          evalmat[i,nodeindex[indi,6]] = 8*M1[indi,2]%*%M3[indi,2];
        }
        else if(nderivs[1] == 0 && nderivs[2] == 2)
        {
          evalmat[i,nodeindex[indi,1]] = 4*M1[indi,3]%*%M1[indi,3];
          evalmat[i,nodeindex[indi,2]] = 4*M2[indi,3]%*%M2[indi,3];
          evalmat[i,nodeindex[indi,3]] = 4*M3[indi,3]%*%M3[indi,3];
          evalmat[i,nodeindex[indi,4]] = 8*M2[indi,3]%*%M1[indi,3];
          evalmat[i,nodeindex[indi,5]] = 8*M3[indi,3]%*%M2[indi,3];
          evalmat[i,nodeindex[indi,6]] = 8*M1[indi,3]%*%M3[indi,3];
        }
        else
        {
          stop('Inadmissible derivative orders.')
        }
      }
      else
      {
        if(sum(nderivs) == 0)
        {
          evalmat[i,nodeindex[indi,1]] = baryc1;
          evalmat[i,nodeindex[indi,2]] = baryc2;
          evalmat[i,nodeindex[indi,3]] = baryc3;
        }
        else if(nderivs[1] == 1 && nderivs[2] == 0)
        {
          evalmat[i,nodeindex[indi,1]] = M1[indi[1],2];
          evalmat[i,nodeindex[indi,2]] = M1[indi[2],2];
          evalmat[i,nodeindex[indi,3]] = M1[indi[3],2]; 
        }
        else if(nderivs[1] == 0 && nderivs[2] == 1)
        {
          evalmat[i,nodeindex[indi,1]] = M1[indi[1],3];
          evalmat[i,nodeindex[indi,2]] = M1[indi[2],3];
          evalmat[i,nodeindex[indi,3]] = M1[indi[3],3]; 
        }
        else
        {
          stop('Inadmissible derivative orders.')
        }
      }
    }
  }
  return(evalmat)
}