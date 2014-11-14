smooth.FEM.basis <- function(Xvec, Yvec, data, basisobj, lambda = 1e-12, wtvec = NULL , covariates = NULL)
{
  if (is.null(basisobj))
  {
    stop('Less than four arguments supplied.')
  }
  
  N = dim(data);
  
  if( length(Xvec) != N[1] || dim(Xvec)[2] != 1 )
  {
    stop('XVEC is not a column vector of n X-coordinate values.')
  }
  
  if( length(Yvec) != N[1] || dim(Yvec)[2] != 1 )
  {
    stop('YVEC is not a column vector of n X-coordinate values.')
  }
  
  
  # Check BASISOBJ
  
  if( !is.basis(basisobj) )
  {
    stop('BASISOBJ is not a functional basis object.')
  }
  
  if( !basisobj$type == 'FEM')
  {
    stop('The basis object for BASISOBJ is not of type FEM.')
  }
  
  
  # LAMBDA
  
  if( !is.double(lambda) && dim(lambda) == NULL ) 
  {
    stop('LAMBDA is not numeric.')
  }
  else if(!all(dim(lambda) == cbind(1,2)))
  {
    error('LAMBDA is not a scalar.')
  }
  
  # check WTVEC
  
  wtvec = wtcheck(N[1],wtvec)
  
  # check covariates
  
  q = 0
  
  if( !is.null(covariates) )
  {
    if( !is.numeric(covariates))
    {
      stop('smooth.FEM.basis:covariates Optional argument COVARIATES is not numeric.')
    }
    
    if(dim(covariates)[1] != n)
    {
      stop('smooth.FEM.basis:covariates Optional argument COVARIATES has incorrect number of rows.')
    }
    
    q = dim(covariates)[2]
  }
  
  # Contruct projection matrix on the space spanned by the columns of the design matrix
  # DESMAT and remove projection of data on DESMAT from DATA
  
  if( matwt )
  {
    wtmat = vtvec
    wtfac = chol(wtmat)
  }
  else
  {
    wtmat = diag(wtvec)
    wtfac = sqrt(wtmat)
  }
  
  if( ! is.null(covariates) )
  {
    QR = qr(wtfac %*% covariates,)
    beta = solve(qr.R(QR), t(qr.Q(QR)) %*% (wtfac %*% data))
    H = qr.Q(QR) %*% qr.Q(QR)
    data = data - H %*% data
  }
  else
  {
    beta = NULL
  }
    
  # Set up n by nbasis matrix of basis values at observation points
  basismat = eval.FEM.basis(Xvec, Yvec, basisobj)
  
  # Set up the linear equations for smoothing
  params = basismat$params
  
  numnodes = dim(params$nodes)[1]
  
  indnodes = 1:numnodes
  
  # Extract quantities required for setting up mass and stiffness matrices
  
  nodeStruct$order      = params$order
  nodeStruct$nodes      = params$nodes
  nodeStruct$nodeindex  = params$nodeindex
  nodeStruct$J          = params$J
  nodeStruct$metric     = params$metric
  
  # construct mass matrix K0
  
  K0mat = mass(nodeStruct)
  
  # construct stiffness matrix K1
  
  K1mat = stiff1(nodeStruct)
  
  # construct the block diagonal matrix L, having upper left block given by I-H and zero otherwise
  
  # construct vector b for system Ax = b
  
  Bmat = matrix(0,numnodes*2, N[2])
  Bmat[indnodes,] = t(basismat) %*% wtmat %*% data
  
  # construct matrix A for system Ax = b
  
  Lmat = t(basismat) %*% wtmat %*% basismat
  Amat = rbind(cbind(Lmat,-lambda*K1mat),cbind(K1mat,K0mat))
  
  # solve system
  
  coefmat = solve(Amat,Bmat)
  
  # Make output function data objects SMOOTH_FD and LAPLACE_FD
  
  coef1 = coefmat[indnodes,]
  coef2 = coefmat[indnodes+numnodes,]
  
  smooth_fd = fd(coef1,basisobj)
  
  laplace_fd = fd(coef2, basisobj)
  
  # Compute SSE, yhat, GCV ect..
  
  # compute map from y to c
  
  Mmat = matrix(0,numnodes*2,N[1])
  Mmat[indnodes,] = t(basismat) %*% wtmat
  
  temp = solve(Amat, Mmat)
  
  y2cMap = temp[indnodes,]
  
  # compute degrees of freedom of smooth
  
  df = trace(y2cmat %*% basismat)
  
  # compute error sum of squares
  
  datahat = basismat %*% coef1
  
  SSE = sum(wtfac %*% (data - datahat) ^ 2)
  
  # compute GCV index
  
  if (df < N[1])
  {
    gcv = (SSE/N[1])/( (n - df)/n)^2
  }
  else
  {
    gcv = NaN
  }

  return(list(smoothFd = smooth_fd, laplaceFd = laplace_fd))
  
  
}