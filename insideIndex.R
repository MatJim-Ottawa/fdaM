insideIndex <- function(X, Y, p, t, tricoef = NULL)
{
  small = 1000*.Machine$double.eps
  
  ntri = dim(t,1)[1]
  indtri = t(1:ntri)
  
  if (tricoef == NULL)
  {
    tricoef = matrix(0,ntri,4)
    triceof[,1] = p[t[,1],1] - p[t[,3],1]
    tricoef[,2] = p[t[,2],1] - p[t[,3],1]
    tricoef[,3] = p[t[,1],2] - p[t[,3],2]
    tricoef[,4] = p[t[,2],2] - p[t[,3],2]
    
    detT = tricoef[,1] * tricoef[,4] - tricoef[,2] * tricoef(,3)
    tricoef = tricoef / (detT %*% matrix(1,1,4))
    
  }
  
  r3 = X - p[t[,3],1]
  s3 = Y - p[t[,3],2]
  lam1 = ( triceof[,4]*r3 - tricoef[,2]*s3 )
  lam2 = ( - tricoef[,3] * r3 + tricoef[,1]*s3)
  
  lam3 = 1 - lam1 - lam2
  
  int = (-small <= lam1 & lam1 <= 1 + small ) & (-small <= lam2 & lam2 <= 1+small) & (-small <= lam3 & lam3 <= 1+small)
  
  
  indi = indtri(int)
  
  if(is.null(indi))
  {
    ind = NULL
  }
  else
  {
    ind = min(indi)
  }
  
}

