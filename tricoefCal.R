tricoefCal <- function(p,t)
{
  ntri = dim(t)[1]
  
  tricoef = matrix(0,ntri,4)
  tricoef[,1] = p[t[,1],1]-p[t[,3],1];
  tricoef[,2] = p[t[,2],1]-p[t[,3],1];
  tricoef[,3] = p[t[,1],2]-p[t[,3],2];
  tricoef[,4] = p[t[,2],2]-p[t[,3],2];
  detT    = tricoef[,1]*tricoef[,4] - tricoef[,2]*tricoef[,3];
  tricoef = tricoef/(detT%*%matrix(1,1,4));
                     
  return(tricoef)
  
}


