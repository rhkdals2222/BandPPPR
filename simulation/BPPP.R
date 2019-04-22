library(CholWishart)
bPPP=function(X,bandwidth,pn=1000,epsilon=10^(-4),df=p+1,A=diag(epsilon,p),adj=FALSE){
  
  p=dim(X)[2]
  n=dim(X)[1]
  k=bandwidth
  psample=alply(rInvWishart(pn, n+df, t(X)%*%X + A),3)
  pppsample=lapply(psample,function(x){banding(x,k);})
  if(adj){
    pppsample=lapply(pppsample,function(x){adjust_pd(x,outlist=FALSE);})  
  }
  
  
  return(list(psample,pppsample))
}

