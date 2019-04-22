library(CholWishart)
library(plyr)
bPPP=function(X,bandwidth,pn=1000,epsilon=10^(-4),df=p+1,A=diag(epsilon,p),adj=FALSE){
  
  p=dim(X)[2]
  n=dim(X)[1]
  k=bandwidth
  psample=alply(rInvWishart(pn, n+df, t(X)%*%X + A),3)
  pppsample=lapply(psample,function(x){banding(x,k);})
  if(adj){
    pppsample=lapply(pppsample,function(x){adjust_pd(x,outlist=FALSE);})  
  }
  
  
  return(list(IW=psample,BPPP=pppsample))
}

adjust_pd=function(Sigma,epsilon=10^(-4),outlist=TRUE){
  p=dim(Sigma)[1]
  ##성능향상위해 필요한 부분
  emin=min(eigen(Sigma)$values)
  Sigmaa=Sigma+diag((emin<0)*(-emin+epsilon),p)
  if(!outlist){
    return(Sigmaa)
  }
  tt=chol(Sigmaa)
  D=diag((diag(tt))^2,p)
  L=t(solve(sqrt(D))%*%tt)
  return(list(Sigma=Sigmaa,L=L,D=D))
}
