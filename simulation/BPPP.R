##밴드 사후처리 

library(CholWishart)
library(plyr)
library(RSpectra)
library(purrr)


adjust_pd=function(Sigma,epsilon=10^(-4),outlist=FALSE){
  p=dim(Sigma)[1]
  ##성능향상위해 필요한 부분
  emin=eigs(Sigma,1,"SR")$values
  Sigmaa=Sigma+diag((emin<0)*(-emin+epsilon),p)
  if(!outlist){
    return(Sigmaa)
  }
  tt=chol(Sigmaa)
  D=diag((diag(tt))^2,p)
  L=t(solve(sqrt(D))%*%tt)
  return(list(Sigma=Sigmaa,L=L,D=D))
}

bPPP=function(X,bandwidth,pn=1000,epsilon=10^(-4),df=p,A=diag(epsilon,p),adj=FALSE){
  p=dim(X)[2]
  n=dim(X)[1]
  k=bandwidth
  psample=alply(rInvWishart(pn, n+df, t(X)%*%X + A),3)
  pppsample=lapply(psample,function(x){banding(x,k);})
  if(adj){
    pppsample=map(pppsample,adjust_pd)  
  }
  return(list(IW=psample,BPPP=pppsample))
}
