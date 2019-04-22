library(mvtnorm)
library(FinCovRegularization)

##시뮬레이션을 위해 실제값과 그로부터 추출된 샘플값을 만든다.

simulsetting=function(n,p,k=NULL,alpha=0.1,rho=0.6,method="band",epsilon=10^(-4)){
  Sigma0=paramgen(p,k,alpha,rho,method,epsilon)
  X=rmvnorm(n,rep(0,p),Sigma0)
  return(list(X=X,Sigma0=Sigma0,k=k))
}

paramgen=function(p,k=NULL,alpha=0.1,rho=0.6,method="band",epsilon=10^(-4)){
  ##
  Sigma0=diag(1,p)
  Sigma0=(rho*(abs(row(Sigma0)-col(Sigma0)))^(-alpha-1))
  diag(Sigma0)=1
  if(method=="band"){
    Sigma0=as.matrix(banding(Sigma0,k))
  }
  if(min(eigen(as.matrix(Sigma0))$value)<epsilon){
    Sigma0=Sigma0 + diag(-min(eigen(as.matrix(Sigma0))$value)+epsilon,dim(Sigma0)[1])
  }
  return(Sigma0)
  
}

system("pwd",intern = T)
sys.frame(1)
