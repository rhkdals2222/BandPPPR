library(MCMCpack)
library(Rfast)
library(bdsmatrix)
library(purrr)
library(FinCovRegularization)


band_Khare=function(X,bandwidth,pn=1000,keepchol=FALSE,U=diag(epsilon,p),epsilon=10^(-4),L=diag(1,p),D=diag(1,p),alpha=rep(1,p)){
  p=dim(X)[2]
  n=dim(X)[1]
  k=bandwidth
  tild_U=t(X)%*%X+U
  tild_alpha=n+alpha
  
  Sigma=L%*%D%*%t(L)
  
  psamples=list()
  
  ##index matrix 만들기
  
  Lind=lower.tri(matrix(1,p,p)) & banding(matrix(1,p,p),k)
  Lindw1=!Lind & lower.tri(matrix(1,p,p))
  for(it in 1:pn){
    cat(it,"\n")
    Linv=forwardsolve(L, x = diag(ncol(L)))
    Lindw2=(t(Linv)==0)&upper.tri(matrix(1,p,p)) 
    
    temp1=mat.mult(Linv,tild_U) ##L^{-1}*\tilde{U}
    temp2=mat.mult(temp1,backsolve(r = t(L), x = diag(ncol(t(L))))) ##L^{-1}*\tilde{U}*(L^{T})^{-1}
    temp3=spdinv(Sigma) ##(LDL^{T})^{-1}
    
    
    ##M_vG 만들기
    Mv_G=map(1:(p-1),~temp3[Lind[,.x],Lind[,.x]]*temp2[.x,.x]) %>% map(spdinv)
    Lindw=Lindw2|Lindw1
    
    #w반영하기 
    wfun=function(x,y){
      uind=Lind[,x];
      wind=Lindw[,x];
      
      temp1[x,uind]/temp2[x,x]+y%*% temp3[uind,wind] %*% t(temp1)[wind,x]
    }
    ##w가 빈 index인 경우 처리해줘야함
    mu_uv=map2(1:(p-1),Mv_G,~wfun(.x,.y))
    L[Lind]=unlist(map2(mu_uv,Mv_G,~rmvnorm(1,.x,.y)))

    ##D matrix
    D=diag(map2_dbl(tild_alpha/2-1,diag(temp2)/2,~rinvgamma(1,shape=.x,scale=.y)))
    
    Sigma=mat.mult(L*outer(rep(1,dim(D)[1]),diag(D)), t(L)) ##L%*%D%*%t(L)
    
    #Sigma=L%*%D%*%t(L)
    if(keepchol){
      psamples[[it]]=list(L=L,D=D,Sigma=Sigma)
    }else{
      psamples[[it]]=Sigma
    }
  }
  return(psamples)
}




