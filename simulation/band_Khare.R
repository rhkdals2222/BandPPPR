library(MCMCpack)
library(Rfast)
library(bdsmatrix)
library(purrr)
library(FinCovRegularization)


band_Khare=post_khare=function(X,bandwidth,pn=1000,keepchol=FALSE,U=diag(epsilon,p),epsilon=10^(-4),L=diag(1,p),D=diag(1,p),alpha=rep(1,p)){
  p=dim(X)[2]
  n=dim(X)[1]
  k=bandwidth
  tild_U=t(X)%*%X+U
  tild_alpha=n+alpha
  
  Sigma=L%*%D%*%t(L)
  
  psamples=list()
  
  ##index matrix 만들기
  
  Lind=as.data.frame(lower.tri(matrix(1,p,p)) & banding(matrix(1,p,p),2))
  
  for(it in 1:pn){
    Linv=forwardsolve(L, x = diag(ncol(L)))
    temp1=mat.mult(Linv,tild_U) ##L^{-1}*\tilde{U}
    temp2=mat.mult(temp1,backsolve(r = t(L), x = diag(ncol(t(L))))) ##L^{-1}*\tilde{U}*(L^{T})^{-1}
    temp3=spdinv(Sigma) ##(LDL^{T})^{-1}
    
    ##M_vG_inv 만들기
    
    Mv_G=map(1:(p-1),temp3[Lind[,.x],Lind[,.x]]*temp2[.x,.x]) %>% map(spdinv)
    mu_uv=map(1:(p-1),temp1[.x,]/temp2[.x,.x])
    #w반영하기 
    
    map2(mu_uv,Mv_G,rmvnorm(1,.x,.y))
    
    #reduce to L
    
    
    
    for(v in 1:(p-1)){
      mind=seq(from=(v+1),by=1,length=k)
      mind=mind[mind<=p]
      mindp=mind-min(mind)+1
      M_vG_inv=temp3[mind,mind]*temp2[v,v]
      M_vG=spdinv(M_vG_inv)
      w=1:p
      wind= ((w>v) & (w-v)>k)|((w<v) & abs(Linv[v,w])<epsilon)
      tmu_uv=temp1[v,]/temp2[v,v]
      tmu_uv[!abs(Linv[v,w])<epsilon]=0
      
      mu_uv=tmu_uv[mind]
      
      for(pi in mindp){
        u=v+pi
        #tmu_uv
        #temp=mu_uv[pi]
        for(pip in mindp){
          up=v+pip
          if(sum(wind)>0){
            #cat(pi,"-",pip)
            mu_uv[pi]=mu_uv[pi]+M_vG[pi,pip]*temp2[v,v]*sum(temp3[up,wind]*tmu_uv[wind])
          }
        }
      }
      L[mind,v]=rmvnorm(1,mu_uv,M_vG)
    }
    
    
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




