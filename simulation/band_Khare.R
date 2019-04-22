library(MCMCpack)


band_Khare=function(X,bandwidth,pn=1000,keepchol=FALSE,U=diag(epsilon,p),epsilon=10^(-4),L=diag(1,p),D=diag(1,p),alpha=rep(p,p)){
  p=dim(X)[2]
  n=dim(X)[1]
  k=bandwidth
  tild_U=t(X)%*%X+U
  tild_alpha=n+alpha
  
  psamples=list()
  for(it in 1:pn){
    Linv=forwardsolve(L, x = diag(ncol(L)))
    temp1=Linv%*%tild_U ##L^{-1}*\tilde{U}
    temp2=temp1%*%backsolve(r = t(L), x = diag(ncol(t(L)))) ##L^{-1}*\tilde{U}*(L^{T})^{-1}
    #temp3=solve(L%*%D%*%t(L)) ##(LDL^{T})^{-1}
    temp3=solve(gchol(L%*%D%*%t(L))) ##(LDL^{T})^{-1}, 이부분 속도 향상시키기 
    
    for(v in 1:(p-1)){
      mind=seq(from=(v+1),by=1,length=k)
      mind=mind[mind<=p]
      mindp=mind-min(mind)+1
      M_vG_inv=temp3[mind,mind]*temp2[v,v]
      M_vG=solve(M_vG_inv)
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
    
    ##D_{ii} sampling
    map_dbl(tild_alpha/2-1)
    
    
    D=diag(map2_dbl(tild_alpha/2-1,diag(temp2/2),~rinvgamma(1,shape=.x,scale=.y)))
    if(keepchol){
      psamples[[it]]=list(L=L,D=D,Sigma=L%*%D%*%t(L))
    }else{
      psamples[[it]]=L%*%D%*%t(L)
    }
    
    
    
    
  }
  return(psamples)
}



