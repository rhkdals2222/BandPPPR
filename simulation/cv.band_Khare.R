
source("simulation/band_Khare.R")
cv.band_Khare=function(X,kvec,n.cv,epsilon=10^(-4),nu0=p,A0=diag(epsilon,p),normtype="F",initSig){
  n <- dim(X)[1]
  n1 <- ceiling(n*(1 - 1/log(n)))
  n2 <- floor(n/log(n))
  
  choldecomp=chol(initSig)
  D=diag((diag(choldecomp))^2,dim(initSig)[1])
  L=t(solve(sqrt(D))%*%choldecomp)
  
  tempfun=function(index){
    S2=mat.mult(t(X[-index,]),X[-index,])/n2
    
    
    draw_all =band_Khare(datawt$X,k,L=L,D=D)
    
    #return(c(  draw_all %>% 
    #    map(function(sp){kvec %>% map(~banding(sp,k=.x)) %>% map_dbl(function(x){norm(x-S2,type=normtype)})}) %>% 
    #    do.call("rbind",.) %>% colMeans,kvec %>% map(~banding(S1,k=.x)) %>% map_dbl(function(x){norm(x-S2,type=normtype)})))
    return(draw_all %>% 
             map(function(sp){kvec %>% map(~banding(sp,k=.x)) %>% map(adjust_pd) %>% map_dbl(function(x){norm(x-S2,type=normtype)})}) %>% 
             do.call("rbind",.) %>% colMeans)
    
  }
  
  tres=map(1:n.cv,function(x){set.seed(x);sample(1:n, size = n1, replace = FALSE)}) %>%
    map(tempfun) %>% do.call("rbind",.) %>% colMeans
  
  return(tres)
  
}