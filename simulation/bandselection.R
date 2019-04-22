##band selection

##map reduce로 구현하기
bandselection=function(X,kvec,n.cv,epsilon=10^(-4)){
  n <- dim(X)[1]
  n1 <- ceiling(n*(1 - 1/log(n)))
  n2 <- floor(n/log(n))
  nu0=p+1
  A0=diag(epsilon,p)
  
  kvec=1:10
  n.cv=10
  
  tempfun=function(index){
    S1=t(X[index,])%*%X[index,]/n1
    S2=t(X[-index,])%*%X[-index,]/n2
    
    draw_all = alply(rInvWishart(piter,nu0+n1,A0+n1*S1),3)
    
    tt=draw_all %>% 
      map(function(sp){kvec %>% map(~banding(sp,k=.x)) %>% map_dbl(function(x){norm(x-S2,type="2")})})
    return(c(colMeans(do.call("rbind",tt)),kvec %>% map(~banding(S1,k=.x)) %>% map_dbl(function(x){norm(x-S2,type="2")})))
    
  }
  
  res=map(1:n.cv,function(x){set.seed(x);sample(1:n, size = n1, replace = FALSE)}) %>%
    map(tempfun)
  
  
  
}
