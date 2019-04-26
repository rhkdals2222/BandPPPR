##cv.bandEST


##cross validation
library(Rfast)
library(purrr)
library(plyr)
library(RSpectra)

source("simulation/BPPP.R") ## adjust를 만들기 위한것

cv.BPPP=function(X,kvec,n.cv=100,epsilon=10^(-4),nu0=p,A0=diag(epsilon,p),normtype="F",filename=NULL){
  n <- dim(X)[1]
  n1 <- ceiling(n*(1 - 1/log(n)))
  n2 <- floor(n/log(n))
  
  tempfun=function(index){
    S1=mat.mult(t(X[index,]),X[index,])/n1
    S2=mat.mult(t(X[-index,]),X[-index,])/n2
    
    return(kvec %>% map(~banding(S1,k=.x))%>% map(adjust_pd) %>% map_dbl(function(x){norm(x-S2,type=normtype)}))
  }
  res=map(1:n.cv,function(x){set.seed(x);sample(1:n, size = n1, replace = FALSE)}) %>%
    map(tempfun) %>% do.call("rbind",.) %>% colMeans
  if(!is.null(filename)){
    saveRDS(res,file = filename)  
  }else{
    return(res)  
  }
  
}
