##cross validation
library(Rfast)
library(purrr)
library(CholWishart)
library(plyr)
library(RSpectra)

source("simulation/BPPP.R") ## adjust를 만들기 위한것

cv.BPPP=function(X,kvec,n.cv,epsilon=10^(-4),nu0=p,A0=diag(epsilon,p),normtype="F",filename=NULL){
  n <- dim(X)[1]
  n1 <- ceiling(n*(1 - 1/log(n)))
  n2 <- floor(n/log(n))
  
  tempfun=function(index){
    S2=mat.mult(t(X[-index,]),X[-index,])/n2
    draw_all = alply(rInvWishart(piter,nu0+n1,A0+n1*mat.mult(t(X[index,]),X[index,])/n1),3)
    return(draw_all %>% 
               map(function(sp){kvec %>% map(~banding(sp,k=.x)) %>% map(adjust_pd) %>% map_dbl(function(x){norm(x-S2,type=normtype)})}) %>% 
               do.call("rbind",.) %>% colMeans)
  }
  res=map(1:n.cv,function(x){set.seed(x);sample(1:n, size = n1, replace = FALSE)}) %>%
           map(tempfun) %>% do.call("rbind",.) %>% colMeans
  
  if(!is.null(filename)){
    saveRDS(res,file = filename)  
  }else{
    return(res)  
  }
  
}
