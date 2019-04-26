##cross validation
library(Rfast)
library(purrr)
library(CholWishart)
library(plyr)
library(RSpectra)


cv.BPPP=function(X,kvec,n.cv,piter=1000,epsilon=10^(-4),nu0=p,A0=diag(epsilon,p),normtype="F",wd,dir,filename=NULL){
  n <- dim(X)[1]
  p=dim(X)[2]
  n1 <- ceiling(n*(1 - 1/log(n)))
  n2 <- floor(n/log(n))
  
  source(paste0(wd,"/BPPP.R"))
  
  
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
    saveRDS(res,file = paste(wd,dir,filename,sep = "/"))  
  }else{
    return(res)  
  }
  
}
