##comparison


library(rslurm)


MLEcomp=function(n,p,k=NULL,method="band",filename=NULL,wd=NULL,dir=NULL){
  ##input : n,p,k,method
  ##output : Sigma0, MLE1, MLE2
  source("MLEexp/simulsetting.R")
  datawt=simulsetting(n=n,p=p,k=k,method=method)
  source("MLEexp/MLE_bandconst.R")
  #res=MLE_bandconst(datawt$X,k=k)
  if(!is.null(filename)){
    saveRDS(list(Sigma0=datawt$Sigma0,const_MLE=MLE_bandconst(datawt$X,k=k)$Sigma,MLE=t(datawt$X)%*%datawt$X/n),
            file = paste(wd,dir,filename,sep = "/"))  
  }else{
    return(list(Sigma0=datawt$Sigma0,const_MLE=MLE_bandconst(datawt$X,k=k)$Sigma,MLE=t(datawt$X)%*%datawt$X/n))  
  }
}


n=10
p=10
k=5
method="band"
mydir="testMLE"
#rr=MLEcomp(n,p,k)
siter=25
for(i in 1:siter){
  sjob<- slurm_call(MLEcomp, list(n=n,p=p,k=k,method=method,
                                  filename=paste0("MLE",i,".rds"),wd=getwd(),dir=mydir))
}
