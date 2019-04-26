#simulation
library(rslurm)

##datagen
source(paste0(getwd(),"/simulsetting.R"))

n=100
p=100
k=5
siter=25
X_list=list()
for(i in 1:siter){
  data=simulsetting(n=n,p=p,k=k)
  X_list[[i]]=data$X
  source(paste0(getwd(),"/cv.BPPP.R"))
  sjob<- slurm_call(cv.BPPP, list(X=data$X,kvec=1:15,n.cv=100,filename=paste0("cvres",i,".rds"),normtype="2",wd=getwd(),dir="test"))
}

