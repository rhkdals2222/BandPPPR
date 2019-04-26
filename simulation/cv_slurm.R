#simulation
library(rslurm)

##datagen
source("simulation/simulsetting.R")

n=20
p=20
k=5
#X_list=
data=simulsetting(n=n,p=p,k=k)

source("simulation/cv.bandEST.R")

sjob<- slurm_call(cv.BPPP, list(X=data$X,kvec=1:15,filename=paste0("cvres",1,".rds")))

