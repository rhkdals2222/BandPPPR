##comparison

source("MLEexp/simulsetting.R")
n=30
p=30
k=5
datawt=simulsetting(n=n,p=p,k=k)


source("MLEexp/MLE_bandconst.R")
res=MLE_bandconst(datawt$X,k=k)
