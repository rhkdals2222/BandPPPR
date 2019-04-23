

##main 함수
n=50
p=50
k=5

##참값과 자료를 만드는 부분
source("simulation/simulsetting.R")
datawt=simulsetting(n=n,p=p,k=k)

source("simulation/cv.band_Khare.R")
kharecv=cv.band_Khare(datawt$X,1:10,n.cv=2,normtype = "2",initSig = datawt$Sigma0)

source("simulation/cv.BPPP.R")

