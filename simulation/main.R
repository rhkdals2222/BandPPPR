
##main 함수
n=100
p=200
k=10

##참값과 자료를 만드는 부분
source("simulation/simulsetting.R")
datawt=simulsetting(n=n,p=p,k=k)

##밴드크기 알고있다는 가정하에 추론하기
source("simulation/BPPP.R")
bPPPsample=bPPP(datawt$X,bandwidth=k,adj=TRUE)
S=t(datawt$X)%*%datawt$X/n
bandingEST=adjust_pd(banding(S,k),outlist=FALSE)
bandingEST=banding(S,k)
##추론된 값에 대해 성능 계산하기 

BPPPmean=Reduce("+",bPPPsample$BPPP)/length(bPPPsample$BPPP)  
c(norm(BPPPmean-datawt$Sigma0,"2"),norm(BPPPmean-datawt$Sigma0,"F"))

c(norm(bandingEST-datawt$Sigma0,"2"),norm(bandingEST-datawt$Sigma0,"F"))


##밴드 선택하는 부분

#cross validation
source("simulation/bandselection.R")

bandres2=bandselection(datawt$X,kvec=1:10,n.cv=100,normtype = "2")
bandresF=bandselection(datawt$X,kvec=1:10,n.cv=100,normtype = "F")

plot(bandres$BPPP)
plot(bandres$BandEST)

plot(bandresF$BPPP)
plot(bandresF$BandEST)
