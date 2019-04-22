##Khare 방법과 BPPP 비교

##main 함수
n=100
p=200
k=5

##참값과 자료를 만드는 부분
source("simulation/simulsetting.R")
datawt=simulsetting(n=n,p=p,k=k)

##MCMC의 초기값 설정(참값으로 설정한다)

choldecomp=chol(datawt$Sigma0)
D=diag((diag(choldecomp))^2,p)
L=t(solve(sqrt(D))%*%choldecomp)

source("simulation/band_Khare.R")
kharesample=band_Khare(datawt$X,k,L=L,D=D)

source("simulation/BPPP.R")
bPPPsample=bPPP(datawt$X,bandwidth=k,adj=TRUE)

BPPPmean=Reduce("+",bPPPsample$BPPP)/length(bPPPsample$BPPP)  
kharemean=Reduce("+",kharesample)/length(kharesample)  

c(norm(BPPPmean-datawt$Sigma0,"2"),norm(BPPPmean-datawt$Sigma0,"F"))
c(norm(kharemean-datawt$Sigma0,"2"),norm(kharemean-datawt$Sigma0,"F"))

S=t(datawt$X)%*%datawt$X/n
bandingEST=adjust_pd(banding(S,k),outlist=FALSE)
bandingEST=banding(S,k)

c(norm(bandingEST-datawt$Sigma0,"2"),norm(bandingEST-datawt$Sigma0,"F"))
