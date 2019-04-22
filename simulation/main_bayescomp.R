##main 함수
n=100
p=90
k=10

##참값과 자료를 만드는 부분
source("simulation/simulsetting.R")
datawt=simulsetting(n=n,p=p,k=k)

##MCMC의 초기값 설정(참값으로 설정한다)

choldecomp=chol(datawt$Sigma0)
D=diag((diag(choldecomp))^2,p)
L=t(solve(sqrt(D))%*%choldecomp)

source("simulation/band_Khare.R")
kharesample=band_Khare(datawt$X,k,L=L,D=D)
