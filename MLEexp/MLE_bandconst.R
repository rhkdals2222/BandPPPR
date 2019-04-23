library(CVXR)
library(Rfast)
library(FinCovRegularization)


MLE_bandconst=function(X,k,Tlen=1000,epsilon=10^(-4),tol=epsilon){
  
  n=dim(X)[1]
  p=dim(X)[2]
  S=mat.mult(t(X),X)/n
  
  Sigmat=banding(S,k)
  
  if(min(eigen(as.matrix(Sigmat))$value)<epsilon){
    Sigmat=as.matrix(Sigmat + diag(-min(eigen(as.matrix(Sigmat))$value)+epsilon,dim(Sigmat)[1]))
  }
  

  emS=max(eigen(S)$value)
  statvec=rep(0,Tlen)
  Sigma <- Semidef(p)
  for(t in 1:Tlen){
    
    ##기울기 계산
    Sigmatinv=spdinv(Sigmat)
    
    der=-0.5*(Sigmatinv-mat.mult(mat.mult(Sigmatinv,S),Sigmatinv))
    obj <- Maximize(sum(der*Sigma))
    ##밴드 바깥에 있는 원소를 지칭하는 행렬
    const_mtx=1-banding(matrix(1,p,p),k)
    
    ##밴드 바깥의 원소의 절대값의 합이 tol보다 작다는 조건과 대각행렬의 합이 너무 크지 않음을 나타냄.
    constraints <- list(sum(abs(Sigma*const_mtx)) <= tol,sum(abs(Sigma*diag(1,p)))<=10*p*emS)
    #constraints <- list(sum(abs(Sigma*const_mtx)) <= tol)
    
    ## Form and solve optimization problem
    prob <- Problem(obj, constraints)
    result2 <- solve(prob)
    statvec[t]=result2$status
    gamma=2/(2+t)
    Sigmat=as.matrix((1-gamma)*Sigmat+gamma*result2$getValue(Sigma))
  }
  return(list(Sigma=Sigmat,status_opt=statvec))
}


