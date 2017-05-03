
rm(list=ls())

####################################################
# Macro
####################################################

diag2 <- function(x){
  if(length(x)==1){
    return(as.matrix(x))
  } else {
    return(diag(x))
  }
}
ginv<-function(A){
  A_svd<-fast.svd(A)
  if(length(A_svd$d)==1){
    A_inv<-A_svd$v%*%as.matrix(1/A_svd$d)%*%t(A_svd$u)
  }else{
    A_inv<-A_svd$v%*%diag(1/A_svd$d)%*%t(A_svd$u)
  }
  return(A_inv)
}
positive <- function(x,turn=F){
  if(turn){
    for(i in 1:ncol(x)){
      x[,i] <- x[,i] * sign(max(x[,i]))
    }
  }
  ifelse(x>0,x,0)
}
qpca <- function(A,scale=T,rank=0){
  if(scale){A <- scale(A)}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-10]
  r <- length(d)
  prop <- d^2; prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop)
  return(rlt)
}

qnmf <- function(A,K=3,lambda=100,a=0.9,maxitn=1000){
  #adjustment?
  K <- K+1
  #initialization
    A.svd <- svd(qpca(A,F,K+1)$Z)
    # A.svd <- svd(scale(A))
    U <- A.svd$u[,1:K,drop=F]
    V <- A.svd$v[,1:K,drop=F]
    D <- diag2(A.svd$d[1:K])
  
    X <- U %*% sqrt(D)
    for(i in 1:ncol(X)){
      X[,i] <- X[,i] * sign(max(X[,i]))
    }
    Y <- sqrt(D) %*% t(V)

  #Loops
    i <- 0
    while(TRUE){
      #print(i <- i+1)
      if(i>=maxitn){break}
      Iy <- (Y<0)
      Y2 <- positive(ginv(t(X)%*% X) %*% t(X) %*% A - lambda * Iy)
      X2 <- positive(A %*% t(Y2) %*% ginv(Y2 %*% t(Y2)))
      Xf <- matrixcalc::frobenius.norm(X2-X)
      Yf <- matrixcalc::frobenius.norm(Y2-Y)
      X <- X2; Y <- Y2  
      lambda <- lambda * a
      if(Xf<=1e-8 & Yf<=1e-8){break}
    }
  #Finalize
    X <-X[,apply(X,2,var)>0]
    Y2 <- positive(ginv(t(X)%*% X) %*% t(X) %*% A)
    X2 <- positive(A %*% t(Y2) %*% ginv(Y2 %*% t(Y2)))
    
  #output
  return(list(A=X2%*%Y2,X=X2[,],Y=Y2,itn=i))
}
