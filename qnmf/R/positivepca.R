
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
scale2 <- function(x){
  apply(x,2,function(xi){
    (xi-min(xi))/max(xi-min(xi))
  })
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
positive <- function(x){
  x * (x>0)
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

qnmf <- function(A,K=3,lambda=0.5,a=0.9,maxitn=1000,X=NULL){
  if(!is.null(X)){
    Y <- positive(ginv(t(X)%*% X) %*% t(X) %*% A)
    i <- 0
    while(TRUE){
      if(i>=maxitn){break}
      i <- i+1
      Y2 <- positive(ginv(t(X)%*% X) %*% t(X) %*% A - lambda*max(Y))
    }
    X2 <- X
  } else {
    #initialization
    A.svd <- svd(qpca(A,F,K+1)$Z)
    # A.svd <- svd(scale(A))
    U <- A.svd$u[,1:K,drop=F]
    V <- A.svd$v[,1:K,drop=F]
    D <- diag2(A.svd$d[1:K])
    for(i in 1:K){
      if(max(U[,i]) < 0){
        U[,i] <- -U[,i]
        V[,i] <- -V[,i]
      }
    }
    X <- U %*% sqrt(D)
    Y <- sqrt(D) %*% t(V)
    #Loops
    i <- 0
    X <- positive(X)
    # X <- scale2(X); Y <- scale2(Y)
    while(TRUE){
      i <- i+1
      if(i>=maxitn){break}
      Y2 <- positive(ginv(t(X)%*% X) %*% t(X) %*% A - lambda*max(Y))
      X2 <- positive(A %*% t(Y2) %*% ginv(Y2 %*% t(Y2)))
      # X2 <- scale2(X2); Y2 <- scale2(Y2)
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
  }

  #output
  rlt <- (list(A=X2%*%Y2,X=X2[,],Y=Y2,itn=i))
  rlt
  # diag(cor(A,rlt$A))
}

