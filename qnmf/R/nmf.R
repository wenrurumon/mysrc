#NMF
nmf2 <- function(A,K){
  model <- NMF::nmf(A,K)
  X <- basis(model)
  Y <- coef(model)
  A <- X %*% Y
  list(A=A,X=X,Y=Y)
}

qnmf_wx <- function(A,X,lambda,a,maxitn){
  Y <- positive(ginv(t(X)%*% X) %*% t(X) %*% A)
  i <- 0
  while(TRUE){
    if(i>=maxitn){break}
    i <- i+1
    Y2 <- positive(ginv(t(X)%*% X) %*% t(X) %*% A - lambda*max(Y))
  }
  X2 <- X
  rlt <- (list(A=X2%*%Y2,X=X2[,],Y=Y2,itn=i))
  rlt
}

qnmf_wox <- function(A,K,lambda,a,maxitn){
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
  #output
  rlt <- (list(A=X2%*%Y2,X=X2[,],Y=Y2,itn=i))
  rlt
}

qnmf <- function(A,K=3,lambda=0.1,a=0.5,maxitn=1000,X=NULL){
  if(is.null(K)&is.null(X)){"either X or K has to be in the model")}
  if(is.null(X)){
    qnmf_wox(A,K,lambda,a,maxitn)
  } else {
    qnmf_wx(A,X,lambda,a,maxitn)
  }
}
