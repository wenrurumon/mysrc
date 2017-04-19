library(GenABEL)
qpca <- function(A,rank=0){
  A <- A[,apply(A,2,var)>0,drop=F]
  A <- scale(A)[,]
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

qpca2 <- function(A,ifscale=T,l1=0.9,l2=0.9){
  A.pca <- qpca(A)
  A.qpca <- qpca(A,rank=which(A.pca$prop>=l1)[1])
  A.qpca$X[,1:which(A.qpca$prop>=l2)[1],drop=F]
}
