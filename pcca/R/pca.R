#PCA
pca <- function(x,...){UseMethod('pca')}
pca1 <- function(A,rank=0,ifscale=TRUE){
  A <- A[,colMeans(is.na(A))<1,drop=F]
  A <- apply(A,2,function(x){
    x[is.na(x)] <- median(x,na.rm=T)
    x
  })
  if(ncol(A)>0){
    A <- A[,apply(A,2,var)>0,drop=F]
  }
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0|rank==ncol(A)){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
    d <- d[d > 1e-8]
  }
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}
pca.matrix <- function(A,prop=0.99,ifscale=TRUE){
  rlt <- pca1(A,rank=0,ifscale=ifscale)
  if(prop==1){
    return(rlt)
  }
  rlt <- pca1(A,rank=which(rlt$prop>=prop)[1],ifscale=TRUE)
  return(rlt$X)
}
pca.list <- function(X,prop=0.99,ifscale=T){
  lapply(X,pca,ifscale=ifscale,prop=prop)
}
