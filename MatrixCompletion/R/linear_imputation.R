imlm <- function(y,x){
  #Linear prediction from x to y
  sel.lm <- coef(lm(y~x))
  sel.lm[1]+sel.lm[2]*x
}
imcor <- function(i,j,test=F){
  #COR(I,J,NA.RM=T)
  #if test then the significance level would be provided
  sel <- !(is.na(i))&!(is.na(j))
  if(test){
    cor.test(i[sel],j[sel])$p.value
  } else {
    cor(i[sel],j[sel])
  }
}
imlms <- function(y,xs){
  #INPUTE data with linear combination of specific dependent variables
  temp <- apply(xs,2,function(x){
    imlm(y,x)
  })
  rowMeans(temp,na.rm=T)
}
qpca <- function(A,rank=0){
  #INPUTE data with QPCA
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
