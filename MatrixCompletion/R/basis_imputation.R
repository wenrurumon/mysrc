imspline <- function(x,y,pos=max(x),thred=0.1,ifplot=F,ifscale=T){
#INPUTE data with spline basis
  pos <- ceiling(pos)
  if(ifscale){
    x <- x/pos * 100
    pos <- 100
  }
  if(length(x)==0){
    return(rep(NA,pos+1))
  }
  if(length(x)==1){
    return(rep(x,pos+1))
  }
  f <- splinefun(x,y,method='monoH.FC')
  f <- f(0:pos)
  f <- ifelse(f>max(y)*(1+thred),max(y)*(1+thred),f)
  f <- ifelse(f<min(y)*(1-thred),min(y)*(1-thred),f)
  if(ifplot){
    plot((0:pos),f,type='l',col=2)
    lines(x,y,type='p')
  }
  f
}
imfourier <- function(x,y,pos=max(x),ifplot=F){
#INPUTE data with fourier basis
  y <- y[order(x)]
  x <- x[order(x)]
  pos <- ceiling(pos)
  fbasis <- create.fourier.basis(range(x)/pos,length(x)*2-1)
  phi <- eval.basis(x/pos,fbasis)
  fcoef <- myinv(t(phi)%*%phi)%*%t(phi)%*%cbind(x/pos)
  phi2 <- eval.basis((0:pos)/pos,fbasis)
  f <- phi2 %*% fcoef
  if(ifplot){
    plot((0:pos)/pos,f,type='l',col=2)
    lines(x/pos,y,type='p')
  }
  f
}
myinv<-function(A){
#MYINV
  A_svd<-fast.svd(A)
  if(length(A_svd$d)==1){
    A_inv<-A_svd$v%*%as.matrix(1/A_svd$d)%*%t(A_svd$u)
  }else{
    A_inv<-A_svd$v%*%diag(1/A_svd$d)%*%t(A_svd$u)
  }
  return(A_inv)
}
