#Rmatrix

rmat <- function(nrow,ncol,k,sparse=0.3,error=1){
  x <- matrix(runif(nrow*k,0,10),nrow,k)
  y <- matrix(runif(k*ncol,0,10),k,ncol)
  y[sample(length(y),length(y)*sparse)] <- 0
  a <- x %*% y
  error <- rnorm(length(a),mean=0,sd=sd(as.vector(a))) * error
  a <- a + error
  a <- a * (a>0)
  list(A=a,X=x,Y=y)
}

fita <- function(A,a){
  e <- (A-a)
  sum(e^2)/sum(A^2)
}
