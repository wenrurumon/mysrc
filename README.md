# mysrc

devtools::install_github("wenrurumon/mysrc/cca",force=T)<br />
devtools::install_github("wenrurumon/mysrc/lrm",force=T)<br />
devtools::install_github("wenrurumon/mysrc/qnmf",force=T)<br />


```bash
library(qnmf)
library(NMF)
library(corpcor)
testi <- function(){
  rmatx <- matrix(runif(150,10,20),50,3)
  rmaty <- matrix(runif(60,5,10),3,20)
  rmaty[sample(rmaty,length(rmaty)/3)] <- 0
  rmat_raw <- rmat <- rmatx %*% rmaty
  error <- matrix(rnorm(length(rmat),mean=0,sd=sd(as.vector(rmat))),ncol=ncol(rmat),nrow=nrow(rmat))
  rmat <- rmat + error
  rmat <- rmat * (rmat>0)
  rmat.nmf <- nmf(rmat,3)
  rmat.nmf <- basis(rmat.nmf) %*% coef(rmat.nmf)
  rmat.qnmf <- qnmf(rmat,3,0.4,0.5)$A
  m1 <- mean((rmat_raw-rmat.nmf)^2)
  m2 <- mean((rmat_raw-rmat.qnmf)^2)
  c(m1,m2)
}
rmaxs <- sapply(1:1000,function(i){
  try(testi())
})
summary(t(rmaxs))
t.test(rmaxs[1,],rmaxs[2,])
```
