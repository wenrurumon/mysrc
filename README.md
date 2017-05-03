# mysrc

devtools::install_github("wenrurumon/mysrc/cca",force=T)<br />
devtools::install_github("wenrurumon/mysrc/lrm",force=T)<br />
devtools::install_github("wenrurumon/mysrc/qnmf",force=T)<br />


```bash
# devtools::install_github("wenrurumon/mysrc/qnmf",force=T)
library(qnmf)
library(NMF)
library(corpcor)
```
```bash
#create random matrix
rmat <- matrix(runif(1500,10,20),50,30)
rmat.nmf <- nmf(rmat,5)
rmat.nmf <- basis(rmat.nmf) %*% coef(rmat.nmf)
rmat.qnmf <- qnmf(rmat,5)$A
mean((rmat-rmat.nmf)^2);mean((rmat-rmat.qnmf)^2)
```
```bash
#create random matrix for deconvolution simulation
rmatx <- matrix(runif(150,10,20),50,3)
rmaty <- matrix(runif(60,5,10),3,20)
rmat <- rmatx %*% rmaty
error <- matrix(rnorm(length(rmat),mean=0,sd=sd(as.vector(rmat))),ncol=ncol(rmat),nrow=nrow(rmat))
rmat <- rmat + error
rmat.nmf <- nmf(rmat,3)
rmat.nmf <- basis(rmat.nmf) %*% coef(rmat.nmf)
rmat.qnmf <- qnmf(rmat,3)$A
mean((rmat-rmat.nmf)^2);mean((rmat-rmat.qnmf)^2)
```
```bash
#simulation for sparsed deconvolution
testi <- function(){
  rmatx <- matrix(runif(150,10,20),50,3)
  rmaty <- matrix(runif(60,5,10),3,20)
  rmat_raw <- rmat <- rmatx %*% rmaty
  error <- matrix(rnorm(length(rmat),mean=0,sd=sd(as.vector(rmat))),ncol=ncol(rmat),nrow=nrow(rmat))
  rmat <- rmat + error
  rmat <- ifelse(rmat>0,rmat,0)
  rmat.nmf <- nmf(rmat,3)
  rmat.nmf <- basis(rmat.nmf) %*% coef(rmat.nmf)
  rmat.qnmf <- qnmf(rmat,3)$A
  m1 <- mean((rmat_raw-rmat.nmf)^2)
  m2 <- mean((rmat_raw-rmat.qnmf)^2)
  # print(paste(m1,m2,m1>m2))
  c(m1,m2)
}
rmaxs <- sapply(1:1000,function(i){
  print(i)
  try(testi())
})
```
