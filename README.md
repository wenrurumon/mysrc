# mysrc

devtools::install_github("wenrurumon/mysrc/xmndata",force=T)<br />
devtools::install_github("wenrurumon/mysrc/qfcca",force=T)<br />

```bash
rm(list=ls())
library(qfcca)
library(GenABEL)
library(flare)
library(fda)
library(corpcor)

A <- iris[,1:4]
B <- cbind(as.numeric(iris[,5]))
cca(A,B)

A.qpca <- qpca(iris[,1:4],rank=which(qpca(A)$prop>=0.8)[1])$X
cca(A.qpca,B)

A.fpca <- fpca(A,1:4,7,0)
A.fpca <- lapply(A.fpca,function(x){
  x$score[,which(x$prop>=0.8)[1],drop=F]
})
sapply(A.fpca,cca,B=B)
```

devtools::install_github("wenrurumon/mysrc/cca",force=T)<br />
devtools::install_github("wenrurumon/mysrc/lrm",force=T)<br />
devtools::install_github("wenrurumon/mysrc/qnmf",force=T)<br />


```bash
rm(list=ls())
library(e1071)
library(parallel)
library(preprocessCore)
library(qnmf)
library(NMF)
library(corpcor)
test <- function(nrow,ncol,k,sparse,error){
  # r <- rmat(100,20,5,.5,1,T)
  r <- rmat(nrow,ncol,k,sparse,error,T)
  rlt.qnmf <- qnmf(r$A,K=5)
  rlt.nmf <- nmf2(r$A,5)
  rlt.qnmfx <- qnmf(r$A,X=r$X,deconv=T)
  rlt.stf <- stf.deconv(r$A,x=r$X)
  sapply(list(rlt.qnmf$A,rlt.nmf$A,rlt.qnmfx$A,rlt.stf$A),fita,A=r$raw)
  # cbind(om=diag(cor(t(r$Y),t(rlt.qnmfx$Y))),stf=diag(cor(t(r$Y),t(rlt.stf$Y))))
}
simu <- lapply(1:1000,function(i){
  test(100,20,5,.5,1)
})
```

devtools::install_github("wenrurumon/mysrc/jointnet",force=T)<br />

```bash
rm(list=ls())
library(jointnet)
library(corpcor)
library(abind)
raw <- jn_example()
jointnet(raw)
```

