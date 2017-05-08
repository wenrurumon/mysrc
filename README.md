# mysrc

devtools::install_github("wenrurumon/mysrc/cca",force=T)<br />
devtools::install_github("wenrurumon/mysrc/lrm",force=T)<br />
devtools::install_github("wenrurumon/mysrc/qnmf",force=T)<br />


```bash
library(e1071)
library(parallel)
library(preprocessCore)
library(qnmf)
library(NMF)
library(corpcor)

r <- rmat(100,20,5,.5,.3,T)
rlt.qnmf <- qnmf(r$A,5)
rlt.nmf <- nmf2(r$A,5)
rlt.qnmfx <- qnmf(r$A,X=r$X,deconv=T)
rlt.stf <- stf.deconv(r$A,x=r$X)
sapply(list(rlt.qnmf$A,rlt.nmf$A,rlt.qnmfx$A,rlt.stf$A),fita,A=r$A)
cbind(om=diag(cor(t(r$Y),t(rlt.qnmfx$Y))),stf=diag(cor(t(r$Y),t(rlt.stf$Y))))
```
