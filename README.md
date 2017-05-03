# mysrc

devtools::install_github("wenrurumon/mysrc/cca",force=T)<br />
devtools::install_github("wenrurumon/mysrc/lrm",force=T)<br />
devtools::install_github("wenrurumon/mysrc/qnmf",force=T)<br />


'''bash
library(qnmf)
library(NMF)
library(corpcor)
#create random matrix
rmat <- matrix(runif(1500,10,20),50,30)
rmat.nmf <- nmf(rmat,5)
rmat.nmf <- basis(rmat.nmf) %*% coef(rmat.nmf)
rmat.qnmf <- qnmf(rmat,5)$A
mean((rmat-rmat.nmf)^2);mean((rmat-rmat.qnmf)^2)
#create random matrix for deconvolution simulation
rmatx <- matrix(runif(150,10,20),50,3)
rmaty <- matrix(runif(60,5,10),3,20)
rmat <- rmatx %*% rmaty
rmat.nmf <- nmf(rmat,3)
rmat.nmf <- basis(rmat.nmf) %*% coef(rmat.nmf)
rmat.qnmf <- qnmf(rmat,3)$A
mean((rmat-rmat.nmf)^2);mean((rmat-rmat.qnmf)^2)
#create random parse matrix for deconvolution simulation
rmatx <- matrix(runif(150,10,20),50,3)
rmaty <- matrix(runif(60,5,10),3,20)
rmat <- rmatx %*% rmaty
rmat.nmf <- nmf(rmat,3)
rmat.nmf <- basis(rmat.nmf) %*% coef(rmat.nmf)
rmat.qnmf <- qnmf(rmat,3)$A
mean((rmat-rmat.nmf)^2);mean((rmat-rmat.qnmf)^2)
'''
