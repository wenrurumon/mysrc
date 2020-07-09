
rm(list=ls())
devtools::install_github("wenrurumon/mysrc/ce.test.v1",force=T)
library(ce.test.v1)

data <- read.csv2("ukb_t2d_0702.csv",header=T,sep=',')
data <- data[,!colnames(data)%in%c('Oestradiol','Rheumatoid_factor')]
data <- sapply(na.omit(data),function(x){as.numeric(paste(x))})
colnames(data) <- tolower(colnames(data))

system.time(DAG <- init.DAG(data=data,n=10,rate=0.8,alpha=0.05,replace=F))
plot.dag(DAG$DAG)

causal.infer.continue(eid~bmi0,as.data.frame(data),DAG,unobserved=F) 
causal.infer.discrete(albumin~coffee_intake0,as.data.frame(data),DAG,0.8,100,T,F)

ce.func('eid','dm',DAG,F)
ce.func('pork_intake0','beef_intake0',DAG,F)
ce.func('eid','bmi0',DAG,F)
ce.func('albumin','coffee_intake0',DAG,F)
