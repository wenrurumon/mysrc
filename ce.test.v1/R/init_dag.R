
plot.dag <- function(x){
  plot(graph.adjacency((x)),
       edge.arrow.size=.5,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=1)
}

score <- function(nodei,DAG,data,add.parent=NULL){
  data <- as.matrix(data)
  nodes.X <- parents(DAG,nodei)
  X <- data[,c(nodes.X, add.parent),drop=F]
  Y <- data[,nodei,drop=F]
  s <- s.Y <- (t(Y)%*%(diag(nrow(data)))%*%Y)
  if(length(nodes.X)>0){
    s <- t(Y) %*% (diag(nrow(data)) - X%*%solve(t(X)%*%X)%*%t(X)) %*% Y
  }
  s.delta <- ifelse(length(nodes.X)==0,0,(s.Y-s)/length(nodes.X))
  c(score=s,delta_score=s.delta,delta_score_frac=s.delta/s.Y,n_parents=length(nodes.X))
}

scores <- function(data,DAG){
  rlt <- do.call(rbind,lapply(colnames(data),score,DAG=DAG,data=data,add.parent=NULL))
  rownames(rlt) <- colnames(data)
  rlt
}

sample.score <- function(n,rate,data,DAG){
  sapply(1:n,function(i){
    s <- scores(data[sample(nrow(data),rate*nrow(data),replace=F),],DAG)
    sum(s[,3]*s[,4],na.rm=T)/sum(s[,4],na.rm=T)
  })
}

add.pa.point <- function(DAG,data,CDF,alpha=0.05){
  map.ori <- cbind(apply(DAG$arcs,2,function(x){match(x,colnames(data))}),score=NA,pvalue=NA)
  map.new <- do.call(rbind,lapply(colnames(data),function(i){
    nodes.p <- parents(DAG,i)
    nodes.c <- children(DAG,i)
    do.call(rbind,lapply(nodes.c,function(j){
      d.score.frac <- 1-(score(i,DAG,data,j)/score(i,DAG,data))[1]
      p.value <- mean(CDF > d.score.frac)
      c(from=which(colnames(data)==j),to=which(colnames(data)==i),d.score.frac,pvalue=p.value)
    }))
  }))
  map <- as.data.frame(rbind(map.ori,map.new[map.new[,4]<alpha,])) %>% arrange(from,to)
  rlt <- sapply(1:ncol(data),function(i){(1:ncol(data)) %in% map[map[,1]==i,2]}) + 0
  dimnames(rlt) <- list(colnames(data),colnames(data))
  list(mix.matrix=t(rlt),map=map)
}

dag2mat <- function(DAG,data){
  map.DAG <- apply(DAG$arcs,2,function(x){match(x,colnames(data))})
  out <- sapply(1:ncol(data),function(i){
    (1:ncol(data) %in% map.DAG[map.DAG[,1]==i,2])
  })
  dimnames(out) <- list(colnames(data),colnames(data))
  t(out)
}

init.DAG <- function(data, n=10, rate=0.8, alpha = 0.05, replace = F){
  DAG <- rsmax2(as.data.frame(data)) 
  trait.score <- scores(data,DAG)
  CDF <- sample.score(n,rate,data,DAG)
  CCM <- add.pa.point(DAG,data,CDF,alpha)
  acm <- dag2mat(DAG, data)
  arrow <- apply(CCM$map,1,function(x){paste(colnames(data)[x[1]],colnames(data)[x[2]],sep='-+')})
  un.path <- apply(CCM$map,1,function(x){paste(colnames(data)[x[2]],colnames(data)[x[1]],sep='-+')})
  un.path <- as.vector(t(cbind(un.path[un.path%in%arrow],arrow[un.path%in%arrow])))
  index <- match(unique(un.path),arrow)
  arrow.ccm <- paste(arrow,collapse=',')
  arrow.ori <- paste(arrow[is.na(CCM$map[,3])],collapse=',')
  list(DAG=acm,trait.score=trait.score,CDF=CDF,mix.DAG=CCM$mix.matrix,map=CCM$map,
       formula=arrow.ori,formula.unobserved=arrow.ccm,unobserved.path.index=index)
}

ce.func <- function(y,x,DAG,unobserved=T){
  if(unobserved){
    dag <- DAG$mix.DAG
    dag[,colnames(dag)==x] <- F
    g <- set.edge.attribute(graph=graph.adjacency(dag),
                            name='description',value='U',index=DAG$unobserved.path.index)
  } else {
    dag <- DAG$DAG
    dag[,colnames(dag)==x] <- F
    g <- set.edge.attribute(graph=graph.adjacency(dag),name='description',value='U',index=NULL)
  }
  causal.effect(y=y,x=x,z=NULL,G=g,expr=T)
}

causal.infer.continue <- function(f, data, DAG, unobserved = T){
  y <- paste(f[2])
  x <- paste(f[3])
  ce <- ce.func(y, x, DAG, unobserved)
  pl <- do.calculus(ce, as.data.frame(data))
  out <- continue.var.test(y,x,pl)
  list(fomula=(f),latex=ce,rlt=out)
}

causal.infer.discrete <- function(f, data, DAG, rate, n, replace = F, unobserved = T){
  y <- paste(f[2])
  x <- paste(f[3])
  ce <- ce.func(y, x, DAG, unobserved = unobserved)
  out <- entrop.test(x, ce, data, rate, n, replace = replace)
  list(fomula=(f),latex=ce,rlt=out)
}
