#Pcca
pcca <- function(x,...){UseMethod('pcca')}
pcca.list <- function(...,prop=0.99,ifscale=T){
  tmp <- list(...)
  if(length(tmp)==2){
    x1 <- tmp[[1]]
    x2 <- tmp[[2]]
  } else {
    x1 <- x2 <- tmp[[1]]
  }
  x1 <- pca(x1)
  x2 <- pca(x2)
  cca(x1,x2)
}

#Cluster
fc <- function(x){
  w<-as.vector(t(x))[t(x)>0]
  x <- igraph::graph_from_adjacency_matrix(x>0,mode='undirected')
  fc <- igraph::membership(igraph::fastgreedy.community(x,weight=w))
  fc[] <- match(fc,unique(fc))
  fc
}
fc2 <- function(x,thred=NULL){
  if(is.null(thred)){
    thred <- 0.01/length(x)
  }
  x.mat <- (x<thred)+0
  diag(x.mat) <- 0
  x.score <- -log(x)
  x.score[x.mat==0] <- 0
  x.score[x.score==Inf] <- max(x.score[x.score!=Inf]*2)
  x.score <- x.score/max(x.score)
  x.g <- igraph::graph_from_adjacency_matrix(x.mat,mode='directed')
  igraph::E(x.g)$weight <- as.vector(x.score)[x.mat>0]
  x.g <- igraph::as.undirected(x.g)
  list(network = x.mat,
       cluster = igraph::fastgreedy.community(x.g)$membership)
}
plotclust <- function(x,membership=NULL,main=NULL){
  G <- igraph::graph_from_adjacency_matrix(x>0,mode='undirected')
  if(is.null(membership)){membership=rep(1,ncol(x))}
  plot(igraph::create.communities(G, membership),
       igraph::as.undirected(G),
       layout=igraph::layout.kamada.kawai(igraph::as.undirected(G)),
       edge.arrow.size=.1,
       vertex.size=.3,
       vertex.label.cex=1,
       edge.width=.1,
       main=main)
}
