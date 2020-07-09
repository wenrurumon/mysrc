getnsteps <- function(v1,v2,g){length(igraph::shortest_paths(g,v1,v2)$vpath[[1]])}
global.test <- function(DAG,var.discrete=NULL,n=10,rate=0.8,replace=T,unobserved=T){
  g <- igraph::graph_from_adjacency_matrix(DAG$DAG)
  g.achieve <- lapply(V(g),function(v1){
    rlt <- sapply(V(g),getnsteps,v1=v1,g=g)
    names(rlt[rlt>1])
  })
  g.achieve <- g.achieve[sapply(g.achieve,length)>0]
  map2test <- do.call(rbind,lapply(1:length(g.achieve),function(i){
    data.table(v1=names(g.achieve)[i],v2=g.achieve[[i]]) %>% mutate(method=(v1%in%var.discrete)|(v2%in%var.discrete))
  }))
  test <- lapply(1:nrow(map2test),function(i){
    f <- formula(paste0(map2test[i,2],'~',map2test[i,1]))
    # print(paste(i,Sys.time()))
    if(map2test[i,3]){
      rlt <- causal.infer.discrete(f,as.data.frame(data),DAG,rate,n,replace,unobserved)
    } else {
      rlt <- causal.infer.continue(f,as.data.frame(data),DAG,unobserved)
    }
  })
  test
}
