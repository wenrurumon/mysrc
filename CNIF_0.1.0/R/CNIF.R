library(Rcpp)
library(RcppArmadillo)
library(igraph)
library(hash)
library(Rsymphony)
require(rstackdeque)

#Rcpp::sourceCpp("./src/initial_sem.cpp")
#Rcpp::sourceCpp("./src/score_function_regression.cpp")

CNIF <- function(data,level=NULL, init.adj=NULL,max_parent =5, score_function=score_function_regression,
                 max_cycle_depth=5,lambda=0.001,rho=0.5){

  data = data.matrix(data)

  if( is.null(init.adj) ) {
    init.adj =  (initial_sem_nofix(data,1:ncol(data),lambda=lambda,rho=rho) != 0 )
  }

  if( is.null(level) ){
    level = rep(1,ncol(data))
  }

  colnames( init.adj ) = colnames(data)
  rownames( init.adj ) = colnames(data)

  scores.list = load.score( init.adj, data, max_parent, level,score_function )

  rlt = init.adj
  rlt[,] = 0

  for( tmp in scores.list){
    vets = tmp$vets
    score.frame = tmp$score.frame
    parms = length(vets)


    obj = score.frame[,parms+2]
    types = rep("I",parms)


    rhs   = hash()
    sense   = hash()
    A     = hash()

    working.node = score.frame[,parms+1]
    working.parent = score.frame[,1:parms]
    for( i in 1:parms ){
      cons.name = paste("one parent set",i)
      rhs[[cons.name]] = 1
      sense[[cons.name]] = "=="
      A[[cons.name]]   = (working.node == i)
    }


    cons.names = keys(A)

    A = t(values(A,cons.names))
    sense = values(sense,cons.names)
    rhs = values(rhs,cons.names)

    csc.A = make.constraints(tmp$tri.mat,score.frame )

    A = rbind(A, csc.A)
    sense = c(sense, rep("<=",nrow(csc.A)))
    rhs = c(rhs, rowSums(csc.A)-1)

    result = Rsymphony_solve_LP(  obj=obj,mat=A,dir=sense,rhs=rhs,types=types,verbosity = 0  )
    solution = score.frame[result$solution == 1,]

    for( i in 1:nrow(solution)){
      rlt[   vets[ solution[i,parms+1] ],
             vets[ as.logical( solution[i,1:parms]) ]
             ] = 1
    }
  }
  rlt
}

#data =read.table("C:\\Users\\winon\\Downloads\\zhu\\zhu\\pheno.txt")

