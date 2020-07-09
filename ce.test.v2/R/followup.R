
P.cond <- function(name, value, data, ld = F){
  D <- data[,name]
  if (length(name) == 1){
    D.tar <- D[which(D == value[1])]
    P <- length(D.tar)/length(D) 
  }
  else{
    for (i in rev(c(2:length(name)))) {
      na <- name[i]
      D <- D[which(D[,which(names(D) == na)] == value[i]),]
    }
    D.tar <- D[which(D[,which(names(D) == name[1])] == value[1]),]
    P <- nrow(D.tar)/nrow(D)
    
  }
  if (ld == T){
    return(list(
      P = P,
      D = D
    ))
  }
  else{
    return(
      P = P
    )
  }
}

P.do.base <- function(var.list, data){
  P.cond <- function(name, value, data, ld = F){
    D <- data[,name]
    if (length(name) == 1){
      D.tar <- D[which(D == value[1])]
      P <- length(D.tar)/length(D) 
    }
    else{
      for (i in rev(c(2:length(name)))) {
        na <- name[i]
        D <- D[which(D[,which(names(D) == na)] == value[i]),]
      }
      D.tar <- D[which(D[,which(names(D) == name[1])] == value[1]),]
      P <- nrow(D.tar)/nrow(D)
      
    }
    if (ld == T){
      return(list(
        P = P,
        D = D
      ))
    }
    else{
      return(
        P = P
      )
    }
  }
  Pvariable <- var.list$P.var
  Sumvariable <- var.list$sum.var
  PPP <- c()
  unique.Ar <- as.data.frame(unique(data[,unique(unlist(Pvariable))]))
  names(unique.Ar) <- unique(unlist(Pvariable))
  for (i in c(1:nrow(unique.Ar))) {
    M <- 1
    for (j in c(1:length(Pvariable))) {
      K <- P.cond(unlist(Pvariable[j]), as.matrix(unique.Ar[i,unlist(Pvariable[j])]), data)
      M <- M*K
    }
    PPP[i] <- M
  }
  P.Ar <- cbind(unique.Ar,PPP)
  pl <- P.Ar
  
  return(P.Ar)
  
}

dis.col <- function(pl, var.list, data){
  si <- unlist(var.list$sum.var)
  si.index <- which(names(pl) %in% si)
  PPP.index <- ncol(pl)
  pl.sum <- pl[,-c(si.index, PPP.index)]
  pl.sum.name <- names(pl)[which(!(c(1:ncol(pl)) %in% c(si.index, PPP.index)))]
  pl.sum <- as.data.frame(pl.sum)
  names(pl.sum) <- pl.sum.name
  eq.ma <- matrix(0, nrow = nrow(pl.sum), ncol = nrow(pl.sum))
  
  
  for (i in c(1:nrow(pl.sum))) {
    for (j in c(1:nrow(pl.sum))) {
      if (sum(pl.sum[i,] == pl.sum[j,]) == ncol(pl.sum)){
        eq.ma[i,j] <- 1
      }
    }
  }
  new.P <- pl$PPP
  pl.new <- cbind(pl.sum, new.P)
  
  for (i in c(1:nrow(pl.new))){
    pl.new[i,ncol(pl.new)] <- sum(new.P[which(eq.ma[i,] == 1)])
  }
  pl.new <- unique(pl.new)
  return(pl.new)
}

dis.col.after.sum <- function(var.list, data){
  
  if(length(var.list$sum.var) == 0){
    pl <- P.do.base(var.list, data)
    out <- dis.col(pl, var.list, data)
    return(out)
  }
  else{
    if(var.list$sum.var == 'devide'){
      pl <- P.do.base(var.list, data)
      pl$PPP <- 1/pl$PPP
      pl <- rename(pl, c('new.P'='PPP'))
      return(pl)
    }
    else{
      pl <- P.do.base(var.list, data)
      out <- dis.col(pl, var.list, data)
      return(out)
    }
  }
}

merge.do <- function(q12, q3, stage.sum){
  #stage.sum <- c('petal.length','petal.width')
  n12 <- names(q12)[which(names(q12) %in% stage.sum)]
  n3 <- names(q3)[which(names(q3) %in% stage.sum)]
  stage.sum <- n12[which(n12 %in% n3)]
  dup1 <- names(q12)[which(!(names(q12) %in% stage.sum))]
  dup2 <- names(q3)[which(!(names(q3) %in% stage.sum))]
  dup.name <- dup1[which(dup1 %in% dup2)]
  dup.name <- dup.name[-which(dup.name =='new.P')]
  stage.sum <- c(stage.sum, dup.name)
  m <- merge(q3, q12, by = stage.sum)
  #m <- merge(d1, d2, by = s.sum)
  new.P <- (m$new.P.x)*(m$new.P.y)
  xy.index <- which(names(m) %in% c('new.P.x', 'new.P.y'))
  pl <- cbind(m[,-xy.index], new.P)
  return(pl)
}

sum.after.merge <- function(pl, stage.sum = NULL, sum = T){
  if (sum == T){
    si <- stage.sum
    si.index <- which(names(pl) %in% si)
    PPP.index <- ncol(pl)
    pl.sum <- pl[,-c(si.index, PPP.index)]
    eq.ma <- matrix(0, nrow = nrow(pl.sum), ncol = nrow(pl.sum))
    for (i in c(1:nrow(pl.sum))) {
      for (j in c(1:nrow(pl.sum))) {
        if (sum(pl.sum[i,] == pl.sum[j,]) == ncol(pl.sum)){
          eq.ma[i,j] <- 1
        }
      }
    }
    new.P <- pl$new.P
    pl.new <- cbind(pl.sum, new.P)
    for (i in c(1:nrow(pl.new))){
      pl.new[i,ncol(pl.new)] <- sum(new.P[which(eq.ma[i,] == 1)])
    }
    pl.new <- unique(pl.new)
    return(pl.new)
  }
  else{
    new.P <- pl$new.P
    new.P <- 1/new.P
    pl$new.P <- new.P
    return(pl)
  }
}

getraw <- function(x,remove=F){
  if(remove){
    return(gsub('\\\\','',x))
  }
  sapply(strsplit(x,''),function(x){
    x[x%in%c('(',')','|','{','}')] <- paste0('\\',x[x%in%c('(',')','|','{','}')])
    paste(x,collapse='')
  })
}

getsyn <- function(x){
  x <- (gsub('\\\\','',tolower(x)))
  p <- getraw(unique(unlist(str_extract_all(x,'p\\(.+?\\)'))))
  for(i in 1:length(p)){x <- gsub(p[i],paste0('(P',i,')'),x)}
  p <- getraw(p,T)
  f <- getraw(unique(unlist(str_extract_all(x,'sum_\\{.+?\\}'))))
  for(i in 1:length(f)){x <- gsub(f[i],paste0('@F',i+2,'@'),x)}
  f <- getraw(f,T)
  x <- gsub('right\\)','\\}',gsub('left\\(','\\{',x))
  frac <- lapply(str_extract_all(x,'frac\\{.+?\\}\\{.+?\\}'),function(x){
    x1 <- getraw(str_extract(x,'frac\\{.+?\\}'),T)
    x2 <- getraw(gsub('frac\\{.+?\\}','',x),T)
    c(getraw(x),paste0(gsub('frac','',x1),'\\{@F1@',x2,'\\}'))
  })
  for(i in 1:length(frac)){x <- gsub(frac[[i]][1],frac[[i]][2],x)}
  list(syn=x,p=p,f=c('devide','',f))
}

splitsyn <- function(x){
  x <- strsplit(x,'')[[1]]
  if(sum(cumsum(x=='{') == cumsum(x=='}'))==1){x <- x[-c(1,length(x))]}
  k <- 1
  j <- 0
  rlt <- list()
  for(i in 1:length(x)){
    if(x[i]%in%c('{')){
      if(j==0){
        k <- k+1
        rlt[k][[1]] <- c(rlt[k][[1]],x[i])
      }
      j <- j+1
    } else if(x[i]%in%c('}')){
      j <- j-1
      if(j==0){
        rlt[k][[1]] <- c(rlt[k][[1]],x[i])
        k <- k+1
      }
    } else {
      rlt[k][[1]] <- c(rlt[k][[1]],x[i])
    }}
  rlt <- rlt[!sapply(rlt,is.null)]
  rlt <- sapply(rlt,paste,collapse='')
  fun <- do.call(rbind,lapply(unlist(str_extract_all(rlt[1],'@.+?@')),function(x){
    c(0,nchar(x))+regexpr(x,rlt[1])
  }))
  if(is.null(fun)){
    rlt <- c('@F2@',rlt)
  } else if(min(fun)!=1){
    rlt <- c('@F2@',rlt)
  } else if(nrow(fun)==1){
    fun <- range(fun)
    rlt <- unlist(list(substr(rlt[1],fun[1],fun[2]-1),substr(rlt[1],fun[2],nchar(rlt[1])),(rlt[-1])))
  } else {
    fun <- range(fun[1,])
    rlt <- unlist(list(substr(rlt[1],fun[1],fun[2]-1),substr(rlt[1],fun[2],nchar(rlt[1])),(rlt[-1])))
  }
  rlt <- rlt[sapply(rlt,function(x){x!=''})]
  rlt <- list(fun=rlt[1],cell=rlt[-1],branch=grepl('@',rlt[-1]))
  rlt
} 

splitsyns <- function(x){
  rlt <- splitsyn(x)
  if(sum(rlt$branch)==0){
    return(rlt)
  } else {
    rlt$cell[rlt$branch] <- lapply(rlt$cell[rlt$branch],splitsyns)
    rlt$cell[!rlt$branch] <- lapply(rlt$cell[!rlt$branch],splitsyn)
  }
  return(rlt)
}

getloop <- function(x,x1){
  x.names <- strsplit(names(x),'\\.')
  x <- data.table(do.call(rbind,lapply(x.names,function(x){
    c(x,rep('NA',max(sapply(x.names,length))-length(x)))
  })),n=sapply(x.names,length),x)
  x <- x[!grepl('branch',sapply(x.names,paste,collapse='')),,drop=F]
  x <- data.table(module=sapply(1:nrow(x),function(i){x[[x$n[i]]][i]}),stage=x$n,value=x$x)
  x$value <- sapply(1:nrow(x),function(i){
    if(x$module[i]=='fun'){
      x1$f[as.numeric(gsub('@F|@','',x$value[i]))]
    } else if(x$module[i]=='cell'){
      paste(
        x1$p[as.numeric(gsub('\\(P|\\)','',str_extract_all(x$value[i],'\\(P.+?\\)')[[1]]))],
        collapse='*'
      )
    }
  })
  x
}

parseletax <- function(x){
  x1 <- getsyn(x)
  x1$f <- gsub('\\}','\\)',gsub('\\{','\\(',x1$f))
  x2 <- (splitsyns(x1$syn))
  getloop(unlist(x2),x1)
}

extract <- function(ce1, sum = T){
  if (ce1 != 'devide'){
    split.str <- unlist(strsplit(ce1,''))
    l.rb <- which(split.str == "(")
    r.rb <- which(split.str == ')')
    rb <- rbind(l.rb, r.rb)
    P.variable <- c()
    for (i in c(1:ncol(rb))){
      P.variable[i] <- substr(ce1, rb[1,i] + 1, rb[2,i]-1)
    }
    P.variable.char <- c()
    for (i in c(1:length(P.variable))) {
      if (("|"%in%unlist(strsplit(P.variable[i],''))) == TRUE){
        unchar <- unlist(strsplit(P.variable[i],''))
        unchar[which(unchar == "|")] <- ','
        new <- NULL
        for (k in c(1:length(unchar))){
          new <- paste(new, unchar[k], sep = "")
        }
        P.variable[i] <- new
        P.variable.char[i] <- strsplit(P.variable[i],',')
      }
      else{
        P.variable.char[i] <- P.variable[i]
      }
    }
    if (sum == T){
      sum.variable <- c()
      
      sum.variable[i] <- substr(ce1, l.rb +1, r.rb-1)
      
      sum.variable.char <- c()
      
      if ((","%in%unlist(strsplit(sum.variable[i],''))) == TRUE){
        sum.variable.char[i] <- strsplit(sum.variable[1],',')
      }
      else{
        sum.variable.char[i] <- sum.variable[i]
      }
      
      return(list(sum.var = sum.variable.char))
    }
    else{
      return(list(P.var = P.variable.char))
    }
  }
  else{
    return(list(sum.var = 'devide'))
  }
  
}

cal.base.parselatex <- function(parse, data){
  max.stage <- max(unique(parse$stage))
  
  if(max.stage > 1){
    DC.var <- NULL
    for (i in rev(c(2:max.stage))){
      subtable <- parse[which(parse$stage == i),]
      if (nrow(subtable)%% 2 != 0){
        subtable <- subtable[-nrow(subtable)]
      }
      DC <- list()
      k <- 1
      for (j in c(1:nrow(subtable))) {
        
        if (subtable[j,1] == 'fun'){
          if (as.character(subtable[j,3]) != ''){
            Sum.list <- extract(as.character(subtable[j,3]), sum = T)
            P.list <- extract(as.character(subtable[j+1,3]), sum = F)
            DC[[k]] <- c(Sum.list, P.list)
            k <- k + 1
          }
          else{
            P.list <- extract(as.character(subtable[j+1,3]), sum = F)
            DC[[k]] <- P.list
            k <- k + 1
          }
          
        }
      }
      
      fun.index <- min(which(parse$stage == i)) - 1
      if (parse[fun.index, 1] == 'fun'){
        if (substr(as.character(parse[fun.index,3]), 1, 3) == 'sum'){
          stage.sum <- as.character(unlist(extract(as.character(parse[fun.index,3]), sum = T)))
          for (b in c(1:(length(DC)))){
            m.dc.2 <- dis.col.after.sum(DC[[b]], data)
            if (length(DC.var) == 0){
              DC.var <- m.dc.2
            }
            else{
              DC.var <- merge.do(DC.var, m.dc.2, stage.sum = stage.sum)
            }
          }
          DC.var <- sum.after.merge(DC.var, stage.sum = stage.sum, sum = T)
        }
        else{
          
          for (b in c(1:(length(DC)))){
            m.dc.2 <- dis.col.after.sum(DC[[b]], data)
            if (length(DC.var) == 0){
              DC.var <- m.dc.2
            }
            else{
              DC.var <- merge.do(DC.var, m.dc.2, stage.sum = stage.sum)
            }
          }
          DC.var <- sum.after.merge(DC.var, stage.sum = stage.sum, sum = F)
        }
        
        
        
      }
    }
  }
  else{
    Sum.list <- extract(as.character(parse[1,3]), sum = T)
    P.list <- extract(as.character(parse[2,3]), sum = F)
    DC <- c(Sum.list, P.list)
    DC.var <- dis.col.after.sum(DC, data)
  }
  
  return(DC.var)
}
do.calculus <- function(ce, data){
  names(data) <- tolower(names(data))
  if (substr(ce, 1,1) == 'P'){
    var.list <- extract(ce, sum = F)
    DC.var <- P.do.base(var.list, data)
    names(DC.var) <- c(names(DC.var)[-length(names(DC.var))], 'new.P')
    
    return(DC.var)
  }
  else{
    parse.ce <- parseletax(ce)
    DC.var <- cal.base.parselatex(parse.ce, data)
    return(DC.var)
  }
}

entropy <- function(X){
  en <- -sum(X*log(X))
  return(en)
}

entropy.sample.list <- function(X, ce, data, sample.size, sample.time, replace = F){
  En.list <- c()
  for ( i in c(1:sample.time)){
    sample.index <- sample(nrow(data), replace = replace, size = (sample.size*nrow(data)))
    sample.data <- data[sample.index,]
    P <- do.calculus(ce, sample.data)
    L <- list()
    k <- 1
    if (X %in% names(P)){
      for (i in as.character(unique(P[,X]))) {
        L[[k]] <- P[which(P[,X] == i),]
        k <- k + 1
      }
    }
    else{
      L[[1]] <- P
    }
    En <- c()
    for (i in c(1:length(L))) {
      En[i] <- entropy(L[[i]]$new.P)
    }
    
    En.list <- as.data.frame(rbind(En.list, En))
  }
  return(En.list)
}

entropy.cor <- function(en.list){
  if (ncol(en.list) > 1){
    p <- kruskal.test(en.list)
    p.value <- p$p.value
    stat <- as.numeric(p$statistic)
    effective.size <- (stat-ncol(en.list)+1)/(nrow(en.list)*ncol(en.list) - ncol(en.list))
  }
  else{
    p.value <- 1
  }
  cor.size <- length(names(en.list))
  en.cor.matrix <- matrix(0, nrow = cor.size, ncol = cor.size)
  for (i in c(1:cor.size)) {
    for (j in c(1:cor.size) ) {
      en.cor.matrix[i,j] <- sum(abs(en.list[,i] - en.list[,j]))/nrow(en.list)
    }
    
  }
  row.names(en.cor.matrix) <- names(en.list)
  colnames(en.cor.matrix) <- names(en.list)
  if (ncol(en.list) > 1){
    delta.entropy <- sum(en.cor.matrix)/(2*(length(names(en.list))*(length(names(en.list)) - 1)/2))
  }
  else{
    delta.entropy <- 0
    return(list(delta.entropy = delta.entropy, p.value = p.value))
  }
  
  return(list(delta.entropy = delta.entropy, p.value = p.value, effective.size = effective.size))
}

entrop.test <- function(X, ce, data, sample.size, sample.time, replace = F){
  
  en.list <- entropy.sample.list(X, ce, data, sample.size, sample.time, replace)
  cor.p.value <- entropy.cor(en.list)
  return(cor.p.value)
}# X is the do variable

integrate.entropy <- function(y1, y2, delta.y){
  return(0.5*(entropy(y1) + entropy(y2))*delta.y)
}


continue.var.test <- function(Y, X, pl){
  if (!(X %in% names(pl))){
    return(list(delta.entropy = 0, p.value = 1))
  }
  else{
    pl <- pl[,c(Y, X, 'new.P')]  
    L <- list()
    pl <- sapply(pl, as.numeric)
    pl <- as.data.frame(pl)
    k <- 1
    for (i in as.character(unique(pl[,X]))) {
      L[[k]] <- pl[which(pl[,X] == i),]
      k <- k + 1
    }
    I.P <- NULL
    for (j in c(1:length(L))) {
      p.list <- L[[j]]
      Y.name <- Y
      p.list <- p.list[order(p.list[,Y.name]),]
      I <- 0
      if (nrow(p.list) == 1){
        I <- entropy(p.list[,'new.P'])
      }
      else{
        for (k in c(1:(length(p.list[,Y.name]) -1 ))) {
          I <- I + integrate.entropy(p.list[k,'new.P'], p.list[k+1, 'new.P'], as.numeric(p.list[k+1, Y.name])- as.numeric(p.list[k, Y.name]))
        }
      }
      I.P <- rbind(I.P, as.data.frame(cbind(I,p.list[1,X])))
      
    }
    names(I.P) <- c('Integrate','x')
    I.P <- I.P[order(I.P[,'x']),]
    
    delta.entropy <- c()
    for (i in c(1:(length(I.P$Integrate) - 1))) {
      delta.entropy[i] <- (as.numeric(I.P$Integrate[i+1]) - as.numeric(I.P$Integrate[i]))/(as.numeric(I.P$x[i +1])- as.numeric(I.P$x[i]))
    }
    b <- jitter(rep(0, (length(delta.entropy)))) 
    delta.jitter <- (b - rep(0, (length(delta.entropy))))
    a <- delta.entropy + delta.jitter
    
    p.value.mean <- wilcox.test(a, b, exact = T)$p.value
    p.value.var <- ansari.test(a,b, exact = T)$p.value
    s <- abs(delta.entropy) - mean(abs(delta.entropy))
    Z.score.var <- mean(abs(s))/sd(abs(s))
    effect.size.var <- (Z.score.var)/(sqrt(length(s)))
    SE <- sd(abs(delta.entropy))/(sqrt(length(s)))
    delta.entropy <- mean(abs(delta.entropy))
    
    return(list(delta.entropy = delta.entropy, p.value.mean = p.value.mean,  p.value.var = p.value.var, effect.size = effect.size.var, SE = SE))
  }
  
  causal.infer.continue <- function(y, x, data, sample.size.CDF, sample.time.CDF, replace.CDF = F, unobserved = T){
    DAG <- DAG.with.without.unob(data, sample.size.CDF, sample.time.CDF, replace.CDF)
    ce <- ce.func(y,x,DAG, unobserved)
    pl <- do.calculus(ce, data)
    out <- continue.var.test(x, pl)
    return(out)
  }
}