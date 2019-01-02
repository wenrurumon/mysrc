imrate <- function(pidi,ifscale=T){
  #Specific for FL data
  #INPUTE data with constant value within specific range
  i <- filter(rate,pid==pidi) %>% arrange(rate)
  pos <- filter(base,pid==pidi)$hours
  pos <- ceiling(pos)
  if(ifscale){
    i <- i %>% mutate(hour1=hour1/pos*100+1,hour2=hour2/pos*100+1)
    pos <- 100
  }
  out <- rep(0,length=pos+1)
  if(nrow(i)==0){
    return(out)
  }
  for(j in 1:nrow(i)){
    out[(floor(i$hour1[j])*(i$hour1[j]>0)):min(ceiling(i$hour2[j]),(pos+1))] <- i$rate[j]
  }
  out <- out[1:(pos+1)]
  return(out)
}
imrange <- function(i,pos,ifscale=T){
  #Specific for FL data
  #INPUTE data with constant value within specific range
  pos <- ceiling(pos)
  if(ifscale){
    i <- i %>% mutate(hour1=hour1/pos*100+1,hour2=hour2/pos*100+1)
    pos <- 100
  }
  out <- rep(0,length=pos+1)
  if(nrow(i)==0){
    return(out)
  }
  for(j in 1:nrow(i)){
    out[(floor(i$hour1[j])*(i$hour1[j]>0)):min(ceiling(i$hour2[j]),(pos+1))] <- 1
  }
  out <- out[1:(pos+1)]
  return(out)
}
