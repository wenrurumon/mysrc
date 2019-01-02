jm <- function(pidi){
  # print(k<<-k+1)
  trate <- imrate(pidi)
  i.base <- filter(base,pid==pidi)
  i.idx <- filter(idx,pid==pidi)
  i.idx <- sapply(tag.idx,function(i){
    i <- filter(i.idx,item==i) %>%
      group_by(pid,hours,item) %>% summarise(value=mean(value))
    imspline(i$hours,i$value,i.base$hours,0.1,F,T)
  })
  i.exe <- filter(fexe,pid==pidi)
  i.exe <- sapply(tag.exe,function(i){
    i <- filter(i.exe,item==i)
    imrange(i,i.base$hours)
  })
  rlt <- data.table(pid=pidi,com=i.base$com_or_not,diag_no=i.base$diag_no,
                    orate=i.base$orate,hours=i.base$hours,
                    i.idx,i.exe,trate=trate)
  rlt
}
