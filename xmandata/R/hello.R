# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}

library(rjson)
rj0 <- function(x){
  x <- fromJSON(paste(readLines(x),collapse=' '))
  return(x)
}
rj1 <- function(x){
  x <- fromJSON(paste(readLines(x),collapse=' '))
  xtitle <- sapply(x[[1]],paste,collapse='::')
  x <- do.call(rbind,x[[2]])
  colnames(x) <- xtitle
  x
}
rj <- function(x){try(rj1(x))}
