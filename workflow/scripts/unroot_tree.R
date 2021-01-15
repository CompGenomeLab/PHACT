#!/usr/bin/env Rscript
library(ape)
args = commandArgs(trailingOnly=TRUE)
unroot_tree <- function(tree_file) {
  vert.tree<-read.tree(tree_file)
  unro_tree<- unroot(vert.tree)
  unroot_tree_name<- deparse(substitute(tree_file))
  write.tree(unro_tree,file=paste(tree_file,"_unrooted",sep=""))}
  
unroot_tree(tree_file=args[1])
