i=commandArgs(trailingOnly = TRUE)
load(file="vars0")
source("manytree.R")
source("treeCrawler.R")

trees=manytree(vGmix, Gmix)
save(trees, file=paste("trees100","part",i,sep="_",collapse=""))
