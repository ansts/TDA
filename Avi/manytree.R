manytree <- function(vGmix, Gmix, n=16) {
  require(parallel)
  
  cl=makeCluster(2)  # Apply n iterations of the treeCrawler iterative 
  
  clusterExport(cl, list("vGmix","Gmix","treeCrawler","rocrack"))
  trees=parSapply(cl, 1:n,function(i){
    trcr=treeCrawler(vGmix, Gmix, Iterations=3, s=10000)
    return(trcr)
  })
  stopCluster(cl)
  return(trees)
}