treeCrawler=function(L, G, s=100, Iterations=5, ini=NULL, mxc=100, nn="C"){
  require(igraph)
  require(stringi)
  
  g=induced.subgraph(G,L)
  fl1=T     
  cl=lapply(1:Iterations, function(i){
    proct=proc.time()
    print(c("Iteration ",i))
    clu=cluster_label_prop(g, initial = ini)
    print(proc.time()-proct)
      if (max(clu$membership)==1) return (list(L)) else return(communities(clu))
    })
  l=lengths(cl)
  if (max(l)==1) {
    print(c(length(L),"single"))
    mcl=1
    iniv=ini    
    x=rocrack(mcl, mxc, g, iniv)
    mcl=x[[1]]
    cl=x[[2]]
    if(mcl==1) {
        print("Failed split")
        fl1=F           
      }
  }
  else {
    cl=cl[[which.max(l)]] 
  }
  ni=paste(nn,seq_along(cl), sep="_")

  x=lapply(seq_along(cl),function(i){
    co=cl[i]
    if (is.null(ini)) ini1=NULL else {
      x=ini[co[[1]]]
      ini1=factor(x)
      ini1=as.double(ini1)
      names(ini1)=names(x)
    }
    if (length(co[[1]])>=s & fl1) {                  
        nii=ni[i]
        y=treeCrawler(co[[1]],g, Iterations=Iterations, ini=ini1, s=s, nn=nii)
        return(y)
      }
      else {
        names(co)=ni[i]
        # if (length(co[[1]])==1) {
        #   co=unlist(co,recursive = F)
        # }
        # else{
          co=lapply(co[[1]],function(y) list(Seq=y))
          names(co)=paste(ni[i],seq_along(co), sep="_")
        # }
          print(names(co))
        return(co)
      }
  })

  names(x)=ni
  return(x)
}
rocrack <- function(mcl, mxc, g, iniv) {
  i=1
  while (mcl==1 & i<mxc){
    proct=proc.time()
    clu=cluster_label_prop(g, initial = iniv)
    cl=communities(clu)
    mcl=max(clu$membership)
    print(c("Rocrac iteration ",i,mcl))
    print(proc.time()-proct)
    i=i+1
  }
  return(list(mcl, cl))
}