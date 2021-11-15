# Aligns a 7-mer to cooc model

aligncc<- function(p, cc) {
  p=unlist(strsplit(p, split=""))
  x=unlist(lapply(1:6, function(i){
            lapply((i+1):7, function(j){
              paste(p[c(i,j)], sep="", collapse="")
            })
  }))
  i=c(1,2,3,4,5,6,1,2,3,4,5,1,2,3,4,1,2,3,1,2,1)
  xy=aggregate(1:21, by=data.frame(x,i), "length")
  xy[,2]=as.numeric(xy[,2])
  xy[,3]=as.numeric(xy[,3])
  res=sum(sapply(seq_along(xy$x),function(i) cc[rownames(cc)==xy[i,1],xy[i,2]]*xy[i,3]))
  return(res)
}