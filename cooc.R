
#  A function creating a cooccurrence matrix of the AA at different gap lengths
#  as a statistical model of a cluster. Takes as input a list of 7-aa sequences.  



cooc=function(L){
  require(Biostrings)
  aa=AA_ALPHABET[1:20]

  AA=strsplit(c(sapply(aa,function(a1){
      sapply(aa, function(a2){
       paste(a1,a2, sep="") 
      })
  })), split="")
  AA=as.matrix(t(as.data.frame(AA)))
  rownames(AA)=NULL
  AA=data.frame(AA,array(rep(0,2400), dim=c(400,6)))
  colnames(AA)=c("A1","A2","D1","D2","D3","D4","D5","D6")
  for (p in L){
    n=unlist(strsplit(p, split=""))
    ij=0
    for (i in 1:6){
        for (j in (i+1):7){
          AA[AA$A1==n[i]&AA$A2==n[j],(2+j-i)]=AA[AA$A1==n[i]&AA$A2==n[j],(2+j-i)]+1
        }
    }
  }
  AA[,3:8]=AA[,3:8]/sum(AA[,3:8])
  AA[,1]=paste(AA$A1,AA$A2, sep="")
  AA=AA[,-2]
  rownames(AA)=AA[,1]
  AA=AA[,-1]
  return(AA)
}