require(igraph)
require(parallel)
require(reshape2)
require(data.table)
require(corrplot)
require(pbapply)
require(MASS)
require(umap)
require(rgl)
require(PST)
require(Biostrings)
require(seqinr)
require(ggseqlogo)
require(matrixStats)
require(stringi)

load("Gmix") # Graph of 553634 7-aa sequences linked at longest common subsequence > 4 aa
load("fl") #labels of the sequences 1 - Control, 2 - Patients, 3 - common for 1&2, 4 - Public IgM repertoire, 5 - 1&4, 6 - 2&4, 7 - 3&4

vGmix=V(Gmix)

save(Gmix, vGmix, file-"vars0")
rm(Gmix,fl,vGmix)

#  To Avitohol - scripts in folder Avi ------------------------------------



# From Avitohol -----------------------------------------------------------
N=380000  # Threshold number of singlet clusters
seqs=lapply(1:7,function(i){ 
  print(i)
  fnm=paste("Avi\\trees100_part_",i,sep="") # Load a group of 16 trees, each in the form of a nested list
  load(fnm)
  print(fnm)
  trees=sapply(trees, function(x0){  # Cut each tree at a suitable level
    maxc=c()
    x=x0
    repeat{ # Determine number of clusters ( size > 3 sequences) at each cut level
      x=unlist(x, recursive = F)
      xi=lapply(x, unlist)
      if (length(xi[lengths(xi)==1])>N) break
      n=length(xi[lengths(xi)>3])
      maxc=c(maxc,n)
    }
    x=x0
    maxc=which.max(maxc) # Determine the cut level producing the greatest number 
    print(maxc)
    for (i in 1:maxc) x=unlist(x, recursive = F) # Cut at max level
    x=lapply(x, unlist)
    x=x[lengths(x)>3] # Take the clusters of >3
    print(length(x))
    names(x)=seq_along(x)
    tr=melt(x)
    x=tr$L1
    names(x)=tr$value
    return(x) # Return a list of the sequences with cluster membership
  })
})
seqs=unlist(seqs, recursive = F)
save(seqs, file="Seqs")

x=lapply(seqs,as.matrix) # reshape data
seqsm=melt(x)
rm(seqs,x)
seqsm=seqsm[,-2]
colnames(seqsm)=c("Seq","Cl", "Tr")

l0=rep(0,112)

sdt=data.table(seqsm, key = "Seq")
rm(seqsm)
sdta=sdt[,.(list(c(as.matrix(rbind(as.numeric(Cl),as.numeric(Tr)))))), by="Seq"]
rm(sdt)

proct=proc.time()
ClTr=t(sapply(sdta$V1, function(l){
  n=length(l)
  Tr=(1:(n/2))*2
  Cl=Tr-1
  l0[l[Tr]]=l[Cl]
  return(l0)
}))
print(proc.time()-proct)

# The most stringent scenario when a consensus between all trees is required
ClTra=aggregate(as.character(sdta$Seq), by=as.data.frame(ClTr), "list")
rm(ClTr)
tClTra=table(lengths(ClTra$x))
i=as.numeric(names(tClTra))>100
sum(tClTra[i]) # 528/213178 clusters are larger than 100 sequences
xi=as.numeric(names(tClTra))
sum(xi[xi>100]*tClTra[xi>100])/sum(lengths(ClTra$x))  # and they contain ~20% of the sequences
x=as.matrix(ClTra[lengths(ClTra$x)>100,1:112])
L=nrow(x)
m0=x!=0

j=cut(seq_along(x[,1]), L%/%50, labels = F)
cl=makeCluster(4)
clusterExport(cl, list("x","L","m0","j"), envir = environment())
proct=proc.time()
ALseqscl=pbsapply(1:10,function(ii){
  re=sapply(which(j==ii), function(i){
    l1=x[i,]
    l1all=t(array(rep(l1,L), dim=c(112,L)))
    res=which(rowSums((x==l1all) & m0)>0.67*rowSums(m0)) 
    return(res)
  })
  names(re)=which(j==ii)
  return(re)
}, cl=cl)
stopCluster(cl)
print(proc.time()-proct)

ALseqscl=unlist(ALseqscl, recursive = F)
table(lengths(ALseqscl))
nms=rownames(x)
names(ALseqscl)=nms[as.numeric(names(ALseqscl))]
ALseqscl=lapply(ALseqscl,names)
ELseqscl=melt(ALseqscl)
Gseqscl=graph_from_edgelist(as.matrix(ELseqscl), directed = F)
Gseqscl=simplify(Gseqscl)
plot(Gseqscl, vertex.size=0, vetex.label=NULL)
Gseqsclcmp=components(Gseqscl)

x=ClTra$x[lengths(ClTra$x)>100]
names(x)=nms[seq_along(x)]
mmbrshp=Gseqsclcmp$membership
names(mmbrshp)=as.character(names(mmbrshp))
mmbrshp=mmbrshp[names(x)]
ultimateCl=aggregate(as.array(x), by=list(mmbrshp), "c")
ultimateCl=lapply(ultimateCl$x, unlist)
proct=proc.time()
ultAdjMx=sapply(ultimateCl, function(x1){
            sapply(ultimateCl, function(x2){
              x=c(x1,x2)
              fi=c(rep(1, length(x1)), rep(2, length(x2)))
              names(fi)=x
              g=induced.subgraph(Gmix, x)
              modularity(g, fi[names(V(g))])
            })
})
print(proc.time()-proct)
hist(ultAdjMx[ultAdjMx>0], breaks=100)

K=c(1:50)
stress=t(sapply(K, function(k){
  x=isoMDS(as.dist(ultAdjMx), k=k)
  return(c(k,x$stress))
}))
ultAdjemb=isoMDS(as.dist(ultAdjMx), k=15)
ultAdjemb=ultAdjemb$points

x=stda(ultAdjemb)

uini=umap.defaults
uini$verbose=T
uini$n_neighbors=10
uini$n_components=2
uultemb=umap(ultAdjemb, config=uini)

plot(uultemb$layout)
text(uultemb$layout, texts=seq_along((uultemb$layout[,1])))
cy=c(22,24,32,31,25,35,26,40,22)
lines(uultemb$layout[cy,])

t0=1:7
names(t0)=t0
ultfl=t(sapply(ultimateCl, function(l){
  t=t0
  t1=table(fl[l])
  t[names(t1)]=t1
  return(t)
}))

as=sum(ultfl)
rss=rowSums(ultfl)/as
css=colSums(ultfl)/as
chsqult=apply(ultfl,1,function(l){
  chisq.test(l,p=css, simulate.p.value = T)
})
pchsqult=sapply(chsqult, function(ch) ch$p.value)
adj.pchsqult=p.adjust(pchsqult)
resqult=t(sapply(chsqult, function(ch) ch$stdres))
jo=rowSums(resqult[,c(1,3,5)])-rowSums(resqult[,c(2,3,6)])
rownames(resqult)=seq_along(chsqult)
resqult=resqult[order(jo),]
ultflrcol1=apply(resqult,1,function(l) rgb(sum(10^(l[c(2,3,6)]))/sum(10^l), sum(10^(l[c(1,3,5)]))/sum(10^l), sum(10^l[4])/sum(10^l),1))
plot(uultemb$layout, col=ultflrcol1, pch=16, cex=2, xlab="UMAP D1", ylab="UMAP D2")
lines(uultemb$layout[cy,], lwd=2)
prsclt=pnorm(abs(resqult),lower.tail = F)
prsclt=array(p.adjust(prsclt), dim=dim(resqult))
resqultsig=resqult[(prsclt[,1]<0.05 | prsclt[,2]<0.05)&!(prsclt[,1]<0.05 & prsclt[,2]<0.05),]
posqlt=prsclt[(prsclt[,1]<0.05 | prsclt[,2]<0.05)&!(prsclt[,1]<0.05 & prsclt[,2]<0.05),1:2]
resqultsig[order(posqlt[,2]),]

clAPL=c(10,3,36,9,34)
clCntr=c(18,42,28,17,23)
clsall=c(clAPL,clCntr)
pultsig=unlist(ultimateCl[clsall])


# Co-occurrence ------------------------------------------------------------

coocs=lapply(seq_along(ultimateCl),function(i){
  print(i)
  l=ultimateCl[[i]]
  cooc(l)
})
valigncc=Vectorize(aligncc, vectorize.args="p")

cl=makeCluster(4)
clusterExport(cl,list("ultimateCl","coocs", "valigncc"), envir = environment())
ccomp10=pblapply(seq_along(ultimateCl), function(i){
  l=ultimateCl[[i]]
  v=valigncc(l, coocs[[10]])
  return(v)
}, cl=cl)
stopCluster(cl)

ccomp42=lapply(seq_along(ultimateCl), function(i){
  proct=proc.time()
  l=ultimateCl[[i]]
  print(c(i, length(l)))
  v=valigncc(l, coocs[[42]])
  print(proc.time()-proct)
  return(v)
})

boxplot(ccomp42, outline=F,notch=T, main="Cooccurrence model comparison of cluster 42 to the rest", ylab="Cooc score of individual peptides",xlab="Cluster #" )
boxplot(ccomp10, outline=F,notch=T, main="Cooccurrence model comparison of cluster 10 to the rest", ylab="Cooc score of individual peptides",xlab="Cluster #" )


