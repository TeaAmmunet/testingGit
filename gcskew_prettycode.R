library(seqinr)#something here
#something more here to test git!
library(ape)
eclDNAbin<-read.GenBank("NC_008563")
attributes(eclDNAbin)
eclFasta<-write.dna(eclDNAbin,file ="eclFasta.fasta", format = "fasta")
eclSeq<-read.fasta("eclFasta.fasta")
nclseq<-sapply(mySeq,function(x) getSequence(x))
outgc<-list()
gc<-function(myseq){
  for (i in 1:length(myseq)){
    if (myseq[i] == "g"){
      outgc[i]<-1
    } else if (myseq[i] =="c"){
      outgc[i]<-(-1)
    } else{
      outgc[i]<-0
    }
  }
  return(outgc)
}
gclist<-unlist(gc(nclseq),recursive = TRUE)
xpos<-c(1:length(gclist))
#Plotting gc3 skew with sliding window
plotgc<-function(mygc, mypos){
  windowlen<-1000
  step<-100
  start<-1
  end<-start+windowlen
  gcskew<-list()
  x<-list()
  y<-list()
  for (i in 1:(length(mygc)/step)){
    if (end <= length(mygc)){
      gcskew[i]<-sum(mygc[start:end])
      x[i]<-mypos[median(c(start:end))]
      start<-start+step-1
      end<-start+windowlen
    } else {
      gcskew[i]<-sum(sum(mygc[start]:mygc[length(mygc)]), sum(mygc[1]:mygc[(end-length(mygc))]))
      x[i]<-mypos[median(c(start:end))]
    }
  }
  #print(gcskew[1:10])
  #y<-cumsum(unlist(gcskew, recursive = TRUE))
  y<-unlist(gcskew,recursive = TRUE)
  y2<-cumsum(y)
  #print(y[length(y)/3:length(y)/3+10])
  #print(y2[length(y2)/3:length(y2)/3+10])
  x<-unlist(x, recursive = TRUE)
  #print(tail(x))
  if (length(y2)!=length(x)){
    print('y and x not of the same length!')
  } else{
    plot(x,y, type='l', xlim=c(0,max(x)), ylim = c(min(y),max(y)))
    par(new=T)
    plot(x, y2, type='l', xlim =c(0, max(x)), ylim = c(min(y2),max(y2)), axes=F)
    axis(side=4)
    ter<-x[match(max(y2),y2)]
    ori<-x[match(min(y2),y2)]
    cat('The origin is at position', ori);cat('  and the terminus at position', ter )
  }
}
plotgc(gclist, xpos)
