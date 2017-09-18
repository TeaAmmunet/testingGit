library(seqinr)
library(ape)
eclDNAbin<-read.GenBank("NC_008563")
attributes(eclDNAbin)
eclFasta<-write.dna(eclDNAbin,file ="eclFasta.fasta", format = "fasta")
eclSeq<-read.fasta("eclFasta.fasta")
system ('prodigal -h')
system('prodigal -i eclFasta.fasta -o eclGenCoord.fasta -d eclGenSeqs.fasta')
#eclgenfasta2<-read.fasta('eclGenCoord.fasta')
eclgenseq<-read.fasta('eclGenSeqs.fasta')
eclannot<-sapply(eclgenseq,function(x) strsplit(getAnnot(x),'#'))
#annotstr<-sapply(eclannot, function(x) strsplit(x,"#"))
annotinfo<- lapply(lapply(eclannot,'[', 1:4),function(x) as.integer(x[2:4]))
threes<-function(mydata){
  #This function gets the third codons and their positions for the strand marked (+1) in prodigal
  dna<-getSequence(mydata)
  #print('getSeq works')
  everythird<-dna[seq(1,length(dna),3)]
  #append(outhird,everythird, after=length(outhird))
  #outhird
  as.vector(everythird)
}
eclgc3s<-lapply(eclgenseq,FUN=threes)
outgc3<-list()
gc3<-function(mycodons){
  for (j in 1:length(mycodons)){
    for (i in 1:length(mycodons[[j]])){
      if (mycodons[[j]][i] == "g"){
        outgc3<-c(outgc3,1)
      } else if (mycodons[[j]][i] =="c"){
        outgc3<-c(outgc3,(-1))
      } else{
        outgc3<-c(outgc3,0)
      }
    }
  }
  return(outgc3)
}
mocCG<-list(c("g", "c", "a", "t", "g"), c("a","t","g","c","c"))
shortgc3<-function(mydata){
  #outgc3<-list()
  for (i in 1:length(mydata[[1]])){
    if (mydata[[1]][i] == "g"){
       outgc3<-c(outgc3,1)}else {outgc3<-c(outgc3,0)}
    #return(outgc3)
  }
  return(outgc3)
}
gclist12<-gc3(eclgc3s[1:2])
plotgc<-function(mygc){
  windowlen<-1000
  step<-100
  start<-1
  end<-start+windowlen
  for (i in ceiling(length(mygc)/step)){
    if (start==1){
      gcskew<-c(gcskew,sum(mygc[start]:mygc[end]))
      x[i]<-(end-start)+1/2
      start<-end+1
    } else {
      if (end > length(mygc)){
        gcskew<-c(gcskew,sum(sum(mygc[start]:mygc[length(mygc)]), sum(mygc[1]:mygc[(end-length(mygc))])))
        x[i]<-(end-start)+1/2
      } else {
        gcskew<-c(gcskew,sum(mygc[start]):mygc[end])
        x[i]<-(end-start)+1/2
      }
    }
  }
  y<-cumsum(gcskew)
  x<-
}
