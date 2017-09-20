library(seqinr)#something here
#something more here to test git!
library(ape)
#eclDNAbin<-read.GenBank("NC_008563")
myDNAbin<-read.GenBank("CP011051.1")
myDNAbin2<-read.GenBank("AP014658.1")
myDNAbin3<-read.GenBank("NC_012967")
myDNAbin4<-read.GenBank("NC_009719.1")
myDNAbin5<-read.GenBank("NC_016639.1")
#attributes(eclDNAbin)
getfastafiles<-function(accession){
  myDNAbin<-read.GenBank("accession")
  myFasta<-write.dna(myDNAbin,file ="myFasta.fasta", format = "fasta")
  mySeq<-read.fasta("myFasta.fasta")
}
getDNAbin<-function(mybinfile){
  myFasta<-write.dna(mybinfile,file ="myFasta.fasta", format = "fasta")
  mySeq<-read.fasta("myFasta.fasta")
}
#eclFasta<-write.dna(eclDNAbin,file ="eclFasta.fasta", format = "fasta")
#eclSeq<-read.fasta("eclFasta.fasta")
#------------------------GCskew on genomic region-----------------------------------
#system ('prodigal -h')
system('prodigal -i eclFasta.fasta -o eclGenCoord.fasta -d eclGenSeqs.fasta')
system('prodigal -i myFasta.fasta -o myGenCoord.fasta -d myGenSeqs.fasta')
#eclgenfasta2<-read.fasta('eclGenCoord.fasta')
#eclgenseq<-read.fasta('eclGenSeqs.fasta')
mygenseq<-read.fasta('myGenSeqs.fasta')
#eclannot<-sapply(eclgenseq,function(x) strsplit(getAnnot(x),'#'))
myannot<-sapply(mygenseq,function(x) strsplit(getAnnot(x),'#'))
#annotstr<-sapply(eclannot, function(x) strsplit(x,"#"))
annotinfo<- lapply(lapply(myannot,'[', 1:4),function(x) as.integer(x[2:4]))
positcod<-list()
negcod<-list()
seqpos<-list()
codonseqpos<-list()
seqpos_neg<-list()
codonseqpos_neg<-list()
threes<-function(mydata, myannot){
  #This function gets the third codons and their positions for the strand marked (+1)(or -1) in prodigal
  #it returns the strand vectors
  for (i in 1:length(mydata)){
    if (names(mydata[i])==names(myannot[i])){
      mydna<-getSequence(mydata[[i]])
      #getting third nucleotides from sequence on the +1 strand
      if (myannot[[i]][3]==1){
        #print('plus strand')
        #positcod list creation and concatenation works, tried with two elements
        positcod<-c(positcod,mydna[seq(3,length(mydna),3)])
        #Setting the positions of the nucleotides to another list
        seqpos<-c(myannot[[i]][1]:myannot[[i]][2])
        codonseqpos<-c(codonseqpos,seqpos[seq(3,length(seqpos),3)])
      } else {
        #print('minus strand')
        #Getting third nucleotides from sequence on the -1 strand
        negcod<-c(negcod,mydna[seq(3,length(mydna),3)])
        seqpos_neg<-c(myannot[[i]][1]:myannot[[i]][2])
        codonseqpos_neg<-c(codonseqpos_neg,seqpos_neg[seq(3,length(seqpos_neg),3)])
        #print(unlist(codonseqpos_neg, recursive=TRUE)) #works!
      }
    } else {  #If names do not match
      print('Datasets do not match!')
    }
  }
  positcod<-unlist(positcod,recursive=TRUE)
  codonseqpos<-unlist(codonseqpos, recursive=TRUE)
  negcod<-unlist(negcod, recursive = TRUE)
  codonseqpos_neg<-unlist(codonseqpos_neg, recursive=TRUE)
  thirds<-list(plusncl=positcod,pluspos=codonseqpos,minusncl=negcod,minuspos=codonseqpos_neg)
  #thirds<-list(plusncl=positcod,pluspos=codonseqpos)
  return(thirds)
}
thirds<-threes(mygenseq,annotinfo)
gc3<-function(mycodons){
  for (i in 1:length(mycodons)){
      if ((mycodons[i] == "g")||(mycodons[i]=='G')){
        outgc3[i]<-1
      } else if ((mycodons[i] =="c") ||(mycodons[i]=='C')){
        outgc3[i]<-(-1)
      } else{
        outgc3[i]<-0
      }
    if (i %%1000==0){
      cat('processing nucleotide', i, '\n')
    }
    }
  return(outgc3)
}
gc3list<-gc3(thirds$plusncl)
gc3list<-unlist(gc3list, recursive = TRUE)
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
  #print(length(x)/2)
  if (length(y2)!=length(x)){
    print('y and x not of the same length!')
  } else{
    plot(x,y, type='l', xlim=c(0,max(x)), ylim = c(min(y),max(y)))
    par(new=T)
    plot(x, y2, type='l', xlim =c(0, max(x)), ylim = c(min(y2),max(y2)), axes=F)
    axis(side=4)
    ter<-x[match(max(y2),y2)]
    ori<-x[match(min(y2),y2)]
    cat('The origin is at position', ori);cat(' and the terminus at position', ter )
  }
}
plotgc(gc3list, thirds$pluspos)
