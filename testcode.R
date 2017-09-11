library(seqinr)
install.packages("micropan")
source("http://bioconductor.org/biocLite.R")
biocLite()
library(BiocInstaller)
biocVersion()
biocValid()
biocLite(c("GenomicRanges","biomaRt", "AnnotationDbi" ))
biocLite("Biostrings")
biocLite("BSgenome.Ecoli.NCBI.20080805")
library(GenomicRanges)
library(biomaRt)
library(BSgenome.Ecoli.NCBI.20080805)
library(Biostrings)
#help("useMart")
#mart=useMart("ensembl")
#listDatasets(mart)
#Importing ecoli genome to test
ecl<-Ecoli$NC_008563
#counting nucleotides (this one counts AAs as well)
alphabetFrequency(ecl)
#counting only nucleotides
letterFrequency(ecl, letters = "ACGT", OR=0)
#trying to save as DNA only and counting w alphabets
eclDna<-DNAString(ecl)
eclDna
alphabetFrequency(eclDna) #No, still counts all
alphabetFrequency(eclDna, baseOnly=TRUE) #This works! 
window<-1000
gc<-rowSums(letterFrequencyInSlidingView(ecl, window, c("G", "C")))/window
plot(gc, type='l')
plot(1:length(gc),gc)
lines(lowess(x = 1:length(gc), y= gc, f = 0.10), col = 12, lwd = 2)
#Cumulative GC skew curve
##First to create the vector y and x (x = length of sequence)

sq<-as.character(eclDna)
#eclDNA[seq]
#ecl
#sq[[1]]
splitsq<-strsplit(sq,"")[[1]]
splitsq[2]=="A"
y<- vector(,length=length(splitsq))
#trying first on shorter seq
for (i in 1:length(splitsq)){
  if (splitsq[i] == "G"){
  y[i]=1
  } else if (splitsq[i] =="C"){
    y[i]=-1
  } else
    y[i]=0
}
#length(y)
cumy<-vector()
#cumy[1]+y2[2]
#y2[2][1]cumy[2]<-cumy[1]+y2[2]
#instead of whole genome, using overlapping windows
#wind<-window(y,width=1000) #does not work, meant for time series?
#lengths(wind)
wind<-1000
step<-100
windylen<- ceiling((length(y)-window)/step)
windy<-vector(,length=windylen)
start<-1
for (j in 1:windylen){
  windy[j]<-sum(y[start:(start+wind)])
    start<-start+step
}
#windy
cumy<-cumsum(windy)
#x<-vector(,length=length(splitsq))
#x<-rep(1:length(splitsq))
#plot(x,cumy)
windx<-rep(1:windylen)
plot(windx,cumy, type='l')
#lines(windx,cumy)
#length(x)
#length(cumy)
#length(splitsq)
#Testing getting annotations through seqinr
choosebank("genbank")
#getting fasta file for E.coli NC_008563
library(ape)
eclDNAbin<-read.GenBank("NC_008563")
attributes(eclDNAbin)
eclFasta<-write.dna(eclDNAbin,file ="eclFasta.fasta", format = "fasta")
eclSeq<-read.fasta("eclFasta.fasta")
eclSeq
eclAnnot<-getAnnot(eclSeq[[1]])
eclAnnot
?getAnnot
methods(getAnnot)
eclDNAseq<-readBStringSet("eclFasta.fasta", format='fasta')
eclDNAseq<-gsub(" ", "", eclDNAseq, fixed = TRUE)
eclDNAseq['seq']
library(micropan)
extdata.path <- file.path(path.package("micropan"),"extdata")
prodigalPredict("eclFasta.fasta","eclPredGen.fasta")
