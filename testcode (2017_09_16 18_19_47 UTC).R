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
windylen<- ceiling((length(y)-window)/step)#think of changing so that all pieces the same size
#, and last piece overlaps with the beginning
#also check x axis
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
system("prodigal -h")
###------------------Prodigal and GC3 things begin here-------------------
#From the results, we need to 1. read.fasta 
#2. chage it into string and split it 
#3. get the start + stop 
#4. change start and stop into float
#eclGen<-read.fasta(file="eclPredGen.fasta")
#Calling prodigal directly without micropan
system ('prodigal -h')
system('prodigal -i eclFasta.fasta -o eclGenCoord.fasta -d eclGenSeqs.fasta')
#genfile<-system.file('eclPredgen.gbk', package="genbankr")
#Getting genbankr package for parsing gbk files
#biocLite("genbankr")
#library(genbankr)
#gb<-readGenBank(genfile)
#testgenb=readGenBank(system.file('sample.gbk',package='genbankr'))#Doesnt work
#Trying parsing with regular delim import
#delims<-c('\t',';','.')
#eclGFF<-read.delim('eclGenCoord.gff', header=TRUE,sep='\t', skip=2)
#head(eclGFF)
#eclGFF[[1]][3]
#firstg<-eclGFF[[1]][1]
#charg<-as.character(firstg)
#charg
#-> hard to extract gene start and stop
#Trying with GenFeatures instead
#library(GenomicFeatures)
#eclgenabank<-readGenBank(system.file("eclPredgen.gbk", package='genbankr'))
#ecltxdb<-makeTxDbFromGenBank('eclPredgen.gbk')#complicated does not work
#--------THIS IS THE CURRENT WAY-------------------
#trying with fasta again
#eclgenfasta<-read.delim('eclPredGen.fasta', sep='', fill=TRUE)
#head(eclgenfasta)
#eclgenfasta #not easy from here either...weird Xs before numbers
eclgenfasta2<-read.fasta('eclPredGen.fasta')
eclgenseq<-read.fasta('eclGenSeqs.fasta')
head(eclgenfasta2)
head(eclgenseq)
getAnnot(eclgenfasta2)
typeof(eclgenseq)
#getSequence.SeqFastadna(eclgenfasta2)
eclannot<-sapply(eclgenfasta2,function(x) getAnnot(x))
#Try9ing to get 3rd codon positions with one gene first
oneseq<-getSequence(eclgenseq[[1]])
thirdstry<-oneseq[seq(1,length(oneseq),3)]
thirdstry
oneseq[1]
typeof(eclgenseq)
#in principle, for every element (gene) in eclgenseq, take every third element to a new vector
#outhird<-vector(mode='character')
threes<-function(mydata){
  dna<-getSequence(mydata)
  #print('getSeq works')
  everythird<-dna[seq(1,length(dna),3)]
  #append(outhird,everythird, after=length(outhird))
  #outhird
  as.vector(everythird)
}
eclgc3s<-lapply(eclgenseq,FUN=threes)
#headgc3<-head(lapply(eclgenseq,FUN=threes))
head(eclgc3s)
#headgc3
#Trying to get rid of names, does not work
#library(rlang)
#flatgc3<-flatten(eclgc3s)
#flatter<-flatten(flatgc3)
#remove(flatter)
#remove(flatgc3)
#head(flatter)
#length(eclgc3s[1:2])
#input= list with headings
gc3<-function(mycodons){
  outgc3<-vector('list')
  for (j in 1:length(mycodons)){
    for (i in 1:length(mycodons[[j]])){
      if (mycodons[[j]][i] == "g"){
        c(outgc3,1)
      } else if (mycodons[[j]][i] =="c"){
        c(outgc3,(-1))
      } else{
        c(outgc3,0)
      }
    }
  }
  outgc3
}
gc3(eclgc3s[1:2])
head(gc3y)
gc3(eclgc3s[1:2])
#------------old prodigal fasta manipulations from here down------------------
eclannot[[1]]
typeof(eclannot)
eclsplit<-sapply(eclannot, function(x) strsplit(x,"#"))
eclsplit[[1]][2:3]
head(eclsplit)
eclcoord<-lapply(lapply(eclsplit, '[', 1:4),function(x) as.integer(x[2:4]))#taking the start and the stop and converting them into integers
head(eclcoord)
typeof(eclcoord)
eclgenes<-data.frame(t(sapply(eclcoord,c)))
eclgenes[[1]][2]
colnames(eclgenes)<-c('start','stop', 'direction')
head(eclgenes)
eclgenes$len=(eclgenes$stop-eclgenes$start)
#Getting codons from one strand
#calculating codons
#codons within a gene
c<-vector(mode='numeric')
dim(eclgenes)[1]
tail(eclgenes)
length(eclgenes[['start']])
#testfu<-function(mydata){
#  for (i in 1:dim(mydata)[1]){
#    c[i]=i
#  }
#return(c)
#}
#testfu(eclgenes)
(eclgenes$len[1]+1)/3
eclcods<-vector(mode='numeric')
a<-c(1,2,3)
eclcods[1]<-0
append(eclcods,a,after=length(eclcods))
codstart[1]<-eclgenes[['start']][1]
codfun<-function(mydata){
  if(mydata[['direction']]==1){
    print('+1 strand')
    nrcodons<-(mydata[['len']]+1)/3
    for (i in 2:nrcodons){
      codstart[i]<-codstart[i-1]+3
    }
    print('for loop ok')
  }
  append(eclcods,codstart, after=length(eclcods))
  #return(eclcods)
  codstart
  eclcods
}
codfun(eclgenes[1,])
eclgencod<-apply(eclgenes,1,FUN=codfun) #General looping works, but output is empty, and does not return codstart
head(codstart)
