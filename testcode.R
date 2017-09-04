library(seqinr)
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
