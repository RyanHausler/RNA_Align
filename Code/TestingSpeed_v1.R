library(DECIPHER)
library(gplots)

memory.limit(size=50000)
rna <- readRNAStringSet("C:/Users/Ryan/Documents/Wright/Seeds/SILVA_128_LSURef_tax_silva_full_align_trunc.fasta.gz") # n=1000

mask = MaskAlignment(rna, threshold=1)
rna = RNAStringSet(mask)

coreTime=numeric()
# Time stamps for core usage
for(i in seq(1:detectCores())){
  if(i==1){
    coreTime[1] = system.time(PredictDBN(rna, processors = i))
  }
  else{
    coreTime = c(coreTime, system.time(PredictDBN(rna, processors = i)))
  }
  
}

n = numeric(length(seq(2000, length(rna), 2000)))
seq1 = strsplit(as.character(rna[1]), "")
seq1 = matrix(unlist(seq1))
l = numeric(length(seq(50, dim(seq1)[1], 50)))
t = matrix(data=0, nrow=length(n), ncol=length(l))

x=1
y=1
for(i in seq(2000, length(rna), 2000)){
  n[x] = i
  x=x+1
  for(j in seq(50, dim(seq1)[1], 50)){
    l[y] = j
    rna_test = subseq(rna, 1, j)
    t[x, y] = system.time(PredictDBN(rna_test[1:i]))[3]
    y=y+1
  }
  y=1
}

rownames(t) = n
colnames(t) = l

heatmap.2(t, Rowv=NA, Colv=NA, col=heat.colors(256), main='Performance of PredictDBN', xlab='Number of Seqs', ylab='Length of Seqs')