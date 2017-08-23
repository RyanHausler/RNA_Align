library(DECIPHER)

#memory.limit(size=50000)
rna <- readRNAStringSet('/home/eswright/Downloads/SILVA_128_LSURef_tax_silva_full_align_trunc.fasta.gz')

mask <- MaskAlignment(rna,
  type="ranges",
  threshold=0,
  maxFractionGaps=0.5,
  includeTerminalGaps=TRUE,
  showPlot=T)
rna <- replaceAt(rna, mask)

# Time stamps for core usage
nCores <- detectCores()/2 # number of cores = threads/2
coreTime <- numeric(nCores)
for(i in seq_len(nCores)){
  print(i)
  coreTime[i] <- system.time(PredictDBN(rna, processors=i))[3]
}
plot(coreTime,
  xlab="Number of cores",
  ylab="Elapsed time (secs)",
  log="y")

n <- seq(100, length(rna), 20000) # number of sequences
l <- seq(200, width(rna[1]), 200) # length of sequences
t <- matrix(data=0, nrow=length(n), ncol=length(l), dimnames=list(n, l))

for (i in seq_along(n)) {
  for (j in seq_along(l)) {
    cat("\ni =", i, "j =", j)
    rna_test <- subseq(rna[sample(length(rna), n[i])], 1, l[j])
    t[i, j] <- system.time(PredictDBN(rna_test, processors=nCores, v=F))[3]
  }
}

save(nCores, coreTime, n, l, t, file="~/Desktop/TestingSpeed_v1.RData")

dev.new()
pdf("~/Desktop/RNA/Heatmap.pdf")
heatmap(t, Rowv=NA, Colv=NA, col=heat.colors(256), main='Performance of PredictDBN', xlab='Number of Seqs', ylab='Length of Seqs')
dev.off()

dev.new()
pdf("~/Desktop/RNA/Number_Seqs_vs_Time.pdf")
matplot(n,
  t,
  xlab="Number of sequences",
  ylab="Elapsed time (secs)",
  type="l",
  log="y")
dev.off()

dev.new()
pdf("~/Desktop/RNA/Length_Seqs_vs_Time.pdf")
matplot(l,
  t(t),
  xlab="Length of sequences",
  ylab="Elapsed time (secs)",
  type="l",
  log="y")
dev.off()
