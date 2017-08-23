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
dev.new()
png("~/Desktop/RNA/Cores_vs_Time.png")
plot(coreTime,
  xlab="Number of cores",
  ylab="Elapsed time (secs)")
lines(25269*seq_len(nCores)^-.6289586+-1460.803) # best fitting line
dev.off()

# find best fitting line
plot <- FALSE # add lines to plot?
opt <- function(p) {
  if (plot)
    lines(p[3]*(seq_len(nCores)^p[2]) + p[1], col="red")
  # return sum of squared error
  return(sum((coreTime - (p[3]*(seq_len(nCores)^p[2]) + p[1]))^2))
}
o <- optim(c(2000, -0.5, 10000), opt)
o$par[2] # t~ncores^result

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
heatmap(t, Rowv=NA, Colv=NA, col=heat.colors(256), main='Performance of PredictDBN', xlab='Number of Seqs', ylab='Length of Seqs')

dev.new()
# png("~/Desktop/RNA/Number_Seqs_vs_Time.png")
matplot(n,
  t,
  xlab="Number of sequences",
  ylab="Elapsed time (secs)",
  type="l",
  log="xy")
# dev.off()

dev.new()
png("~/Desktop/RNA/Length_Seqs_vs_Time.png")
matplot(l,
  t(t),
  xlab="Length of sequences",
  ylab="Elapsed time (secs)",
  type="l",
  log="xy")
dev.off()

#fit functions
apply(log(t), 1, function(y) coef(lm(y~log(l)))[2]) # t~O(l^result)
apply(log(t), 2, function(y) coef(lm(y~log(n)))[2]) # t~O(n^result)

