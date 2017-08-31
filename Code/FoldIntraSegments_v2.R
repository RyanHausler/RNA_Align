library(DECIPHER)

dna <- readDNAStringSet("~/Desktop/full_PB2.fasta")
cds <- readDNAStringSet('~/Desktop/cds_PB2.fasta')

x = seq(1, length(dna), 10)
clusters <- IdClusters(myXStringSet=dna[x],
	method="inexact",
	cutoff=c(0.01))
a <- apply(clusters, 2, function(x) length(unique(x))) # choose a cutoff
reps <- which(!duplicated(clusters[, 1]))
dna <- dna[reps]

align <- AlignTranslation(dna)

d <- DistanceMatrix(align, verbose=FALSE, includeTerminalGaps=TRUE)
cutoffs <- sort(unique(c(seq(0, 1, 0.01),
                         round(quantile(d, seq(0, 1, 0.01), na.rm=TRUE), 5))))
dimnames(d) <- NULL
# Clusters similar seqs
suppressWarnings(guideTree <- IdClusters(d,
                                         method="UPGMA",
                                         cutoff=cutoffs,
                                         verbose=FALSE))
guideTree <- guideTree
lastSplit <- numeric(dim(guideTree)[1])
weights <- numeric(dim(guideTree)[1])
numGroups <- 1

for (i in (dim(guideTree)[2] - 1):1) {
  if (length(unique(guideTree[, i])) > numGroups) {
    if (numGroups==1) {
      lastSplit[] <- cutoffs[i]
    } else {
      for (j in unique(guideTree[, i + 1])) {
        w <- which(guideTree[, i + 1]==j)
        if (length(unique(guideTree[w, i])) > 1) {
          weights[w] <- weights[w] + (lastSplit[w] - cutoffs[i])/length(w)
          lastSplit[w] <- cutoffs[i]
        }
      }
    }
    numGroups <- length(unique(guideTree[, i]))
  }
}
weights <- weights + lastSplit
weights <- ifelse(weights==0, 1, weights)
weights <- weights/mean(weights)

# down-weight outlier weights
x = which(weights > (mean(weights) + 3*sd(weights)))
weights[x] = 20
weights = weights/mean(weights)

writeXStringSet(align, file='~/Desktop/PB2_align.txt')

states = PredictDBN(align, processors = 4, weight = weights, pseudo = 3)
write(states, file='~/Desktop/PB2_states.txt')
pairs = PredictDBN(align, type='pairs', processors = 4, weight = weights, pseudo = 3)
write.table(pairs, file='~/Desktop/PB2_pairs.txt')


# visualize



