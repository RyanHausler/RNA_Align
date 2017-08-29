library(DECIPHER)

dna <- readDNAStringSet("~/Library/Containers/com.apple.TextEdit/Data/Downloads/GenomicFastaResults.fasta")

cds <- readDNAStringSet("~/Downloads/CdsFastaResults.fasta")

clusters <- IdClusters(myXStringSet=subseq(cds, 1, 200),
	method="inexact",
	cutoff=c(0.01, 0.02, 0.04, 0.08))
a <- apply(clusters, 2, function(x) length(unique(x))) # choose a cutoff
reps <- which(!duplicated(clusters[, 1]))
dna <- dna[reps]
cds <- cds[reps]

starts <- ends <- integer(length(dna))

for (i in seq_along(dna)) {
	v <- matchPattern(cds[[i]], dna[[i]])
	starts[i] <- start(v)
	ends[i] <- end(v)
}

region1 <- subseq(dna, 1, starts - 1)
region3 <- subseq(dna, ends + 1, width(dna))

region1 <- AlignSeqs(region1, terminalGap=c(-5, -1000))
cds <- AlignTranslation(cds)
region3 <- AlignSeqs(region3, terminalGap=c(-1000, -5))

align <- xscat(region1, cds, region3)

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

p <- PredictDBN(align, type="pairs", processors=8, pseudoknots=3)

# visualize



