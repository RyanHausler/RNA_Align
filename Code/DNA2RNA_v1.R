# This function changes 5' DNA to 3' RNA in order to study H1N1 secondary structure

library(DECIPHER)

p_dna  <- readLines("~/Desktop/NS.fasta") # FASTA is 5' DNA
p_rna  <- gsub(pattern = "T", replace = "U", x = p_dna) # Change to 5' RNA
writeLines(p_rna, con="~/Desktop/NS_RNA.fasta") # Write out to new file
p_rna = readRNAStringSet('~/Desktop/NS_RNA.fasta') # Read in RNA FASTA
n_rna = reverseComplement(p_rna) # Change to 3' RNA

align = AlignSeqs(n_rna[1:5000], processors = 4) # Align 3' RNA strands

# Get weights for PredictDBN
# get weights (normally would be supplied)
d <- DistanceMatrix(align[1:5000], verbose=FALSE, includeTerminalGaps=TRUE)
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

states = PredictDBN(align, processors = 4, weight = weights)
save(states, file='~/Desktop/NS_states.RData')