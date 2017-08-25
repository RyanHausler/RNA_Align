# This function changes 5' DNA to 3' RNA in order to study H1N1 secondary structure

library(DECIPHER)

p_dna  <- readLines("~/Desktop/NP.fasta") # FASTA is 5' DNA
p_rna  <- gsub(pattern = "T", replace = "U", x = tx) # Change to RNA
writeLines(tx2, con="~/Desktop/NP_RNA.fasta") # Write out to new file
p_rna = readRNAStringSet('~/Desktop/NP_RNA.fasta') # Read in RNA FASTA
n_rna = reverseComplement(p_rna) # Change to 3' RNA

align = AlignSeqs(n_rna[1:5000], processors = 4) # Align 3' RNA strands

# Get weights for PredictDBN
# get weights (normally would be supplied)
d <- DistanceMatrix(align[1:5000], verbose=FALSE, includeTerminalGaps=TRUE)
cutoffs <- sort(unique(c(seq(0, 1, 0.01),
                         round(quantile(d, seq(0, 1, 0.01), na.rm=TRUE), 5))))
dimnames(d) <- NULL
# Clusters similar seqs
# Makes phylogenetic tree
# Each row is a seq. Each column is a cutoff.
# An index is what group each seq in
suppressWarnings(guideTree <- IdClusters(d,
                                         method="UPGMA",
                                         cutoff=cutoffs,
                                         verbose=FALSE))
guideTree <- guideTree
# Holds the cutoff at each split
lastSplit <- numeric(dim(guideTree)[1])
weights <- numeric(dim(guideTree)[1])
numGroups <- 1

for (i in (dim(guideTree)[2] - 1):1) {
  if (length(unique(guideTree[, i])) > numGroups) {
    if (numGroups==1) { # root of tree
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
save(states, file='~/Desktop/NP_states.RData')