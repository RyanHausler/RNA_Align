# This function changes 5' DNA to 3' RNA in order to study H1N1 secondary structure

library(DECIPHER)

p_dna  <- readDNAStringSet("~/Desktop/Segments/PA.fasta") # FASTA is 5' DNA
p_rna = RNAStringSet(p_dna) # Read in RNA FASTA
n_rna = reverseComplement(p_rna) # Change to 3' RNA
                      # I can't align more than 5000 seqs due to 'STACK OVERFLOW ERROR'
align = AlignSeqs(n_rna[1:5000], processors = 4) # Align 3' RNA strands
ref = readDNAStringSet("~/Desktop/Segments/PA_Correct.fasta")
ref_rna = RNAStringSet(ref)

full_align = AlignProfiles(ref_rna, align)
BrowseSeqs(full_align)

# Get weights for PredictDBN
# get weights (normally would be supplied)
d <- DistanceMatrix(full_align, verbose=FALSE, includeTerminalGaps=TRUE)
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

# Find outlier sequences
x = which(weights > (mean(weights) + 3*sd(weights)))
weights[x] = 20
weights = weights/mean(weights)
save(full_align, file = '~/Desktop/Alignment/PA_align.RData')

states = PredictDBN(full_align, processors = 4, weight = weights, pseudo = 3)
cat(pairs, file='~/Desktop/Predictions/PA_states.txt')
pairs = PredictDBN(full_align, type='pairs', processors = 4, weight = weights, pseudo = 3)
save(pairs, file='~/Desktop/Predictions/PA_pairs.RData')