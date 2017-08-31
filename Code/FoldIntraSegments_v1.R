library(DECIPHER)

dna <- readDNAStringSet("~/Desktop/FullSegments/full_PB2.fasta")

cds <- readDNAStringSet('~/Desktop/FullSegments/cds_PB2.fasta')

ref <- readDNAStringSet('~/Desktop/Segments/PB2_Correct.fasta')
ref <- reverseComplement(ref)

x = seq(1, length(cds), 15)
clusters <- IdClusters(myXStringSet=cds[x],
	method="inexact",
	cutoff=c(0.02))
# a <- apply(clusters, 2, function(x) length(unique(x))) # choose a cutoff
reps <- which(!duplicated(clusters[, 1]))
dna <- dna[reps]
cds <- cds[reps]

y = which(width(cds) < 50)
if(length(y) > 0){
  dna = dna[-y]
  cds = cds[-y]
}

# Map DNA name to CDS name
dna_names = names(dna)
cds_names = names(cds)

for(i in seq_along(dna_names)){
  n = dna_names[i]
  pos = regexpr('Organism', n)[1] - 2
  new = substr(n, 1, pos)
  dna_names[i] = new
}
for(i in seq_along(cds_names)){
  n = cds_names[i]
  pos = gregexpr(':', n)[[1]][2] - 1
  new = substr(n, 1, pos)
  cds_names[i] = new
}
matches = match(dna_names, cds_names)
if(anyNA(matches) == TRUE){
  delete = which(is.na(matches))
  dna = dna[-delete]
  matches = matches[-delete]
}
cds = cds[matches]

# Define start and stop
starts <- ends <- integer(length(dna))
for (i in seq_along(dna)) {
  v <- matchPattern(cds[[i]], dna[[i]])
  starts[i] = start(v)
  ends[i] = end(v)
}

region1 <- subseq(dna, 1, starts - 1)
region3 <- subseq(dna, ends + 1, width(dna))

region1 <- AlignSeqs(region1, terminalGap=c(-5, -1000))
cds <- AlignTranslation(cds)
region3 <- AlignSeqs(region3, terminalGap=c(-1000, -5))

align <- xscat(region1, cds, region3)
align = AlignProfiles(ref, align)

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



states = PredictDBN(align, processors = 4, weight = weights, pseudo = 3)
write(states, file='~/Desktop/Predictions/PB2_states.txt')
pairs = PredictDBN(align, type='pairs', processors = 4, weight = weights, pseudo = 3)
write.table(pairs, file='~/Desktop/Predictions/MP_pairs.txt')

writeXStringSet(align, file='~/Desktop/Align/PB2_align.txt')
# visualize



