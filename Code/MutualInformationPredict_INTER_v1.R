# source("https://bioconductor.org/biocLite.R")
# biocLite("DECIPHER")

library('DECIPHER')

rna_seg1 = readRNAStringSet("~Desktop/")

rna_seg2 = readRNAStringSet("~/Desktop/")

# Input values
thresh <- 0.3
minOccupancy <- 0.5
penalty <- 1 # inconsistent pairs penalty
weight <- c(rep(1, 2), rep(0.4, 2), rep(1.2, 2), 0) # relative weight of A/U, G/U, G/C, -/-
lambda <- 0
coeff <- 1
Beta <- 0

#####################################################
## Calculate sequence weights
#####################################################

# get weights (normally would be supplied)
d <- DistanceMatrix(rna_seg1, verbose=FALSE, includeTerminalGaps=TRUE)    # How dissimilar is on seq from another?
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
# For every cutoff
for (i in (dim(guideTree)[2] - 1):1) {
  # If this cutoff produces a group greater than the previous number of groups
  if (length(unique(guideTree[, i])) > numGroups) {
    if (numGroups==1) { # root of tree
      lastSplit[] <- cutoffs[i]
    } else {
      # assigning weights based on which group they're in
      for (j in unique(guideTree[, i + 1])) {
        w <- which(guideTree[, i + 1]==j)
        if (length(unique(guideTree[w, i])) > 1) {
          weights[w] <- weights[w] + (lastSplit[w] - cutoffs[i])/length(w) #WHAT IS THIS?
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

################################################################# Need weights for each segment individually

d <- DistanceMatrix(rna_seg2, verbose=FALSE, includeTerminalGaps=TRUE)    # How dissimilar is on seq from another?
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
weights2 <- numeric(dim(guideTree)[1])
numGroups <- 1
# For every cutoff
for (i in (dim(guideTree)[2] - 1):1) {
  # If this cutoff produces a group greater than the previous number of groups
  if (length(unique(guideTree[, i])) > numGroups) {
    if (numGroups==1) { # root of tree
      lastSplit[] <- cutoffs[i]
    } else {
      # assigning weights based on which group they're in
      for (j in unique(guideTree[, i + 1])) {
        w <- which(guideTree[, i + 1]==j)
        if (length(unique(guideTree[w, i])) > 1) {
          weights2[w] <- weights2[w] + (lastSplit[w] - cutoffs[i])/length(w) #WHAT IS THIS?
          lastSplit[w] <- cutoffs[i]
        }
      }
    }
    numGroups <- length(unique(guideTree[, i]))
  }
}
weights2 <- weights2 + lastSplit
weights2 <- ifelse(weights2==0, 1, weights2)
weights2 <- weights2/mean(weights2)

# The weights associated with each seq

#####################################################
## Record the total tree length
#####################################################

suppressWarnings(guideTree <- IdClusters(d,
                                         method="UPGMA",
                                         type="dendrogram",
                                         verbose=FALSE,
                                         collapse=-1))
height <- 0
.treeL <- function(dend) {
  if (is.leaf(dend))
    return(NULL)
  h <- attr(dend, "height")
  s <- sapply(dend, function(x) attr(x, "height"))
  height <<- height + sum(h - s)
  for (i in seq_along(dend))
    .treeL(dend[[i]])
  return(NULL)
}
.treeL(guideTree)

#####################################################
## Calculate column base frequencies
#####################################################

a <- colSums(alphabetFrequency(rna, baseOnly=TRUE))[1:4]
a <- a/sum(a)

r <- strsplit(as.character(rna), "", fixed=TRUE)
r <- matrix(unlist(r), length(rna), byrow=TRUE)

f <- matrix(0, ncol(r), 4, dimnames=list(NULL, names(a)))
for (i in 1:ncol(r)) {
  t <- tapply(weights, r[, i], sum)[names(a)]
  if (any(is.na(t))) {
    t[is.na(t)] <- 0
  }
  f[i, ] <- t/nrow(r)
}
w <- which(rowSums(f) > minOccupancy)

# f = frequency by base
# w = only keep L that don't have gaps more than 50% of the time
#####################################################
## Calculate the mutual information matrix
#####################################################

paired <- c("A U", "U A",
            "G U", "U G",
            "G C", "C G",
            "- -")
unpaired <- c("A C", "C A",
              "A G", "G A",
              "C U", "U C",
              "A A", "C C", "G G", "U U")
all <- c(paired, unpaired)

m <- matrix(0, ncol(r), ncol(r)) # want this to be length(w)^2 in practice

for (i in 1:(length(w) - 1)) {
  for (j in (i + 1):length(w)) {
    bg <- c(f[w[i], "A"]*f[w[j], "U"],
            f[w[i], "U"]*f[w[j], "A"],
            f[w[i], "G"]*f[w[j], "U"],
            f[w[i], "U"]*f[w[j], "G"],
            f[w[i], "G"]*f[w[j], "C"],
            f[w[i], "C"]*f[w[j], "G"],
            1) # -/-
    
    pairs <- paste(r[, w[i]], r[, w[j]])
    
    if (Beta > 0) {
      mutations <- .Fitch(guideTree)[8]
    } else {
      mutations <- 0
    }
    
    pairs <- tapply(weights, pairs, sum)[all]
    if (any(is.na(pairs)))
      pairs[is.na(pairs)] <- 0
    pairs <- (pairs + lambda)/(length(pairs)*lambda + nrow(r))
    if (any(is.nan(pairs)))
      next
    
    #		m[w[j], w[i]] <- sum(pairs[paired]*log(pairs[paired]/bg, 2)*weight, na.rm=TRUE) - penalty*sum(pairs[unpaired], na.rm=TRUE)
    m[w[j], w[i]] <- (1 - Beta)*(sum(pairs[paired]*log(pairs[paired]/bg, 2)*weight, na.rm=TRUE) - penalty*(1 - sum(pairs[paired], na.rm=TRUE))) + Beta*(mutations/(length(pairs) - 1))
  }
}
# image(m[w, w], xlab='Sequence', ylab='Sequence', main='Mutual Information Matrix', xaxt="n", yaxt="n")
# axis(1, at=c(0, .085, .17, .254, .339, .424, .508, .593, .678, .763, .847, .932, 1), labels=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 118))
# axis(2, at=c(0, .085, .17, .254, .339, .424, .508, .593, .678, .763, .847, .932, 1), labels=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 118))

# apply Average Product Correction (APC)

n <- m[w, w]
n <- n + t(n)
diag(n) <- NA
means <- rowMeans(n, na.rm=TRUE) # MIavg by row
avg <- mean(n, na.rm=TRUE) # MIavg overall

for (i in 1:(length(w) - 1)) {
  for (j in (i + 1):length(w)) {
    m[w[j], w[i]] <- m[w[j], w[i]] - coeff*means[i]*means[j]/avg
  }
}
#image(m[w], m[w])

#####################################################
## Determine the actual structure
#####################################################

str <- readLines(paste0("~/Wright/Wuss/", "RF00005.txt"))
s <- strsplit(str, "")[[1]]
c1 <- c2 <- integer()
parens <- square <- curly <- straight <- integer()
for (j in seq_along(s)) {
  if (s[j]=="(") {
    parens <- c(parens, j)
  } else if (s[j]==")") {
    c1 <- c(c1, parens[length(parens)])
    c2 <- c(c2, j)
    length(parens) <- length(parens) - 1L
  } else if (s[j]=="[") {
    square <- c(square, j)
  } else if (s[j]=="]") {
    c1 <- c(c1, square[length(square)])
    c2 <- c(c2, j)
    length(square) <- length(square) - 1L
  } else if (s[j]=="{") {
    curly <- c(curly, j)
  } else if (s[j]=="}") {
    c1 <- c(c1, curly[length(curly)])
    c2 <- c(c2, j)
    length(curly) <- length(curly) - 1L
  } else if (s[j]=="<") {
    straight <- c(straight, j)
  } else if (s[j]==">") {
    c1 <- c(c1, straight[length(straight)])
    c2 <- c(c2, j)
    length(straight) <- length(straight) - 1L
  }
}
m[cbind(c1, c2)] <- NA

#m[which(m < thresh)] <- 0
#image(m, zlim=c(min(m, na.rm=TRUE), 2))

#####################################################
## Apply MI to p-value transformation
#####################################################

cutoffs <- seq(0, 2, 0.01)
pos <- tot <- numeric(length(cutoffs) - 1)
x <- m[cbind(c2, c1)]
for (i in 2:length(cutoffs)) {
  pos[i - 1] <- length(which(x < cutoffs[i] & x >= cutoffs[i - 1]))
  tot[i - 1] <- length(which(m < cutoffs[i] & m >= cutoffs[i - 1]))
}
#pdf(paste0("~/Desktop/temp2/", substring(ls[k], 1, nchar(ls[k]) - 3), ".pdf"), ylim=c(0, 1))
#plot(cutoffs[-1], pos/tot)

#hist(m[w,w][which(m[w,w] != 0)], breaks=100, prob=TRUE)
#curve(dnorm(x, mean=median(m[w,w][which(m[w,w] != 0)]), sd=sd(m[w,w][which(m[w,w] != 0)]) - 0.03), add=TRUE, yaxt="n", col="orange")

lower_bound <- quantile(m[w, w][lower.tri(m[w,w])],
                        1 - length(w)/2/((length(w)^2 - length(w))/2)) # fraction that could be paired
shift <- lower_bound*1.3 # empirically fitted
slope <- 2/(shift - lower_bound)
#points(cutoffs, 1/(1 + exp(-slope*(log(cutoffs) - log(shift)))), type="l")
#dev.off()

p <- 1/(1 + exp(-slope*(log(m) - log(shift))))
p[is.nan(p)] <- 0
probs <- matrix(0,
                nrow=3,
                ncol=ncol(r),
                dimnames=list(c(".", "(", ")")))
for (i in 1:(length(w) - 1)) {
  for (j in (i + 1):length(w)) {
    if (p[w[j], w[i]] > 0) {
      if (p[w[j], w[i]] > probs["(", w[i]])
        probs["(", w[i]] <- p[w[j], w[i]]
      if (p[w[j], w[i]] > probs[")", w[j]])
        probs[")", w[j]] <- p[w[j], w[i]]
    }
  }
}
x <- colSums(probs)
for (i in seq_len(ncol(r))) {
  if (x[i] < 1) {
    probs[".", i] <- 1 - x[i]
  } else {
    probs[, i] <- probs[, i]/x[i]
  }
}

probs2 <- PredictDBN(rna, "scores", weight=weights)
#plot(probs, probs2)

#####################################################
## Calculate area under the ROC curve
#####################################################

sensitivity <- PPV <- numeric(length(threshs))
for (i in seq_along(threshs)) {
  TP <- length(which(p[cbind(c2, c1)] > threshs[i]))
  FN <- length(which(p[cbind(c2, c1)] <= threshs[i]))
  FP <- (length(which(p > threshs[i])) - length(which(p[cbind(c2, c1)] > threshs[i])))
  #	sensitivity[i] <- length(which(m[cbind(c2, c1)] > threshs[i]))/length(c1)
  #	PPV[i] <- length(which(m[cbind(c2, c1)] > threshs[i]))/length(which(m > threshs[i]))
  sensitivity[i] <- TP/(TP + FN)
  if (is.nan(sensitivity[i]))
    sensitivity[i] <- 0
  PPV[i] <- TP/(TP + FP) # aka precision
  if (is.nan(PPV[i]))
    PPV[i] <- 0
}
#if (any(is.nan(PPV)))
#	PPV[is.nan(PPV)] <- 0
#plot(sensitivity, PPV, xlim=c(0, 1), ylim=c(0, 1), pch=46, col=rainbow(length(threshs)))

# AUC = delta_x*mean(y)
AUC <- sum(-diff(sensitivity)*(PPV[-1] + PPV[-length(PPV)])/2, na.rm=TRUE)
AUC

best <- which.max(2/(1/sensitivity + 1/PPV))
#points(sensitivity[best], PPV[best])
threshs[best]

#####################################################
## Predict 2ndary structure (dot bracket notation)
#####################################################

m <- p
m[upper.tri(m)] <- 0
zero <- which(m < thresh)
if (length(zero) > 0)
  m[zero] <- 0

for (d in 2:ncol(r)) { # diagonal number
  for (i in 1:(ncol(r) - d + 1)) {
    j <- i + d - 1 # i <= j
    
    if (d > 3) { # not along first diagonals
      M <- m[j - 1, i + 1] + m[j, i]
    } else {
      M <- 0
    }
    
    L <- m[j, i + 1] # i unpaired
    if (i < ncol(r)) {
      prevL <- m[i + 1, j]
      if (prevL < 0 && prevL > -1e9) {
        prevL <- prevL - 1
      } else {
        prevL <- -1
      }
    } else {
      prevL <- -1
    }
    
    R <- m[j - 1, i] # j unpaired
    if (j > 1) {
      prevR <- m[i, j - 1]
      if (prevR < -1e9) {
        prevR <- prevR - 1
      } else {
        prevR <- -1e9 - 1
      }
    } else {
      prevR <- -1e9 - 1
    }
    
    if (M > L && M > R) {
      m[j, i] <- M
    } else if (L > R) {
      m[j, i] <- L
      m[i, j] <- prevL
    } else {
      m[j, i] <- R
      m[i, j] <- prevR
    }
    
    # bifurcation
    k <- i + 3
    while (k <= (j - 4)) {
      if (m[k, i] + m[j, k + 1] > m[j, i]) {
        m[j, i] <- m[k, i] + m[j, k + 1]
        m[i, j] <- k + 1e9
      }
      k <- k + 1
    }
    
    #		if (d==3)
    #			cat("\ni = ", i, " j = ", j, " ", m[j, i], " ", m[i, j], sep="")
  }
}

# traceback
traceback <- function(i, j) {
  while (j > (i + 1)) { # j > i allows 0 base hairpins
    cat("\ni =", i, "j =", j)
    if (m[i, j] > 1e9) { # bifurcation
      cat("\njump to", m[i, j] - 1e9 + 1, j)
      traceback(m[i, j] - 1e9 + 1, j)
      cat("\njump to", i, m[i, j] - 1e9)
      traceback(i, m[i, j] - 1e9)
      break
    } else if (m[i, j] < 0 && m[i, j] > -1e9) {
      i <- i + 1
    } else if (m[i, j] < -1e9) {
      j <- j - 1
    } else { # base pairing
      if (struct[i] !=  ".")
        stop("error i!")
      if (struct[j] !=  ".")
        stop("error j!")
      struct[i] <<- "("
      struct[j] <<- ")"
      #			m[j, i] <<- Inf
      j <- j - 1
      i <- i + 1
    }
  }
}

struct <- character(ncol(r))
struct[] <- "."
traceback(1, ncol(r)) # start at corner
struct
#image(m)

# looping start
AUCs[k] <- AUC
Ws[k] <- length(w)
Qs[k] <- quantile(m[w, w][lower.tri(m[w,w])], 1 - Ws[k]/2/((Ws[k]^2 - Ws[k])/2))
results[[k]] <- pos/tot
heights[k] <- height
threshold[k] <- threshs[best]
sensitivities <- sensitivities + sensitivity
PPVs <- PPVs + PPV
corrs[k] <- cor(as.vector(probs), as.vector(probs2))
}
z <- matrix(unlist(results), nrow=length(cutoffs) - 1)
plot(cutoffs[-1], rowMeans(z, na.rm=TRUE))
p <- 1/(1 + exp(-slope*(log(cutoffs[-1]) - log(shift))))
points(cutoffs[-1], p, type="l")
# looping end

