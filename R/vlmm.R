# unexported functions for estimating a variable length Markov model (VLMM).
#
# Here we implement the 21 bp VLMM used in Cufflinks to estimate read start
# biases. The method is described in:
#
# Roberts et al, "Improving RNA-Seq expression estimates by correcting for fragment bias"
# Genome Biology (2011) doi:101186/gb-2011-12-3-r22

alphafun <- function(x, order) {
  if (order == 0) {
    return(x)
  } else{
    return(as.vector(t(outer( alphafun(x, order-1), x, paste0))))
  }
}
getKmerFreqs <- function(seqs, dna.letters, order, pc=1) {
  alpha <- alphafun(dna.letters, order)
  n <- sum(width(seqs)) - order*length(seqs)
  if (order > 0) {
    out <- sapply(alpha, function(p) sum(vcountPattern(p, seqs)))
  } else {
    out <- colSums(letterFrequency(seqs, dna.letters))
  }
  stopifnot(sum(out) == n)
  out <- out + pc
  out/sum(out)
}
getPositionalKmerFreqs <- function(seqs, dna.letters, order, pos, pc=1) {
  alpha <- alphafun(dna.letters, order)
  out <- as.numeric(table(factor(substr(seqs, pos-order, pos), alpha)))
  names(out) <- alpha
  out <- out + pc
  out/sum(out)
}
getPositionalObsOverExp <- function(seqs, gene.seqs, dna.letters, order, pos) { 
  npos <- length(pos)
  res <- sapply(pos, function(i) getPositionalKmerFreqs(seqs,
                                   dna.letters, order=order, pos=i))
  # 'obs' is a 3 dimensional array:
  #   1st dim: A,C,G,T the current base
  #   2nd dim: the 4^order previous bases
  #   3rd dim: the position, which is within a subset of the full VLMM order
  obs <- array(0, dim=c(4, 4^order, npos), dimnames=
               list(dna.letters, alphafun(dna.letters, order-1), seq_len(npos)))  
  for (i in 1:npos) {
    obs[,,i] <- res[,i]
    obs[,,i] <- sweep(obs[,,i], 2, colSums(obs[,,i]), "/")
  }
  res.gene <- getKmerFreqs(gene.seqs, dna.letters, order)
  alpha <- alphafun(dna.letters, order-1)
  # 'expect' has dims 1 and 2 from above
  expect <- array(0, dim=c(4, 4^order),
                  dimnames=list(dna.letters, alphafun(dna.letters, order-1)))
  for (p in alpha) {
    prob <- res.gene[grep(paste0("^",p), names(res.gene))]
    expect[,p] <- prob/sum(prob)
  }
  list(obs=obs, expect=expect)
}
fitVLMM <- function (seqs, gene.seqs) {
  # fit a VLMM according to Roberts et al. (2011), doi:101186/gb-2011-12-3-r22
  dna.letters <- c("A","C","G","T")
  vlmm.order <- c(0,0,0,0,1,1,1,2,2,2,2,2,2,2,2,2,2,1,1,0,0)
  # sum(table(vlmm.order) * c(4,4^2,4^3)) # the 744 parameters
  # order 0
  order0 <- list()
  order0$obs <- sapply(seq_along(vlmm.order), function(i)
    getPositionalKmerFreqs(seqs, dna.letters, order=0, pos=i))
  colnames(order0$obs) <- seq_along(vlmm.order)
  order0$expect <- getKmerFreqs(gene.seqs, dna.letters, 0)
  # order 1
  order <- 1
  pos1 <- which(vlmm.order >= order)
  order1 <- getPositionalObsOverExp(seqs, gene.seqs, dna.letters, order, pos1)
  # order 2
  order <- 2
  pos2 <- which(vlmm.order >= order)
  order2 <- getPositionalObsOverExp(seqs, gene.seqs, dna.letters, order, pos2)
  list(order0=order0, order1=order1, order2=order2)
}
calcVLMMBias <- function(seqs, vlmm.model, short=FALSE, pseudocount=1) {
  stopifnot(!is.null(vlmm.model))
  dna.letters <- c("A","C","G","T")
  vlmm.order <- if (short) {
    # short = the VLMM when the reads are 8 or less positions from the end of transcript
    c(0,1,2,2,2,2,2,2,2,1,1,0,0)
  } else {
    c(0,0,0,0,1,1,1,2,2,2,2,2,2,2,2,2,2,1,1,0,0)
  }
  # maps from position in the seq to the VLMM matrices
  map <- if (short) {
    list("order0"=c(9:21),
         "order1"=c(rep(NA,1),6:15,rep(NA,2)),
         "order2"=c(rep(NA,2),4:10,rep(NA,4)))
  } else {
    list("order0"=1:21,
         "order1"=c(rep(NA,4),1:15,rep(NA,2)),
         "order2"=c(rep(NA,7),1:10,rep(NA,4)))
  }
  bias <- matrix(NA, length(seqs), length(vlmm.order))
  for (i in seq_along(vlmm.order)) {
    order <- vlmm.order[i]
    alpha <- alphafun(dna.letters, order)
    o <- paste0("order",order)
    kmer <- substr(seqs, i - order, i)
    j <- map[[o]][i]
    if (order == 0) {
      bias.lookup <- vlmm.model[[o]]$obs[,j] / vlmm.model[[o]]$exp
    } else {
      bias.lookup <- as.vector(vlmm.model[[o]]$obs[,,j] / vlmm.model[[o]]$exp)
      names(bias.lookup) <- alpha
    }
    bias[,i] <- bias.lookup[ kmer ]
  }
  bias
}
addVLMMBias <- function(fragtypes, vlmm.fivep, vlmm.threep) {
  fivep <- fragtypes$fivep[fragtypes$fivep.test]
  fivep.short <- fragtypes$fivep[!fragtypes$fivep.test]
  threep <- fragtypes$threep[fragtypes$threep.test]
  threep.short <- fragtypes$threep[!fragtypes$threep.test]   
  fivep.bias <- numeric(nrow(fragtypes))
  fivep.bias[fragtypes$fivep.test] <- rowSums(log(calcVLMMBias(fivep, vlmm.fivep, short=FALSE)))
  fivep.bias[!fragtypes$fivep.test] <- rowSums(log(calcVLMMBias(fivep.short, vlmm.fivep, short=TRUE)))
  fragtypes$fivep.bias <- fivep.bias
  threep.bias <- numeric(nrow(fragtypes))
  threep.bias[fragtypes$threep.test] <- rowSums(log(calcVLMMBias(threep, vlmm.threep, short=FALSE)))
  threep.bias[!fragtypes$threep.test] <- rowSums(log(calcVLMMBias(threep.short, vlmm.threep, short=TRUE)))
  fragtypes$threep.bias <- threep.bias
  fragtypes
}
