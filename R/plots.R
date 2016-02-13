# plot for fragment GC bias over multiple samples
# This takes the natural spline estimated coefficients from the Poisson regression
# and generates a smooth function of log fragment rate over fragment GC.
# This is similar to the plot you would get with, e.g. plot.gam()
# after having fit a generalized additive model.
plotGC <- function(fitpar, knots=c(.4,.5,.6), bk=c(0,1), col, lty, m="GC", ylim) {
  n <- length(knots)
  coefmat <- sapply(fitpar, function(elem) elem[["coefs"]][[m]][1:(n+2)])
  genecoefs <- lapply(fitpar, function(elem) elem[["coefs"]][[m]][ grep("gene", names(elem[["coefs"]][[m]])) ])
  # new intercept: the average of the intercept + gene coefficients
  coefmat[1,] <- coefmat[1,] + sapply(genecoefs, mean)
  z <- seq(from=.2,to=.8,length=101)
  x <- model.matrix(~ ns(z, knots=knots, Boundary.knots=bk))
  logpred <- x %*% coefmat
  if (missing(ylim)) {
    ylim <- c(min(logpred),max(logpred))
  }
  plot(0,0,type="n",xlim=c(.2,.8),ylim=ylim,
       ylab="log fragment rate", xlab="fragment GC content",
       main="fragment sequence effect")
  if (missing(col)) {
    col <- rep("black", ncol(logpred))
  }
  if (missing(lty)) {
    lty <- rep(1, ncol(logpred))
  }
  for (i in 1:ncol(logpred)) {
    lines(z, logpred[,i], col=col[i], lwd=2, lty=lty[i])
  }
}
# plots for VLMM (read start sequence bias) for a single sample
plotOrder0 <- function(order0, ...) {
  dna.letters <- c("A","C","G","T")
  mat <- log(order0$obs/order0$expect)
  xpos <- -8:12
  dna.cols <- c("green3","blue3","orange3","red3")
  plot(0,0,xlim=c(-8,12),ylim=c(-0.3,0.3),type="n",xlab="position",ylab="log(observed / expected)",...)
  for (i in 1:4) {
    points(xpos, mat[i,], col=dna.cols[i], type="b", lwd=2)
  }
  abline(v=0, h=0, col=rgb(0,0,0,.3))
  legend("topright",dna.letters,pch=1,lty=1,col=dna.cols,bg="white")
}
plotOrder1 <- function(order1, pos1) {
  order <- 1
  npos <- length(pos1)
  dna.cols <- c("green3","blue3","orange3","red3")
  mypar(1,npos+1,mar=c(5,1,3,1))
  for (i in 1:npos) {
    plot(as.vector(log(order1$obs[,,i]/order1$expect)), rev(seq_len(4 * 4^order)),
         col=rep(dna.cols,each=4), xlim=c(-1,1),
         ylab="",xlab="",yaxt="n",main=pos1[i] - 10 + 1,cex=2)
    abline(v=0)
  }
  plot(0,0,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  alpha <- alphafun(dna.letters, order-1)
  legend("center",alpha,pch=1,col=dna.cols,cex=2,title="prev")
}
plotOrder2 <- function(order2, pos2) {
  order <- 2
  npos <- length(pos2)
  dna.cols <- c("green3","blue3","orange3","red3")
  mypar(1,npos+1,mar=c(5,1,3,1))
  for (i in 1:npos) {
    plot(as.vector(log(order2$obs[,,i]/order2$expect)), rev(seq_len(4 * 4^order)),
         col=rep(dna.cols,each=4),pch=rep(1:4,each=16), xlim=c(-1,1),
         ylab="",xlab="",yaxt="n",main=pos2[i] - 10 + 1,cex=2)
    abline(v=0)
  }
  plot(0,0,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  alpha <- alphafun(dna.letters, order-1)
  legend("center",alpha,pch=rep(1:4,each=4),col=dna.cols,cex=2,title="prev")
}
plotRelPos <- function(fitpar, knots=c(.25,.5,.75), bk=c(0,1), col, lty, m="GC", ylim) {
  n <- length(knots)
  # assuming same number of knots for GC and for relpos
  coefmat <- sapply(fitpar, function(elem) elem[["coefs"]][[m]][c(1,(3+n):(3+2*n))])
  z <- seq(from=0,to=1,length=101)
  x <- model.matrix(~ ns(z, knots=knots, Boundary.knots=bk))
  logpred <- x %*% coefmat
  logpred <- scale(logpred, scale=FALSE)
  if (missing(ylim)) {
    ylim <- c(min(logpred),max(logpred))
  }
  plot(0,0,type="n",xlim=c(0,1),ylim=ylim,
       ylab="log fragment rate", xlab="5' -- position in transcript -- 3'",
       main="relative position bias")
  if (missing(col)) {
    col <- rep("black", ncol(logpred))
  }
  if (missing(lty)) {
    lty <- rep(1, ncol(logpred))
  }
  for (i in 1:ncol(logpred)) {
    lines(z, logpred[,i], col=col[i], lwd=2, lty=lty[i])
  }
}
