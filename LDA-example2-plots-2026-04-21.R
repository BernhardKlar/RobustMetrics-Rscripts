###############################################
######## LDA and performance metrics

# Mahalanobis Abstand
Delta=1

# del: Threshold f_1/f_0 > del
condprob = function(del, Delta) {
  pi = matrix(0, length(del), 2)
  pi[,1] = pnorm((log(1/del) + Delta^2/2)/Delta)
  pi[,2] = pnorm((log(del) + Delta^2/2)/Delta)
  return(pi)
}

###############################
###### JAC - ratio of derivatives
JAC = function(pi11, pi00, pr) { #pi11, pi00: vector, pr: scalar
  return( pi11*(1-pr) / (1-pi00*(1-pr)) )
}
###### F_beta-score
Fscore = function(pi11, pi00, pr, beta) { #pi11, pi00: vector, pr: scalar
  return( pi11*(1-pr) / (beta^2*pr+(1-pi00)*(1-pr)) )
}
###### MCC
MCC = function(pi11,pi00, pr) { #pi11, pi00: vector, pr: scalar
  return( ( (pi11+pi00-1)*(2*pi11-1)*pr + (1-2*pi11)*pi00 + pi11 - 1 ) / 
            ( (1-2*pi00)*(pi11+pi00-1)*pr + 2*pi00^2 - 2*pi00) )
  }
###### Yule
Yule = function(pi11, pi00, pr) { #pi11, pi00: vector, pr: scalar
  return((pi11*(1-pi11))/(pi00*(1-pi00)))
  }
###### Cohen's kappa
Kappa = function(pi11, pi00, pr) { #pi11, pi00: vector, pr: scalar
  return( (pi11+pr*(1-2*pi11)) / (1-pi00-pr*(1-2*pi00)) )
}

sol.metric = function(del, pr, Delta, metric, beta=1) {
  pi = condprob(del, Delta)
  if (metric == "JAC") erg = del - JAC(pi[,1], pi[,2], pr) else
    if (metric == "Fscore") erg = del - Fscore(pi[,1], pi[,2], pr, beta) else
      if (metric == "MCC") erg = del - MCC(pi[,1], pi[,2], pr) else
        if (metric == "Yule") erg = del - Yule(pi[,1], pi[,2], pr) else
          if (metric == "Kappa") erg = del - Kappa(pi[,1], pi[,2], pr) else
              print("metric not implemented")
  return(erg)
}

sol.metric(seq(0.1,10,0.1), 0.5, 1, "Kappa")

root.s = function(pr, Delta, metric, beta=1) { #pr: scalar
  uniroot(sol.metric, lower=0.01, upper=100, extendInt="yes", pr=pr, Delta=Delta, metric=metric, beta=beta)$root
}
root.v = function(pr, Delta, metric, beta=1) sapply(pr, root.s, Delta, metric, beta=beta) #pr: vector

#
Gmean = function(pi11, pi00) {
  pi11 / pi00
}

sol.Gmean = function(del, pr, Delta) {
  pi = condprob(del, Delta)
  del - Gmean(pi[,1], pi[,2])
}

root.Gmean.s = function(pr, Delta) {
  uniroot(sol.Gmean, lower=0.01, upper=100, extendInt="yes", pr=pr, Delta=Delta)$root
}

root.Gmean = function(pr, Delta) sapply(pr, root.Gmean.s, Delta=Delta)
###################
#tests 
root = root.s(0.4, Delta=1, metric="JAC")
root
condprob(root, Delta=1)

Delta=1
(r1 = root.v( 0.001, Delta, metric="JAC"))
root.v( 0.001, Delta, metric="JAC")
condprob(r1, Delta)
root.v( 0.001, Delta, metric="Fscore", beta=3)

root.v( 1e-5, Delta=1, metric="MCC")
root.v( 1-1e-6, Delta=1, metric="MCC")
root.v( 0.016, Delta=2, metric="MCC")
root.v( 1-1e-8, Delta=2, metric="MCC")

#pr=0.4; del= seq(0.001,5,0.05)
#sol.metric(del, pr, Delta, "JAC")
#root = uniroot(sol.metric, lower=0.001, upper=10, extendInt="yes", pr=pr, Delta=Delta, metric="JAC")$root
#root
#condprob(root, Delta)
###############################
par(mfrow=c(1,2))
Delta = 1
curve(root.v(x, Delta, metric="JAC"), from=0.001, to=0.999, ylim=c(0,10), lwd=3, lty=1, col=1, n=500,
      xlab=expression(pi), ylab=expression(delta^"*"), main=expression(Delta==1) )
curve(root.v(x, Delta, metric="MCC"), add=TRUE, lwd=3, lty=2, col=2, n=500)
curve(root.v(x, Delta, metric="Yule"), add=TRUE, lwd=3, lty=3, col=3, n=500)
curve(root.v(x, Delta, metric="Fscore", beta=1.5), add=TRUE, lwd=3, lty=4, col=4, n=500)
curve(root.v(x, Delta, metric="Fscore", beta=0.5), add=TRUE, lwd=3, lty=5, col=5, n=500)
curve(root.v(x, Delta, metric="Kappa"), add=TRUE, lwd=3, lty=6, col=6, n=500)
#curve(root.Gmean(x, Delta), add=TRUE, lwd=3, lty=7, col=7, n=500)
legend("topright", legend=c("JAC     ", "MCC", "Yule / BA / G-mean", expression(paste(F[1.5], "-score")), 
          expression(paste(F[0.5], "-score")), "Kappa"), lwd=3, lty=1:6, col=1:6)

Delta = 2
curve(root.v(x, Delta, metric="JAC"), from=0.001, to=0.999, ylim=c(0,10), lwd=3, lty=1, col=1, n=500,
      xlab=expression(pi), ylab=expression(delta^"*"), main=expression(Delta==2) )
curve(root.v(x, Delta, metric="MCC"), add=TRUE, lwd=3, lty=2, col=2, n=500)
curve(root.v(x, Delta, metric="Yule"), add=TRUE, lwd=3, lty=3, col=3, n=500)
curve(root.v(x, Delta, metric="Fscore", beta=1.5), add=TRUE, lwd=3, lty=4, col=4, n=500)
curve(root.v(x, Delta, metric="Fscore", beta=0.5), add=TRUE, lwd=3, lty=5, col=5, n=500)
curve(root.v(x, Delta, metric="Kappa"), add=TRUE, lwd=3, lty=6, col=6, n=500)
#curve(root.Gmean(x, Delta), add=TRUE, lwd=3, lty=7, col=7, n=500)
legend("topright", legend=c("JAC     ", "MCC", "Yule / BA / G-mean", expression(paste(F[1.5], "-score")), 
                            expression(paste(F[0.5], "-score")), "Kappa"), lwd=3, lty=1:6, col=1:6)
