###############################################
######## LDA and performance metrics

# del: Threshold f_1/f_0 > del
condprob = function(del, Delta) {
  pi = matrix(0, length(del), 2)
  pi[,1] = pnorm((log(1/del) + Delta^2/2)/Delta)
  pi[,2] = pnorm((log(del) + Delta^2/2)/Delta)
  return(pi)
}

###############################
###### F_gen-score
Fscore = function(pi11, pi00, pr, coeff) { #pi11, pi00: vector, pr: scalar, coeff=c(b,d1,d2,d3)
  return( (coeff[1]+pi11)*(1-pr) / (coeff[2] + (coeff[3]-coeff[1]*coeff[4])*pr + (1-pi00)*(1-pr)) )
}

sol.F = function(del, pr, Delta, coeff) {
  pi = condprob(del, Delta)
  erg = del - Fscore(pi[,1], pi[,2], pr, coeff)
  return(erg)
}

root.F.s = function(pr, Delta, coeff) { #pr: scalar
  uniroot(sol.F, lower=0.01, upper=100, extendInt="yes", pr=pr, Delta=Delta, coeff=coeff)$root
}
root.F = function(pr, Delta, coeff) sapply(pr, root.F.s, Delta=Delta, coeff=coeff) #pr: vector

#test 
coeff = c(0,0.1,4,1)
root = root.F( c(1e-5,0.4), Delta=Delta, coeff=coeff)
root
condprob(root, Delta=1)

###############################
par(mfrow=c(1,2))
Delta = 1
# c=0, d0=0, d1=1: F-score
coeff = c(0,0,1,1)
curve(root.F(x, Delta, coeff), from=0.001, to=0.99, ylim=c(0,4.5), lwd=3, lty=1, col=1, n=500,
      xlab=expression(pi), ylab=expression(delta^"*"), main=expression(Delta==1) )
# c=0, d0 variiert, d1=1
coeff = c(0,0.1,1,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=2, col=1, n=500)
coeff = c(0,0.4,1,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=3, col=1, n=500)
coeff = c(0,1,1,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=4, col=1, n=500)
# c=0, d0=0.1, d1 variiert
coeff = c(0,0.1,0.25,1) 
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=1, col=2, n=500)
coeff = c(0,0.1,0.5,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=2, col=2, n=500)
coeff = c(0,0.1,2,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=3, col=2, n=500)
coeff = c(0,0.1,4,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=4, col=2, n=500)
# c=1, d0 variiert, d1=1 (d1=4)
coeff = c(1,1,1,1) 
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=1, col=4, n=500)
coeff = c(1,2,1,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=2, col=4, n=500)
coeff = c(1,0.5,1,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=3, col=4, n=500)
coeff = c(1,0.5,4,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=4, col=4, n=500)
#
legend("topright", legend=expression( list(c==0, d[0]==0, d[1]==1), list(c==0, d[0]==0.1, d[1]==1), 
          list(c==0, d[0]==0.4, d[1]==1), list(c==0, d[0]==1, d[1]==1),
          list(c==0, d[0]==0.1, d[1]==0.25), list(c==0, d[0]==0.1, d[1]==0.5), 
          list(c==0, d[0]==0.1, d[1]==2), list(c==0, d[0]==0.1, d[1]==4),
          list(c==1, d[0]==1, d[1]==1), list(c==1, d[0]==2, d[1]==1), 
          list(c==1, d[0]==0.5, d[1]==1), list(c==1, d[0]==0.5, d[1]==4)),
        lwd=3, lty=c(1:4,1:4,1:4), col=c(rep(1,4),rep(2,4),rep(4,4)) )
#
Delta = 2
# c=0, d0=0, d1=1: F-score
coeff = c(0,0,1,1)
curve(root.F(x, Delta, coeff), from=0.001, to=0.99, ylim=c(0,4.5), lwd=3, lty=1, col=1, n=500,
      xlab=expression(pi), ylab=expression(delta^"*"), main=expression(Delta==2) )
# c=0, d0 variiert, d1=1
coeff = c(0,0.1,1,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=2, col=1, n=500)
coeff = c(0,0.4,1,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=3, col=1, n=500)
coeff = c(0,1,1,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=4, col=1, n=500)
# c=0, d0=0.1, d1 variiert
coeff = c(0,0.1,0.25,1) 
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=1, col=2, n=500)
coeff = c(0,0.1,0.5,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=2, col=2, n=500)
coeff = c(0,0.1,2,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=3, col=2, n=500)
coeff = c(0,0.1,4,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=4, col=2, n=500)
# c=1, d0 variiert, d1=1 (d1=4)
coeff = c(1,1,1,1) 
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=1, col=4, n=500)
coeff = c(1,2,1,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=2, col=4, n=500)
coeff = c(1,0.5,1,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=3, col=4, n=500)
coeff = c(1,0.5,4,1)
curve(root.F(x, Delta, coeff), add=TRUE, lwd=3, lty=4, col=4, n=500)
#
legend("topright", legend=expression( list(c==0, d[0]==0, d[1]==1), list(c==0, d[0]==0.1, d[1]==1), 
                                      list(c==0, d[0]==0.4, d[1]==1), list(c==0, d[0]==1, d[1]==1),
                                      list(c==0, d[0]==0.1, d[1]==0.25), list(c==0, d[0]==0.1, d[1]==0.5), 
                                      list(c==0, d[0]==0.1, d[1]==2), list(c==0, d[0]==0.1, d[1]==4),
                                      list(c==1, d[0]==1, d[1]==1), list(c==1, d[0]==2, d[1]==1), 
                                      list(c==1, d[0]==0.5, d[1]==1), list(c==1, d[0]==0.5, d[1]==4)),
       lwd=3, lty=c(1:4,1:4,1:4), col=c(rep(1,4),rep(2,4),rep(4,4)) )
#####################################################################################################
#####################################################################################################
###### MCC, Delta = 1, 2
MCC = function(pi11,pi00, pr, beta) { #pi11, pi00: vector, pr: scalar
  return( ( (2*pi11-1)*pi00 - (pi11+pi00-1)*(2*pi11-1)*pr - pi11 + 1 + 2*beta) / 
            ( (2*pi00-1)*(pi11+pi00-1)*pr + 2*pi00*(1-pi00) + 2*beta) )
}

sol.metric = function(del, pr, Delta, metric, beta=1) {
  pi = condprob(del, Delta)
  if (metric == "JAC") erg = del - JAC(pi[,1], pi[,2], pr) else
    if (metric == "Fscore") erg = del - Fscore(pi[,1], pi[,2], pr, beta) else
      if (metric == "MCC") erg = del - MCC(pi[,1], pi[,2], pr, beta) else
        if (metric == "Yule") erg = del - Yule(pi[,1], pi[,2], pr) else
          if (metric == "Kappa") erg = del - Kappa(pi[,1], pi[,2], pr) else
            print("metric not implemented")
          return(erg)
}

root.s = function(pr, Delta, metric, beta=1) { #pr: scalar
  uniroot(sol.metric, lower=0.1, upper=10, extendInt="yes", pr=pr, Delta=Delta, metric=metric, beta=beta)$root
}
root.v = function(pr, Delta, metric, beta=1) sapply(pr, root.s, Delta, metric, beta=beta) #pr: vector

par(mfrow=c(1,2))
Delta = 1
curve(root.v(x, Delta, metric="MCC", beta=0), from=0.001, to=0.999, ylim=c(0,4), lwd=3, lty=1, col=1, n=500,
      xlab=expression(pi), ylab=expression(delta^"*"), main=expression(Delta==1) )
curve(root.v(x, Delta, metric="MCC", beta=0.01), add=TRUE, lwd=3, lty=2, col=2, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.05), add=TRUE, lwd=3, lty=3, col=3, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.1), add=TRUE, lwd=3, lty=4, col=4, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.5), add=TRUE, lwd=3, lty=5, col=5, n=500)
curve(root.v(x, Delta, metric="MCC", beta=1.0), add=TRUE, lwd=3, lty=6, col=6, n=500)
abline(v=0); abline(h=0); 
legend("topright", legend=c("MCC     ", expression(MCC[0.01]), expression(MCC[0.05]),  
  expression(MCC[0.1]), expression(MCC[0.5]), expression(MCC[1.0])), lwd=3, lty=1:6, col=1:6)

Delta = 2
curve(root.v(x, Delta, metric="MCC", beta=0), from=0.001, to=0.999, ylim=c(0,12), lwd=3, lty=1, col=1, n=500,
      xlab=expression(pi), ylab=expression(delta^"*"), main=expression(Delta==2) )
curve(root.v(x, Delta, metric="MCC", beta=0.01), add=TRUE, lwd=3, lty=2, col=2, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.05), add=TRUE, lwd=3, lty=3, col=3, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.1), add=TRUE, lwd=3, lty=4, col=4, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.5), add=TRUE, lwd=3, lty=5, col=5, n=500)
curve(root.v(x, Delta, metric="MCC", beta=1.0), add=TRUE, lwd=3, lty=6, col=6, n=500)
abline(v=0); abline(h=0)
legend("topright", legend=c("MCC     ", expression(MCC[0.01]), expression(MCC[0.05]),  
                            expression(MCC[0.1]), expression(MCC[0.5]), expression(MCC[1.0])), lwd=3, lty=1:6, col=1:6)
#####################################################################################################
#####################################################################################################
###### MCC, Delta = 0.5, 4  
Delta = 0.5
curve(root.v(x, Delta, metric="MCC", beta=0), from=0.001, to=0.999, ylim=c(0,2), lwd=3, lty=1, col=1, n=500,
      xlab=expression(pi), ylab=expression(delta^"*"), main=expression(Delta==0.5) )
curve(root.v(x, Delta, metric="MCC", beta=0.01), add=TRUE, lwd=3, lty=2, col=2, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.05), add=TRUE, lwd=3, lty=3, col=3, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.1), add=TRUE, lwd=3, lty=4, col=4, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.5), add=TRUE, lwd=3, lty=5, col=5, n=500)
curve(root.v(x, Delta, metric="MCC", beta=1.0), add=TRUE, lwd=3, lty=6, col=6, n=500)
abline(v=0); abline(h=0)
legend("topright", legend=c("MCC     ", expression(MCC[0.01]), expression(MCC[0.05]),  
                            expression(MCC[0.1]), expression(MCC[0.5]), expression(MCC[1.0])), lwd=3, lty=1:6, col=1:6)

Delta = 4
curve(root.v(x, Delta, metric="MCC", beta=0), from=0.001, to=0.999, ylim=c(0,50), lwd=3, lty=1, col=1, n=500,
      xlab=expression(pi), ylab=expression(delta^"*"), main=expression(Delta==4) )
curve(root.v(x, Delta, metric="MCC", beta=0.01), add=TRUE, lwd=3, lty=2, col=2, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.05), add=TRUE, lwd=3, lty=3, col=3, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.1), add=TRUE, lwd=3, lty=4, col=4, n=500)
curve(root.v(x, Delta, metric="MCC", beta=0.5), add=TRUE, lwd=3, lty=5, col=5, n=500)
curve(root.v(x, Delta, metric="MCC", beta=1.0), add=TRUE, lwd=3, lty=6, col=6, n=500)
abline(v=0); abline(h=0)
legend("topright", legend=c("MCC     ", expression(MCC[0.01]), expression(MCC[0.05]),  
                            expression(MCC[0.1]), expression(MCC[0.5]), expression(MCC[1.0])), lwd=3, lty=1:6, col=1:6)
#################################################################################################
