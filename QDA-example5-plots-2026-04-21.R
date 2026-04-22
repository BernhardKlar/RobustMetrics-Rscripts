###############################################
######## QDA and performance metrics
library(MASS); library(amen); library(CompQuadForm)
library(ISLR); library(class); library(caret)

condprob = function(del, mu0, mu1, Sigma0, Sigma1) {
  pi = matrix(0, length(del), 2)
  Sigma0.sqrt = mhalf(Sigma0); Sigma1.sqrt = mhalf(Sigma1)
  Sigma0.inv = solve(Sigma0); Sigma1.inv = solve(Sigma1)
  Sigma0.det = det(Sigma0); Sigma1.det = det(Sigma1)
  #
  spec1 = eigen(diag(nrow=2) - Sigma1.sqrt%*%Sigma0.inv%*%Sigma1.sqrt)
  Q1 = spec1$vectors
  lambda1 = spec1$values
  a1 = -Q1%*%Sigma1.sqrt%*%Sigma0.inv%*%(mu1-mu0)
  ncen1 = as.vector(a1^2/lambda1^2)
  b1 = log(Sigma0.det/Sigma1.det) + t(mu1 - mu0)%*%(Sigma0.inv)%*%(mu1 - mu0) + sum(a1^2/lambda1) + 2*log(1/del)
  pi[,1] = 1-davies(b1, lambda1, h = rep(1, length(lambda1)), delta = ncen1)$Qq
  #
  spec0 = eigen(diag(nrow=2) - Sigma0.sqrt%*%Sigma1.inv%*%Sigma0.sqrt)
  Q0 = spec0$vectors
  lambda0 = spec0$values
  a0 = -Q0%*%Sigma0.sqrt%*%Sigma1.inv%*%(mu0-mu1)
  ncen0 = as.vector(a0^2/lambda0^2)
  b0=log(Sigma1.det/Sigma0.det) + t(mu1 - mu0)%*%(Sigma1.inv)%*%(mu1 - mu0) + sum(a0^2/lambda0) + 2*log(del)
  pi[,2] = 1- davies(b0, lambda0, h = rep(1, length(lambda0)), delta = ncen0)$Qq
  #
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
###### G-mean: quotient of derivatives
Gmean = function(pi11, pi00, pr) { # pr dummy, for uniform interface
  return(pi11 / pi00)
}

sol.metric = function(del, pr, mu0, mu1, Sigma0, Sigma1, metric, beta=1) {
  pi = condprob(del, mu0, mu1, Sigma0, Sigma1)
  if (metric == "JAC") erg = del - JAC(pi[,1], pi[,2], pr) else
    if (metric == "Fscore") erg = del - Fscore(pi[,1], pi[,2], pr, beta) else
      if (metric == "MCC") erg = del - MCC(pi[,1], pi[,2], pr) else
        if (metric == "Yule") erg = del - Yule(pi[,1], pi[,2], pr) else
          if (metric == "Kappa") erg = del - Kappa(pi[,1], pi[,2], pr) else
            if (metric == "Gmean") erg = del - Gmean(pi[,1], pi[,2], pr) else
              print("metric not implemented")
            return(erg)
}

root.s = function(pr, mu0, mu1, Sigma0, Sigma1, metric, beta=1) { #pr: scalar
  uniroot(sol.metric, lower=0.01, upper=100, extendInt="yes", pr=pr,
          mu0=mu0, mu1=mu1, Sigma0=Sigma0, Sigma1=Sigma1, metric=metric, beta=beta)$root
}
root.v = function(pr, mu0, mu1, Sigma0, Sigma1, metric, beta=1) {
  sapply(pr, root.s, mu0=mu0, mu1=mu1, Sigma0=Sigma0, Sigma1=Sigma1, metric=metric, beta=beta) #pr: vector
}

#test
Sigma0=matrix(c(2,0.5,0.5,1), ncol=2); Sigma1=matrix(c(1,-0.5,-0.5,2), ncol=2)
#Sigma0=matrix(c(2,0.3,0.3,1), ncol=2); Sigma1=matrix(c(1,-0.9,-0.9,2), ncol=2)
mu0=c(0,0); mu1=c(2.5,2.5)
#mu1=c(1.5,1.5)

(r1 = root.v(0.1, mu0, mu1, Sigma0, Sigma1, metric="JAC"))
condprob(r1, mu0, mu1, Sigma0, Sigma1)

# Scenario 1: 1.03, Scenario 2: 1.12
root.v(0.5, mu0, mu1, Sigma0, Sigma1, metric="Gmean")


###############################
par(mfrow=c(1,2))
Sigma0=matrix(c(2,0.5,0.5,1), ncol=2); Sigma1=matrix(c(1,-0.5,-0.5,2), ncol=2)
mu0=c(0,0); mu1=c(2.5,2.5)
#
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="JAC"), from=0.01, to=0.99, ylim=c(0,10), lwd=3, lty=1, col=1,
      xlab=expression(pi), ylab=expression(delta^"*"), main="Szenario 1")
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="MCC"), add=TRUE, lwd=3, lty=2, col=2, n=500)
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="Yule"), add=TRUE, lwd=3, lty=3, col=3, n=500)
abline(h=1, lwd=3, lty=4, col=8)  # balanced accuracy
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="Gmean"), add=TRUE, lwd=3, lty=5, col=7, n=500)
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="Fscore", beta=1.5), add=TRUE, lwd=3, lty=6, col=4, n=500)
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="Fscore", beta=0.5), add=TRUE, lwd=3, lty=7, col=5, n=500)
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="Kappa"), add=TRUE, lwd=3, lty=8, col=6, n=500)
legend("topright", legend=c("JAC     ", "MCC", "Yule", "BA", "G-mean",
                            expression(paste(F[1.5], "-score")), expression(paste(F[0.5], "-score")), "Kappa"),
       lwd=3, lty=1:8, col=c(1,2,3,8,7,4,5,6))

Sigma0=matrix(c(2,0.3,0.3,1), ncol=2); Sigma1=matrix(c(1,-0.9,-0.9,2), ncol=2)
mu0=c(0,0); mu1=c(1.5,1.5)
#
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="JAC"), from=0.01, to=0.99, ylim=c(0,10), lwd=3, lty=1, col=1,
      xlab=expression(pi), ylab=expression(delta^"*"), main="Scenario 2")
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="MCC"), add=TRUE, lwd=3, lty=2, col=2, n=500)
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="Yule"), add=TRUE, lwd=3, lty=3, col=3, n=500)
abline(h=1, lwd=3, lty=4, col=8)  # balanced accuracy
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="Gmean"), add=TRUE, lwd=3, lty=5, col=7, n=500)
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="Fscore", beta=1.5), add=TRUE, lwd=3, lty=6, col=4, n=500)
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="Fscore", beta=0.5), add=TRUE, lwd=3, lty=7, col=5, n=500)
curve(root.v(x, mu0, mu1, Sigma0, Sigma1, metric="Kappa"), add=TRUE, lwd=3, lty=8, col=6, n=500)
legend("topright", legend=c("JAC     ", "MCC", "Yule", "BA", "G-mean",
                            expression(paste(F[1.5], "-score")), expression(paste(F[0.5], "-score")), "Kappa"),
       lwd=3, lty=1:8, col=c(1,2,3,8,7,4,5,6))
