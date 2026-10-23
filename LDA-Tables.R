###############################################
## LDA and performance metrics - Tables

# del: Threshold f_1/f_0 > del
condprob = function(del, Delta) {
  pi = matrix(0, length(del), 2)
  pi[,1] = pnorm((log(1/del) + Delta^2/2)/Delta)
  pi[,2] = pnorm((log(del) + Delta^2/2)/Delta)
  return(pi)
}

#####################################################################################################
#####################################################################################################
## Table JAC, Yule Kappa

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
  mcc = ( (pi11+pi00-1)*(2*pi11-1)*pr + (1-2*pi11)*pi00 + pi11 - 1 ) / 
    ( (1-2*pi00)*(pi11+pi00-1)*pr + 2*pi00^2 - 2*pi00)
  if ( is.nan(mcc) ) mcc = 1e10
  return(mcc)
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

root.s = function(pr, Delta, metric, beta=1, lower=0.01, upper=1000) { #pr: scalar
  uniroot(sol.metric, lower=lower, upper=upper, extendInt="yes", pr=pr, Delta=Delta, metric=metric, beta=beta)$root
}
root.v = function(pr, Delta, metric, beta=1, lower=0.01, upper=1000) {
  sapply(pr, root.s, Delta, metric, beta=beta) #pr: vector
}

#Table
x = c(1e-4,1e-3,0.01,0.1,0.9,0.99,1-1e-3, 1-1e-4)
m = 2*length(x) + 2
root.m = matrix(nrow=2*4*5, ncol=m)

i=1
for (Delta in c(0.5,1,2,4)) {
    root = signif( root.v(x, Delta, metric="JAC", lower=1e-4, upper=10), 3)
    root.m[i,2] = i %% 5
    root.m[i,3:m] =  as.vector( rbind(root, rep(NA, length(root))) )
    root.m[i+1,3:m] = as.vector( t( condprob(root, Delta=Delta) ) )
    i = i + 2
    root = signif( root.v(x, Delta, metric="Yule"), 3)
    root.m[i,2] = i %% 5
    root.m[i,3:m] =  as.vector( rbind(root, rep(NA, length(root))) )
    root.m[i+1,3:m] = as.vector( t( condprob(root, Delta=Delta) ) )
    i = i + 2
    root = signif( root.v(x, Delta, metric="Fscore", beta=1.5), 3)
    root.m[i,2] = i %% 5
    root.m[i,3:m] =  as.vector( rbind(root, rep(NA, length(root))) )
    root.m[i+1,3:m] = as.vector( t( condprob(root, Delta=Delta) ) )
    i = i + 2
    root = signif( root.v(x, Delta, metric="Fscore", beta=0.5), 3)
    root.m[i,2] = i %% 5
    root.m[i,3:m] =  as.vector( rbind(root, rep(NA, length(root))) )
    root.m[i+1,3:m] = as.vector( t( condprob(root, Delta=Delta) ) )
    i = i + 2
    root = signif( root.v(x, Delta, metric="Kappa"), 3)
    root.m[i,2] = i %% 5
    root.m[i,3:m] =  as.vector( rbind(root, rep(NA, length(root))) )
    root.m[i+1,3:m] = as.vector( t( condprob(root, Delta=Delta) ) )
    i = i + 2
}
print( xtable(root.m), include.rownames=FALSE, digits=3) 

######################
## Test pr -> 0 JAC
x = c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,0.01,0.1)
m = length(x) + 1
root.limit = matrix(nrow=4, ncol=m)

root.limit[1,1] = 1
root.limit[1,2:m] = round( root.v2(x, Delta=1, metric="JAC", beta=0), 3)
root.limit[2,1] = 1.5
root.limit[2,2:m] = round( root.v2(x, Delta=1.5, metric="JAC", beta=0, lower=5), 3)
root.limit[3,1] = 2
root.limit[3,2:m] =round( root.v2(x, Delta=2, metric="JAC", beta=0, lower=5), 3)
root.limit[4,1] = 3
root.limit[4,2:m] = round( root.v2(x, Delta=3, metric="JAC", beta=0, lower=10, upper=1e6), 2)
print( xtable(root.limit), include.rownames=FALSE) 
######################
## Test pr -> 0 F(1.5)
root.limit[1,1] = 1
root.limit[1,2:m] = round( root.v2(x, Delta=1, metric="Fscore", beta=1.5), 3)
root.limit[2,1] = 1.5
root.limit[2,2:m] = round( root.v2(x, Delta=1.5, metric="Fscore", beta=1.5, lower=5), 3)
root.limit[3,1] = 2
root.limit[3,2:m] =round( root.v2(x, Delta=2, metric="Fscore", beta=1.5, lower=5), 3)
root.limit[4,1] = 3
root.limit[4,2:m] = round( root.v2(x, Delta=3, metric="Fscore", beta=1.5, lower=10, upper=1e6), 2)
print( xtable(root.limit), include.rownames=FALSE) 
######################
## Test pr -> 0 F(0.5)
root.limit[1,1] = 1
root.limit[1,2:m] = round( root.v2(x, Delta=1, metric="Fscore", beta=0.5), 3)
root.limit[2,1] = 1.5
root.limit[2,2:m] = round( root.v2(x, Delta=1.5, metric="Fscore", beta=0.5, lower=5), 3)
root.limit[3,1] = 2
root.limit[3,2:m] =round( root.v2(x, Delta=2, metric="Fscore", beta=0.5, lower=5), 3)
root.limit[4,1] = 3
root.limit[4,2:m] = round( root.v2(x, Delta=3, metric="Fscore", beta=0.5, lower=10, upper=1e6), 2)
print( xtable(root.limit), include.rownames=FALSE) 
######################
## Test pr -> 0 Kappa
root.limit[1,1] = 1
root.limit[1,2:m] = round( root.v2(x, Delta=1, metric="Kappa", beta=0.5), 3)
root.limit[2,1] = 1.5
root.limit[2,2:m] = round( root.v2(x, Delta=1.5, metric="Kappa", beta=0.5, lower=5), 3)
root.limit[3,1] = 2
root.limit[3,2:m] =round( root.v2(x, Delta=2, metric="Kappa", beta=0.5, lower=5), 3)
root.limit[4,1] = 3
root.limit[4,2:m] = round( root.v2(x, Delta=3, metric="Kappa", beta=0.5, lower=10, upper=1e6), 2)
print( xtable(root.limit), include.rownames=FALSE) 


#####################################################################################################
#####################################################################################################
## Table F_gen-score

Fscore = function(pi11, pi00, pr, coeff) { #pi11, pi00: vector, pr: scalar, coeff=c(b,d1,d2,d3)
  return( (coeff[1]+pi11)*(1-pr) / (coeff[2] + (coeff[3]-coeff[1]*coeff[4])*pr + (1-pi00)*(1-pr)) )
}

sol.F = function(del, pr, Delta, coeff) {
  pi = condprob(del, Delta)
  erg = del - Fscore(pi[,1], pi[,2], pr, coeff)
  return(erg)
}

root.F.s = function(pr, Delta, coeff) { #pr: scalar
  uniroot(sol.F, lower=0.0001, upper=1000, extendInt="yes", pr=pr, Delta=Delta, coeff=coeff)$root
}
root.F = function(pr, Delta, coeff) sapply(pr, root.F.s, Delta=Delta, coeff=coeff) #pr: vector

#Table
x = c(1e-4,1e-3,0.01,0.1,0.9,0.99,1-1e-3, 1-1e-4)
m = 2*length(x) + 2
root.m = matrix(nrow=2*4*4, ncol=m)

# c=0, d0 varies, d1=1, d2=1 
coeff = c(0,0,1,1)  
i=1
for (Delta in c(0.5,1,2,3)) {
  for (d in c(0,0.1,0.4,1)) {
    coeff[2] = d
    root = signif( root.F(x, Delta=Delta, coeff=coeff), 3) 
    root.m[i,2] = d
    root.m[i,3:m] =  as.vector( rbind(root, rep(NA, length(root))) )
    root.m[i+1,3:m] = as.vector( t( condprob(root, Delta=Delta) ) )
    i = i + 2
  }
}
root.m
print( xtable(root.m), include.rownames=FALSE) 

# c=0, d0=0.1, d1 varies, d2=1 
coeff = c(0,0.1,1,1)  
i=1
for (Delta in c(0.5,1,2,3)) {
  for (d in c(0.25,0.5,2,4)) {
    coeff[3] = d
    root = signif( root.F(x, Delta=Delta, coeff=coeff), 3) 
    root.m[i,2] = d
    root.m[i,3:m] =  as.vector( rbind(root, rep(NA, length(root))) )
    root.m[i+1,3:m] = as.vector( t( condprob(root, Delta=Delta) ) )
    i = i + 2
  }
}
root.m
print( xtable(root.m), include.rownames=FALSE) 

# c=1; d1=1 for the first three settings and d1=4 for the last setting
d0_values = c(1,2,0.5,0.5)
d1_values = c(1,1,1,4)
i=1
for (Delta in c(0.5,1,2,3)) {
  for (j in seq_along(d0_values)) {
    d = d0_values[j]
    coeff = c(1,d,d1_values[j],1) # (c, d0, d1, 1), reset for every setting
    root = signif( root.F(x, Delta=Delta, coeff=coeff), 3) 
    root.m[i,2] = d
    root.m[i,3:m] =  as.vector( rbind(root, rep(NA, length(root))) )
    root.m[i+1,3:m] = as.vector( t( condprob(root, Delta=Delta) ) )
    i = i + 2
  }
}
root.m
print( xtable(root.m), include.rownames=FALSE) 

#####################################################################################################
#####################################################################################################

## Table MCC ##
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

#Table
x = c(1e-4,1e-3,0.01,0.1,0.9,0.99,1-1e-3, 1-1e-4)
m = 2*length(x) + 2
root.m = matrix(nrow=2*4*6, ncol=m)

i=1
for (Delta in c(0.5,1,2,4)) {
  for (beta in c(0,0.01,0.05,0.1,0.5,1)) { #beta = d in paper
    root = signif( root.v(x, Delta=Delta, metric="MCC", beta=beta), 3) 
    root.m[i,2] = beta
    root.m[i,3:m] =  as.vector( rbind(root, rep(NA, length(root))) )
    root.m[i+1,3:m] = as.vector( t( condprob(root, Delta=Delta) ) )
    
    i = i + 2
  }
}
library(xtable)
print( xtable(root.m), include.rownames=FALSE) 

######################
## Test pr -> 0
root.s2 = function(pr, Delta, metric, beta, lower, upper) { #pr: scalar
  uniroot(sol.metric, lower=lower, upper=upper, extendInt="yes", pr=pr, Delta=Delta, metric=metric, beta=beta)$root
}
root.v2 = function(pr, Delta, metric, beta=0, lower=1, upper=1000) {
  sapply(pr, root.s2, Delta, metric, beta=beta, lower=lower, upper=upper) #pr: vector
}

x = c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,0.01,0.1)
m = length(x) + 1
root.limit = matrix(nrow=4, ncol=m)

root.limit[1,1] = 1
root.limit[1,2:m] = round( root.v2(x, Delta=1, metric="MCC", beta=0), 3)
root.limit[2,1] = 1.5
root.limit[2,2:m] = round( root.v2(x, Delta=1.5, metric="MCC", beta=0, lower=5), 3)
root.limit[3,1] = 2
root.limit[3,2:m] =round( root.v2(x, Delta=2, metric="MCC", beta=0, lower=5), 3)
root.limit[4,1] = 3
root.limit[4,2:m] = round( root.v2(x, Delta=3, metric="MCC", beta=0, lower=10, upper=1e6), 2)
print( xtable(root.limit), include.rownames=FALSE) 

