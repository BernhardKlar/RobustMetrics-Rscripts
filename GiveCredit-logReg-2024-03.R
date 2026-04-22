#########################
## function definition ##
#########################
# del: Threshold f_1/f_0 > del
condprob = function(del, Delta) {
  pi = matrix(0, length(del), 2)
  pi[,1] = pnorm((log(1/del) + Delta^2/2)/Delta)
  pi[,2] = pnorm((log(del) + Delta^2/2)/Delta)
  return(pi)
}

#############################################################
###### JAC - ratio of derivatives
JAC = function(pi11, pi00, pr) { #pi11, pi00: vector, pr: scalar
  return( pi11*(1-pr) / (1-pi00*(1-pr)) )
}
###### F_beta-score - ratio of derivatives
Fscore = function(pi11, pi00, pr, beta) { #pi11, pi00: vector, pr: scalar
  return( pi11*(1-pr) / (beta^2*pr+(1-pi00)*(1-pr)) )
}
###### MCC - ratio of derivatives
MCC = function(pi11,pi00, pr) { #pi11, pi00: vector, pr: scalar
  return( ( (pi11+pi00-1)*(2*pi11-1)*pr + (1-2*pi11)*pi00 + pi11 - 1 ) / 
            ( (1-2*pi00)*(pi11+pi00-1)*pr + 2*pi00^2 - 2*pi00) )
}
###### Yule - ratio of derivatives
Yule = function(pi11, pi00, pr) { #pi11, pi00: vector, pr: scalar
  return((pi11*(1-pi11))/(pi00*(1-pi00)))
}
###### Cohen's kappa - ratio of derivatives
Kappa = function(pi11, pi00, pr) { #pi11, pi00: vector, pr: scalar
  return( (pi11+pr*(1-2*pi11)) / (1-pi00-pr*(1-2*pi00)) )
}
###### MCC_gen: beta is here used as robustness parameter
MCC.gen = function(pi11,pi00, pr, beta) { #pi11, pi00: vector, pr: scalar
  return( ( (2*pi11-1)*pi00 - (pi11+pi00-1)*(2*pi11-1)*pr - pi11 + 1 + 2*beta) / 
            ( (2*pi00-1)*(pi11+pi00-1)*pr + 2*pi00*(1-pi00) + 2*beta) )
}

sol.metric = function(del, pr, Delta, metric, beta=1) {
  pi = condprob(del, Delta)
  if (metric == "JAC") erg = del - JAC(pi[,1], pi[,2], pr) else
    if (metric == "Fscore") erg = del - Fscore(pi[,1], pi[,2], pr, beta) else
      if (metric == "MCC") erg = del - MCC(pi[,1], pi[,2], pr) else
        if (metric == "MCC_gen") erg = del - MCC.gen(pi[,1], pi[,2], pr, beta) else
          if (metric == "Yule") erg = del - Yule(pi[,1], pi[,2], pr) else
            if (metric == "Kappa") erg = del - Kappa(pi[,1], pi[,2], pr) else
              print("metric not implemented")
            return(erg)
}

root.s = function(pr, Delta, metric, beta=1) { #pr: scalar
  uniroot(sol.metric, lower=0.01, upper=100, extendInt="yes", pr=pr, Delta=Delta, metric=metric, beta=beta)$root
}
root.v = function(pr, Delta, metric, beta=1) sapply(pr, root.s, Delta, metric, beta=beta) #pr: vector
#############################################################

###### F_gen-score - ratio of derivatives
F.gen = function(pi11, pi00, pr, coeff) { #pi11, pi00: vector, pr: scalar, coeff=c(b,d1,d2,d3)
  return( (coeff[1]+pi11)*(1-pr) / (coeff[2] + (coeff[3]-coeff[1]*coeff[4])*pr + (1-pi00)*(1-pr)) )
}

sol.F = function(del, pr, Delta, coeff) {
  pi = condprob(del, Delta)
  erg = del - F.gen(pi[,1], pi[,2], pr, coeff)
  return(erg)
}

root.F.s = function(pr, Delta, coeff) { #pr: scalar
  uniroot(sol.F, lower=0.01, upper=100, extendInt="yes", pr=pr, Delta=Delta, coeff=coeff)$root
}
root.F = function(pr, Delta, coeff) sapply(pr, root.F.s, Delta=Delta, coeff=coeff) #pr: vector

#############################################################
## Computation of the different metrics ##
#############################################################
## F-Score 
F.metric = function(pi11, pi00, pr, beta) { #pi11, pi00: vector, pr: scalar
  return( (1+beta^2)*pr*pi11 / (beta^2*pr + pr*pi11 + (1-pr)*(1-pi00)) )
}
## MCC
MCC.metric = function(pi11,pi00, pr) { #pi11, pi00: vector, pr: scalar
  return( sqrt(pr*(1-pr)) * ( pi11*pi00 - (1-pi11)*(1-pi00) ) / 
            sqrt( (pr*pi11 + (1-pr)*(1-pi00)) * (pr*(1-pi11)+(1-pr)*pi00) ) )
}
## F_rb-Score, standard: coeff[1]=c=0, coeff[2]=d_0!=0, coeff[3]=d_1=beta^2, coeff[4]=d_2=1,
## with normalization
Fgen.metric = function(pi11, pi00, pr, coeff) { #pi11, pi00: vector, pr: scalar
  norm = (coeff[2]/pr  + coeff[3] + 1) / (1 + coeff[1])
  return( norm * (coeff[1]*pr+pr*pi11) / (coeff[2] + coeff[3]*pr + coeff[4]*pr*pi11 + (1-pr)*(1-pi00)) )
}
## MCC_rb, with normalization
MCCgen.metric = function(pi11,pi00, pr, d) { #pi11, pi00: vector, pr: scalar
  norm = sqrt(d/pr/(1-pr) + 1)
  return( norm * sqrt(pr*(1-pr)) * ( pi11*pi00 - (1-pi11)*(1-pi00) ) / 
            sqrt( d + (pr*pi11 + (1-pr)*(1-pi00)) * (pr*(1-pi11)+(1-pr)*pi00) ) )
}
#
## old versions, without normalization
## F_gen-Score, standard: coeff[1]=c=0, coeff[2]=d_0!=0, coeff[3]=d_1=beta^2, coeff[4]=d_2=1,
#Fgen.metric = function(pi11, pi00, pr, coeff) { #pi11, pi00: vector, pr: scalar
#  return( (coeff[1]*pr+pr*pi11) / (coeff[2] + coeff[3]*pr + coeff[4]*pr*pi11 + (1-pr)*(1-pi00)) )
#}
## MCC_gen
#MCCgen.metric = function(pi11,pi00, pr, d) { #pi11, pi00: vector, pr: scalar
#  return( sqrt(pr*(1-pr)) * ( pi11*pi00 - (1-pi11)*(1-pi00) ) / 
#            sqrt( d + (pr*pi11 + (1-pr)*(1-pi00)) * (pr*(1-pi11)+(1-pr)*pi00) ) )
#}

## JAC
JAC.metric = function(pi11,pi00, pr) { #pi11, pi00: vector, pr: scalar
  return( (pr*pi11) / (1 - (1-pr)*pi00) )
}
## Yule
Yule.metric = function(pi11,pi00, pr) { #pi11, pi00: vector, pr: scalar
  return( (pi11*pi00 - (1-pi11)*(1-pi00)) / (pi11*pi00 + (1-pi11)*(1-pi00)) )
}
## Kappa
Kappa.metric = function(pi11,pi00, pr) { #pi11, pi00: vector, pr: scalar
  return( (2*pr*(1-pr) * (pi11+pi00-1)) / (2*pr*(1-pr) * (pi11+pi00-1) - pr*pi11 - (1-pr)*pi00 + 1) )
}
#############################
## end function definition ##
#############################

################################################################################
################################################################################
## Load packages
library(randomForest); library(caret); library(sampling); library(xtable)

# Importing the datasets
cs_training = read.csv("cs-training.csv")
#str(cs_training)
#head(cs_training)
#summary(cs_training)

# Data Cleaning
cs_training = cs_training[,(-1)] #remove names column

# split into training (80%) and test dataset, retain proportion of positives
set.seed(100)
n1 = sum(cs_training$SeriousDlqin2yrs==1)
n0 = sum(cs_training$SeriousDlqin2yrs==0)
ind = createDataPartition( as.factor(cs_training$SeriousDlqin2yrs), p=0.8)[[1]]
training = cs_training[ind,]
# sum(training$SeriousDlqin2yrs==1) / n1 # 0.8
summary(training)
test = cs_training[-ind,]


####################################################################################
# Fitting Logistic Regression to the Training set
# nodesize: Minimum size of terminal nodes. 
# default values: 1 for classification, 5 for regression (!)
# remove col1 (response SeriousDlqin2yrs), col6 (MonthlyIncome, many NA's), 
# col11 (NumberOfDependents, some NA's)
# treating NA's by multiple imputation does not improve model
set.seed(100)
LR.train = glm(SeriousDlqin2yrs ~ ., family=binomial, data=training[,-c(6,11)])
summary(LR.train)

# Making the Confusion Matrix
train.pred.num = predict(LR.train, type = 'response')
train.pred = factor( ifelse(train.pred.num >= 0.5, 1, 0))
cm.train = confusionMatrix(train.pred, as.factor(training$SeriousDlqin2yrs))
cm.train$table
###############################################################################
# computation of empirical optimal threshold for regression function
cg1 = train.pred.num
yf = factor(training$SeriousDlqin2yrs, levels = c(1,0))
del.vec = seq(0.001,0.999,0.001); nd = length(del.vec)
M.res = Fgen.res = MCCgen.res = matrix(0, nrow=nd, ncol=6)
n = length(training[,1])
n1 = length( training[training[,1]==1, 1] )
n0 = n-n1; pr = n1/n

start_time = Sys.time()
for (i in 1:nd) {
  lr.pred1 = factor( ifelse(cg1 >= del.vec[i], 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  cf1$table
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  #
  M.res[i, 1] = JAC.metric(pi1[1], pi1[2], pr)
  M.res[i, 2] = MCC.metric(pi1[1], pi1[2], pr)
  M.res[i, 3] = Yule.metric(pi1[1], pi1[2], pr)
  M.res[i, 4] = F.metric(pi1[1], pi1[2], pr, beta=1.5)
  M.res[i, 5] = F.metric(pi1[1], pi1[2], pr, beta=0.5)
  M.res[i, 6] = Kappa.metric(pi1[1], pi1[2], pr)
  #
  Fgen.res[i, 1] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.0,1,1))
  Fgen.res[i, 2] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.1,1,1))
  Fgen.res[i, 3] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.2,1,1))
  Fgen.res[i, 4] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.5,1,1))
  Fgen.res[i, 5] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.8,1,1))
  Fgen.res[i, 6] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.2,2.0,1))
  #
  MCCgen.res[i, 1] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.0)
  MCCgen.res[i, 2] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.01)
  MCCgen.res[i, 3] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.05)
  MCCgen.res[i, 4] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.1)
  MCCgen.res[i, 5] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.5)
  MCCgen.res[i, 6] = MCCgen.metric(pi1[1], pi1[2], pr, d=1.0)
  print(i)
}
end_time = Sys.time()
end_time - start_time

tab1 = tab2 = tab3 = matrix(0, nrow=4, ncol=6)
for (j in 1:6) {
  t1 = del.vec[ which.max( M.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( M.res[,j], na.rm=TRUE)
  tab1[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  t1 = del.vec[ which.max( Fgen.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( Fgen.res[,j], na.rm=TRUE)
  tab2[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  t1 = del.vec[ which.max( MCCgen.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( MCCgen.res[,j], na.rm=TRUE)
  tab3[1:4, j] =  c(M, t1, pi1)
}
#
cn1 = c("training", "$\\hat{\\pi}=0.067$", "", "")
cn2 = c("value", "$\\tilde\\delta$", "$\\hat\\pi_{1|1}$", "$\\hat\\pi_{0|0}$")
tab1.df = data.frame(cn1, cn2, tab1)
colnames(tab1.df) = c("","","JAC","MCC","Yule","$F_{1.5}$","$F_{0.5}$","Kappa")
print(xtable(tab1.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

tab2.df = data.frame(cn1, cn2, tab2)
colnames(tab2.df) = c("","$(d_0,d_1)$","$(0,1)$","$(0.1,1)$","$(0.2,1)$","$(0.5,1)$",
                      "$(0.8,1)$","$(0.2,2)$")
print(xtable(tab2.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

tab3.df = data.frame(cn1, cn2, tab3)
colnames(tab3.df) = c("","$d$","$0$","$0.01$","$0.05$","$0.1$","$0.5$","$1.0$")
print(xtable(tab3.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

####################################################################################
####################################################################################
# use subset of positives to increase imbalance
training1 = training
set.seed(2024)
training1.temp = training1[ training[,1]==1,]
training1.temp = training1.temp[sample(n1, round(n1/5)), ] #randomly select n1/5 of positives
training1 = rbind(training1.temp, training1[training1[,1]==0,])
dim(training1)
length(training1[training1$SeriousDlqin2yrs==1,1])

set.seed(100)
LR.train = glm(SeriousDlqin2yrs ~ ., family=binomial, data=training1[,-c(6,11)])

###############################################################################
# computation of empirical optimal threshold for regression function
cg1 = predict(LR.train, type = 'response')
yf = factor(training1$SeriousDlqin2yrs, levels = c(1,0))
del.vec = seq(0.001,0.999,0.001); nd = length(del.vec)
M.res = Fgen.res = MCCgen.res = matrix(0, nrow=nd, ncol=6)
n = length(training1[,1])
n1 = length( training1[training1[,1]==1, 1] )
n0 = n-n1; pr = n1/n

start_time = Sys.time()
for (i in 1:nd) {
  lr.pred1 = factor( ifelse(cg1 >= del.vec[i], 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  cf1$table
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  #
  M.res[i, 1] = JAC.metric(pi1[1], pi1[2], pr)
  M.res[i, 2] = MCC.metric(pi1[1], pi1[2], pr)
  M.res[i, 3] = Yule.metric(pi1[1], pi1[2], pr)
  M.res[i, 4] = F.metric(pi1[1], pi1[2], pr, beta=1.5)
  M.res[i, 5] = F.metric(pi1[1], pi1[2], pr, beta=0.5)
  M.res[i, 6] = Kappa.metric(pi1[1], pi1[2], pr)
  #
  Fgen.res[i, 1] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.0,1,1))
  Fgen.res[i, 2] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.1,1,1))
  Fgen.res[i, 3] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.2,1,1))
  Fgen.res[i, 4] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.5,1,1))
  Fgen.res[i, 5] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.8,1,1))
  Fgen.res[i, 6] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.2,2.0,1))
  #
  MCCgen.res[i, 1] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.0)
  MCCgen.res[i, 2] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.01)
  MCCgen.res[i, 3] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.05)
  MCCgen.res[i, 4] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.1)
  MCCgen.res[i, 5] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.5)
  MCCgen.res[i, 6] = MCCgen.metric(pi1[1], pi1[2], pr, d=1.0)
  print(i)
}
end_time = Sys.time()
end_time - start_time

tab1 = tab2 = tab3 = matrix(0, nrow=4, ncol=6)
for (j in 1:6) {
  t1 = del.vec[ which.max( M.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( M.res[,j], na.rm=TRUE)
  tab1[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  t1 = del.vec[ which.max( Fgen.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( Fgen.res[,j], na.rm=TRUE)
  tab2[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  t1 = del.vec[ which.max( MCCgen.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( MCCgen.res[,j], na.rm=TRUE)
  tab3[1:4, j] =  c(M, t1, pi1)
}
#
cn1 = c("subset.train", "$\\hat{\\pi}=0.014$", "", "")
cn2 = c("value", "$\\tilde\\delta$", "$\\hat\\pi_{1|1}$", "$\\hat\\pi_{0|0}$")
tab1.df = data.frame(cn1, cn2, tab1)
colnames(tab1.df) = c("","","JAC","MCC","Yule","$F_{1.5}$","$F_{0.5}$","Kappa")
print(xtable(tab1.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

tab2.df = data.frame(cn1, cn2, tab2)
colnames(tab2.df) = c("","$(d_0,d_1)$","$(0,1)$","$(0.1,1)$","$(0.2,1)$","$(0.5,1)$",
                      "$(0.8,1)$","$(0.2,2)$")
print(xtable(tab2.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

tab3.df = data.frame(cn1, cn2, tab3)
colnames(tab3.df) = c("","$d$","$0$","$0.01$","$0.05$","$0.1$","$0.5$","$1.0$")
print(xtable(tab3.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function


####################################################################################
####################################################################################
## repeat with test data ##
####################################################################################
####################################################################################

####################################################################################
# Fitting random forest to the test set
set.seed(100)
LR.test = glm(SeriousDlqin2yrs ~ ., family=binomial, data=test[,-c(6,11)])

###############################################################################
# computation of empirical optimal threshold for regression function
cg1 = predict(LR.test, type = 'response')
yf = factor(test$SeriousDlqin2yrs, levels = c(1,0))
del.vec = seq(0.001,0.999,0.001); nd = length(del.vec)
M.res = Fgen.res = MCCgen.res = matrix(0, nrow=nd, ncol=6)
n = length(test[,1])
n1 = length( test[test[,1]==1, 1] )
n0 = n-n1; pr = n1/n

start_time = Sys.time()
for (i in 1:nd) {
  lr.pred1 = factor( ifelse(cg1 >= del.vec[i], 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  cf1$table
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  #
  M.res[i, 1] = JAC.metric(pi1[1], pi1[2], pr)
  M.res[i, 2] = MCC.metric(pi1[1], pi1[2], pr)
  M.res[i, 3] = Yule.metric(pi1[1], pi1[2], pr)
  M.res[i, 4] = F.metric(pi1[1], pi1[2], pr, beta=1.5)
  M.res[i, 5] = F.metric(pi1[1], pi1[2], pr, beta=0.5)
  M.res[i, 6] = Kappa.metric(pi1[1], pi1[2], pr)
  #
  Fgen.res[i, 1] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.0,1,1))
  Fgen.res[i, 2] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.1,1,1))
  Fgen.res[i, 3] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.2,1,1))
  Fgen.res[i, 4] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.5,1,1))
  Fgen.res[i, 5] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.8,1,1))
  Fgen.res[i, 6] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.2,2.0,1))
  #
  MCCgen.res[i, 1] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.0)
  MCCgen.res[i, 2] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.01)
  MCCgen.res[i, 3] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.05)
  MCCgen.res[i, 4] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.1)
  MCCgen.res[i, 5] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.5)
  MCCgen.res[i, 6] = MCCgen.metric(pi1[1], pi1[2], pr, d=1.0)
  print(i)
}
end_time = Sys.time()
end_time - start_time

tab1 = tab2 = tab3 = matrix(0, nrow=4, ncol=6)
for (j in 1:6) {
  t1 = del.vec[ which.max( M.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( M.res[,j], na.rm=TRUE)
  tab1[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  t1 = del.vec[ which.max( Fgen.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( Fgen.res[,j], na.rm=TRUE)
  tab2[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  t1 = del.vec[ which.max( MCCgen.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( MCCgen.res[,j], na.rm=TRUE)
  tab3[1:4, j] =  c(M, t1, pi1)
}
#
cn1 = c("test", "$\\hat{\\pi}=0.067$", "", "")
cn2 = c("value", "$\\tilde\\delta$", "$\\hat\\pi_{1|1}$", "$\\hat\\pi_{0|0}$")
tab1.df = data.frame(cn1, cn2, tab1)
colnames(tab1.df) = c("","","JAC","MCC","Yule","$F_{1.5}$","$F_{0.5}$","Kappa")
print(xtable(tab1.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

tab2.df = data.frame(cn1, cn2, tab2)
colnames(tab2.df) = c("","$(d_0,d_1)$","$(0,1)$","$(0.1,1)$","$(0.2,1)$","$(0.5,1)$",
                      "$(0.8,1)$","$(0.2,2)$")
print(xtable(tab2.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

tab3.df = data.frame(cn1, cn2, tab3)
colnames(tab3.df) = c("","$d$","$0$","$0.01$","$0.05$","$0.1$","$0.5$","$1.0$")
print(xtable(tab3.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

####################################################################################
####################################################################################
# use subset of positives to increase imbalance
test1 = test
set.seed(2024)
test1.temp = test1[ test[,1]==1,]
test1.temp = test1.temp[sample(n1, round(n1/5)), ] #randomly select n1/5 of positives
test1 = rbind(test1.temp, test1[test1[,1]==0,])
length(test1[test1$SeriousDlqin2yrs==1,1])

set.seed(100)
LR.test = glm(SeriousDlqin2yrs ~ ., family=binomial, data=test1[,-c(6,11)])

###############################################################################
# computation of empirical optimal threshold for regression function
cg1 = predict(LR.test, type = 'response')
yf = factor(test1$SeriousDlqin2yrs, levels = c(1,0))
del.vec = seq(0.001,0.999,0.001); nd = length(del.vec)
M.res = Fgen.res = MCCgen.res = matrix(0, nrow=nd, ncol=6)
n = length(test1[,1])
n1 = length( test1[test1[,1]==1, 1] )
n0 = n-n1; pr = n1/n

start_time = Sys.time()
for (i in 1:nd) {
  lr.pred1 = factor( ifelse(cg1 >= del.vec[i], 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  cf1$table
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  #
  M.res[i, 1] = JAC.metric(pi1[1], pi1[2], pr)
  M.res[i, 2] = MCC.metric(pi1[1], pi1[2], pr)
  M.res[i, 3] = Yule.metric(pi1[1], pi1[2], pr)
  M.res[i, 4] = F.metric(pi1[1], pi1[2], pr, beta=1.5)
  M.res[i, 5] = F.metric(pi1[1], pi1[2], pr, beta=0.5)
  M.res[i, 6] = Kappa.metric(pi1[1], pi1[2], pr)
  #
  Fgen.res[i, 1] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.0,1,1))
  Fgen.res[i, 2] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.1,1,1))
  Fgen.res[i, 3] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.2,1,1))
  Fgen.res[i, 4] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.5,1,1))
  Fgen.res[i, 5] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.8,1,1))
  Fgen.res[i, 6] = Fgen.metric(pi1[1], pi1[2], pr, coeff=c(0,0.2,2.0,1))
  #
  MCCgen.res[i, 1] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.0)
  MCCgen.res[i, 2] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.01)
  MCCgen.res[i, 3] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.05)
  MCCgen.res[i, 4] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.1)
  MCCgen.res[i, 5] = MCCgen.metric(pi1[1], pi1[2], pr, d=0.5)
  MCCgen.res[i, 6] = MCCgen.metric(pi1[1], pi1[2], pr, d=1.0)
  print(i)
}
end_time = Sys.time()
end_time - start_time

tab1 = tab2 = tab3 = matrix(0, nrow=4, ncol=6)
for (j in 1:6) {
  t1 = del.vec[ which.max( M.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( M.res[,j], na.rm=TRUE)
  tab1[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  t1 = del.vec[ which.max( Fgen.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( Fgen.res[,j], na.rm=TRUE)
  tab2[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  t1 = del.vec[ which.max( MCCgen.res[,j] ) ]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = max( MCCgen.res[,j], na.rm=TRUE)
  tab3[1:4, j] =  c(M, t1, pi1)
}
#
cn1 = c("subset.test", "$\\hat{\\pi}=0.014$", "", "")
cn2 = c("value", "$\\tilde\\delta$", "$\\hat\\pi_{1|1}$", "$\\hat\\pi_{0|0}$")
tab1.df = data.frame(cn1, cn2, tab1)
colnames(tab1.df) = c("","","JAC","MCC","Yule","$F_{1.5}$","$F_{0.5}$","Kappa")
print(xtable(tab1.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

tab2.df = data.frame(cn1, cn2, tab2)
colnames(tab2.df) = c("","$(d_0,d_1)$","$(0,1)$","$(0.1,1)$","$(0.2,1)$","$(0.5,1)$",
                      "$(0.8,1)$","$(0.2,2)$")
print(xtable(tab2.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

tab3.df = data.frame(cn1, cn2, tab3)
colnames(tab3.df) = c("","$d$","$0$","$0.01$","$0.05$","$0.1$","$0.5$","$1.0$")
print(xtable(tab3.df, digits=3), include.rownames=FALSE, include.colnames=TRUE,
      sanitize.text.function=function(x){x}) #use sanitization function

