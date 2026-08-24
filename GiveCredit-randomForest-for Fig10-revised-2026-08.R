################################################
## Data example in Section 7, , Random Forest ##
################################################

################################################################
## Values for plots saved in
##
## Logistic: tab1.lrf, tab2.lrf, tab3.lrf, pi.test.lrf  (full)
##           tab1.lrs, tab2.lrs, tab3.lrs, pi.test.lrs  (subset)
##
## RF:       tab1.f, tab2.f, tab3.f, pi.test.full (full)
##           tab1.s, tab2.s, tab3.s, pi.test.sub (subset)
##
## GB:       tab1.gbf, tab2.gbf, tab3.gbf, pi.test.gbf (full)
##           tab1.gbs, tab2.gbs, tab3.gbs, pi.test.gbs (subset)
#################################################################


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

# Importing dataset
cs_training = read.csv("cs-training.csv")
str(cs_training)
dim(cs_training)
summary(cs_training)

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
# Fitting random forest to the Training set
# nodesize: Minimum size of terminal nodes. 
# default values: 1 for classification, 5 for regression (!)
# remove col1 (response SeriousDlqin2yrs), col6 (MonthlyIncome, many NA's), 
# col11 (NumberOfDependents, some NA's)
# treating NA's by multiple imputation does not improve model
set.seed(100)
time.fit.full = system.time({
  RF.train = randomForest(training[,-c(1,6,11)], training$SeriousDlqin2yrs,
                sampsize=c(10000), do.trace=FALSE, importance=TRUE, ntree=1000, keep.forest=TRUE)
})["elapsed"]

# Making the Confusion Matrix using all-tree in-sample predictions
## change begin ## 
#time.pred.train.full = system.time({
#  cg.train = predict(RF.train, newdata=training[,-c(1,6,11)])
#})["elapsed"]
#train.pred = factor( ifelse(cg.train >= 0.5, 1, 0))

train.pred = factor( ifelse(RF.train$predicted >= 0.5, 1, 0))
## change end ## 
cm.train = confusionMatrix(train.pred, as.factor(training$SeriousDlqin2yrs))
cm.train$table
###############################################################################
# computation of empirical optimal threshold for regression function
cg1 = cg.train
yf = factor(training$SeriousDlqin2yrs, levels = c(1,0))
del.vec = seq(0.001,0.999,0.001); nd = length(del.vec)
M.res = Fgen.res = MCCgen.res = matrix(0, nrow=nd, ncol=6)
n = length(training[,1])
n1 = length( training[training[,1]==1, 1] )
n0 = n-n1; pr = n1/n

time.threshold.full = system.time({
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
  # print(i)
}
})["elapsed"]
time.threshold.full

tab1 = tab2 = tab3 = matrix(0, nrow=4, ncol=6)
ind.M.full = ind.Fgen.full = ind.MCCgen.full = integer(6)
for (j in 1:6) {
  ind.M.full[j] = which.max( M.res[,j] )
  i1 = ind.M.full[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = M.res[i1,j]
  tab1[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  ind.Fgen.full[j] = which.max( Fgen.res[,j] )
  i1 = ind.Fgen.full[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = Fgen.res[i1,j]
  tab2[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  ind.MCCgen.full[j] = which.max( MCCgen.res[,j] )
  i1 = ind.MCCgen.full[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = MCCgen.res[i1,j]
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
time.fit.sub = system.time({
  RF.train1 = randomForest(training1[,-c(1,6,11)], training1$SeriousDlqin2yrs,
                          sampsize=c(10000), do.trace=FALSE, importance=TRUE, ntree=1000, keep.forest=TRUE)
})["elapsed"]

###############################################################################
# computation of empirical optimal threshold for regression function
## change begin ## 
#time.pred.train.sub = system.time({
#  cg1 = predict(RF.train1, newdata=training1[,-c(1,6,11)])
#})["elapsed"]
cg1 = RF.train1$predicted
## change end ## 
yf = factor(training1$SeriousDlqin2yrs, levels = c(1,0))
del.vec = seq(0.001,0.999,0.001); nd = length(del.vec)
M.res = Fgen.res = MCCgen.res = matrix(0, nrow=nd, ncol=6)
n = length(training1[,1])
n1 = length( training1[training1[,1]==1, 1] )
n0 = n-n1; pr = n1/n

time.threshold.sub = system.time({
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
  # print(i)
}
})["elapsed"]
time.threshold.sub

tab1 = tab2 = tab3 = matrix(0, nrow=4, ncol=6)
ind.M.sub = ind.Fgen.sub = ind.MCCgen.sub = integer(6)
for (j in 1:6) {
  ind.M.sub[j] = which.max( M.res[,j] )
  i1 = ind.M.sub[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = M.res[i1,j]
  tab1[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  ind.Fgen.sub[j] = which.max( Fgen.res[,j] )
  i1 = ind.Fgen.sub[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = Fgen.res[i1,j]
  tab2[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  ind.MCCgen.sub[j] = which.max( MCCgen.res[,j] )
  i1 = ind.MCCgen.sub[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = MCCgen.res[i1,j]
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
# Prediction on the untouched test set using the forest fitted to the training set
time.pred.test.full = system.time({
  cg1 = predict(RF.train, newdata=test[,-c(1,6,11)])
})["elapsed"]

###############################################################################
# evaluation on the test set at thresholds selected on the training set
yf = factor(test$SeriousDlqin2yrs, levels = c(1,0))
del.vec = seq(0.001,0.999,0.001); nd = length(del.vec)
M.res = Fgen.res = MCCgen.res = matrix(0, nrow=nd, ncol=6)
pi.test.full = matrix(0, nrow=nd, ncol=2)
n = length(test[,1])
n1 = length( test[test[,1]==1, 1] )
n0 = n-n1; pr = n1/n

time.eval.test.full = system.time({
for (i in 1:nd) {
  lr.pred1 = factor( ifelse(cg1 >= del.vec[i], 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  cf1$table
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  pi.test.full[i,] = pi1 #save result for ROC curve 
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
  # print(i)
}
})["elapsed"]
time.eval.test.full

tab1 = tab2 = tab3 = matrix(0, nrow=4, ncol=6)
for (j in 1:6) {
  i1 = ind.M.full[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = M.res[i1,j]
  tab1[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  i1 = ind.Fgen.full[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = Fgen.res[i1,j]
  tab2[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  i1 = ind.MCCgen.full[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = MCCgen.res[i1,j]
  tab3[1:4, j] =  c(M, t1, pi1)
}
#
tab1.f = tab1; tab2.f = tab2; tab3.f = tab3 #save results for ROC curve
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

time.pred.test.sub = system.time({
  cg1 = predict(RF.train1, newdata=test1[,-c(1,6,11)])
})["elapsed"]

###############################################################################
# evaluation on the test set at thresholds selected on the training set
yf = factor(test1$SeriousDlqin2yrs, levels = c(1,0))
del.vec = seq(0.001,0.999,0.001); nd = length(del.vec)
M.res = Fgen.res = MCCgen.res = matrix(0, nrow=nd, ncol=6)
pi.test.sub = matrix(0, nrow=nd, ncol=2)
n = length(test1[,1])
n1 = length( test1[test1[,1]==1, 1] )
n0 = n-n1; pr = n1/n

time.eval.test.sub = system.time({
for (i in 1:nd) {
  lr.pred1 = factor( ifelse(cg1 >= del.vec[i], 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  cf1$table
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  pi.test.sub[i,] = pi1 #save result for ROC curve 
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
  # print(i)
}
})["elapsed"]
time.eval.test.sub

## evaluate on test data for comparison
#d.values <- c(0, 0.01, 0.05, 0.1, 0.5, 1)
#ind.test.opt <- vapply(seq_along(d.values), 
#                       function(j) which.max(MCCgen.res[, j]), integer(1))
#test.opt <- data.frame(d = d.values, threshold = del.vec[ind.test.opt],
#  TPR = pi.test.sub[ind.test.opt, 1], TNR = pi.test.sub[ind.test.opt, 2],
#  value = MCCgen.res[cbind(ind.test.opt, seq_along(d.values))])
#test.opt
#

tab1 = tab2 = tab3 = matrix(0, nrow=4, ncol=6)
for (j in 1:6) {
  i1 = ind.M.sub[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = M.res[i1,j]
  tab1[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  i1 = ind.Fgen.sub[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = Fgen.res[i1,j]
  tab2[1:4, j] =  c(M, t1, pi1)
}
#
for (j in 1:6) {
  i1 = ind.MCCgen.sub[j]
  t1 = del.vec[i1]
  lr.pred1 = factor( ifelse(cg1>=t1, 1, 0), levels=c(1,0))
  cf1 = confusionMatrix(lr.pred1, yf)
  pi1 = c(cf1$table[1,1]/n1, cf1$table[2,2]/n0)
  M = MCCgen.res[i1,j]
  tab3[1:4, j] =  c(M, t1, pi1)
}
#
tab1.s = tab1; tab2.s = tab2; tab3.s = tab3 #save results for ROC curve
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
##################################################################################

####################################################################################
# Computational timings (elapsed seconds)
timings = data.frame(
  setting = c("full data", "strongly imbalanced subset"),
  model = "regression forest",
  fit = c(time.fit.full, time.fit.sub),
  train_prediction = c(time.pred.train.full, time.pred.train.sub),
  threshold_search_train = c(time.threshold.full, time.threshold.sub),
  test_prediction = c(time.pred.test.full, time.pred.test.sub),
  test_evaluation = c(time.eval.test.full, time.eval.test.sub)
)
print(timings, row.names=FALSE)

