########################################
## Plot of ROC curves in data example ##
########################################
## using gradient boosting instead of
## random forest
########################################

par(mar=c(4,4,0,2), mfcol=c(1,1), pty="s")

#######################################
## complete test data
#######################################
pr = 0.067 # proportion in full test set
cond.prop = pi.test.gbf
tpr.fpr = cbind(1-cond.prop[,2], cond.prop[,1])

plot(tpr.fpr, type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", lwd=3,
     xlab="FPR / 1-precision", ylab="TPR")
lines(tpr.fpr, lwd=3, col="black")      # ROC curve
abline(0,1,lty=2)

# recall against 1-precision
prec.rec = cbind(cond.prop[,1],
                 cond.prop[,1]*pr / (cond.prop[,1]*pr + (1-cond.prop[,2])*(1-pr)))
lines(1-prec.rec[,2], prec.rec[,1], lwd=3, lty=3, col="black")

# MCC
cond.prop1 = tab1.gbf[3:4,2]
tpr.fpr1 = cbind(1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=16, cex=1.5)
prec.rec1 = cbind(cond.prop1[1],
                  cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)))
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=16, cex=1.5)

# MCC(0.01)
cond.prop1 = tab3.gbf[3:4,2]
tpr.fpr1 = cbind(1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=17, cex=1.5)
prec.rec1 = cbind(cond.prop1[1],
                  cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)))
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=17, cex=1.5)

# MCC(0.05)
cond.prop1 = tab3.gbf[3:4,3]
tpr.fpr1 = cbind(1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=18, cex=1.5)
prec.rec1 = cbind(cond.prop1[1],
                  cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)))
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=18, cex=1.5)

# MCC(0.1)
cond.prop1 = tab3.gbf[3:4,4]
tpr.fpr1 = cbind(1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=15, cex=1.5)
prec.rec1 = cbind(cond.prop1[1],
                  cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)))
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=15, cex=1.5)

#######################################
## subset of test data
#######################################
pr = 0.014 # proportion in subset of test set
cond.prop = pi.test.gbs
tpr.fpr = cbind(1-cond.prop[,2], cond.prop[,1])

lines(tpr.fpr, lwd=3, col="darkgrey")   # ROC curve

# recall against 1-precision
prec.rec = cbind(cond.prop[,1],
                 cond.prop[,1]*pr / (cond.prop[,1]*pr + (1-cond.prop[,2])*(1-pr)))
lines(1-prec.rec[,2], prec.rec[,1], lwd=3, lty=3, col="darkgrey")

# MCC
cond.prop1 = tab1.gbs[3:4,2]
tpr.fpr1 = cbind(1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="darkgrey", pch=16, cex=1.5)
prec.rec1 = cbind(cond.prop1[1],
                  cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)))
points(1-prec.rec1[2], prec.rec1[1], col="darkgrey", pch=16, cex=1.5)

# MCC(0.01)
cond.prop1 = tab3.gbs[3:4,2]
tpr.fpr1 = cbind(1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="darkgrey", pch=17, cex=1.5)
prec.rec1 = cbind(cond.prop1[1],
                  cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)))
points(1-prec.rec1[2], prec.rec1[1], col="darkgrey", pch=17, cex=1.5)

# MCC(0.05)
cond.prop1 = tab3.gbs[3:4,3]
tpr.fpr1 = cbind(1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="darkgrey", pch=18, cex=1.5)
prec.rec1 = cbind(cond.prop1[1],
                  cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)))
points(1-prec.rec1[2], prec.rec1[1], col="darkgrey", pch=18, cex=1.5)

# MCC(0.1)
cond.prop1 = tab3.gbs[3:4,4]
tpr.fpr1 = cbind(1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="darkgrey", pch=15, cex=1.5)
prec.rec1 = cbind(cond.prop1[1],
                  cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)))
points(1-prec.rec1[2], prec.rec1[1], col="darkgrey", pch=15, cex=1.5)

legend("bottomright",
       legend=c("ROC, π = 0.067", "ROC, π = 0.014",
                "RP, π = 0.067",  "RP, π = 0.014"),
       col=c("black","darkgrey","black","darkgrey"),
       lty=c(1,1,3,3), lwd=3, cex=0.8)