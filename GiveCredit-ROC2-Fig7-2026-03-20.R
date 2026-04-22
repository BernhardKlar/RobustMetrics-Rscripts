########################################
## Plot of ROC curves in data example ##
########################################
## Left: complete data, right: subset ##
########################################

#######################################
## Plot ROC curves for complete data ##
#######################################
par(mar=c(3,4,0,2), mfcol=c(1,2), pty="s")

cond.prop = pi.test.lrf 
tpr.fpr = cbind( 1-cond.prop[,2], cond.prop[,1])
plot( tpr.fpr, type="n", xlim=c(0,1), ylim=c(0,1), xaxs = "i", yaxs = "i", lwd=3,
      xlab="FPR / 1-precision", ylab="TPR")             
lines( tpr.fpr, lwd=3)  #ROC curve, i.e. TPR(recall)-FPR-curve
abline(0,1,lty=2)
pr = 0.067 # proportion in full test set
# precision-recall-curve
prec.rec = cbind(cond.prop[,1], cond.prop[,1]*pr / (cond.prop[,1]*pr + (1-cond.prop[,2])*(1-pr)) ) 
lines(1-prec.rec[,2], prec.rec[,1], lwd=3, lty=3, col="black") #recall--1-precision-curve
# MCC
cond.prop1 = tab1.lrf[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=16, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=16, cex=1.5)
# MCC(0.01)
cond.prop1 = tab3.lrf[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=17, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=17, cex=1.5)
# MCC(0.05)
cond.prop1 = tab3.lrf[3:4,3]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=18, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=18, cex=1.5)
# MCC(0.1)
cond.prop1 = tab3.lrf[3:4,4]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=15, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=15, cex=1.5)
############################################
cond.prop = pi.test.full 
tpr.fpr = cbind( 1-cond.prop[,2], cond.prop[,1])
col_rf    <- "grey40"
lines( tpr.fpr, lwd=3, col=col_rf)  #ROC curve, i.e. TPR(recall)-FPR-curve
pr = 0.067 # proportion in full test set
# precision-recall-curve
prec.rec = cbind(cond.prop[,1], cond.prop[,1]*pr / (cond.prop[,1]*pr + (1-cond.prop[,2])*(1-pr)) ) 
lines(1-prec.rec[,2], prec.rec[,1], lwd=3, lty=3, col=col_rf) #recall--1-precision-curve
# MCC
cond.prop1 = tab1.f[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_rf, pch=16, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_rf, pch=16, cex=1.5)
# MCC(0.01)
cond.prop1 = tab3.f[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_rf, pch=17, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_rf, pch=17, cex=1.5)
# MCC(0.05)
cond.prop1 = tab3.f[3:4,3]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_rf, pch=18, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_rf, pch=18, cex=1.5)
# MCC(0.1)
cond.prop1 = tab3.f[3:4,4]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_rf, pch=15, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_rf, pch=15, cex=1.5)
############################################
cond.prop = pi.test.gbf 
tpr.fpr = cbind( 1-cond.prop[,2], cond.prop[,1])
col_gb = "grey60"
lines( tpr.fpr, lwd=3, col=col_gb)  #ROC curve, i.e. TPR(recall)-FPR-curve
pr = 0.067 # proportion in full test set
# precision-recall-curve
prec.rec = cbind(cond.prop[,1], cond.prop[,1]*pr / (cond.prop[,1]*pr + (1-cond.prop[,2])*(1-pr)) ) 
lines(1-prec.rec[,2], prec.rec[,1], lwd=3, lty=3, col=col_gb) #recall--1-precision-curve
# MCC
cond.prop1 = tab1.gbf[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_gb, pch=16, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_gb, pch=16, cex=1.5)
# MCC(0.01)
cond.prop1 = tab3.gbf[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_gb, pch=17, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_gb, pch=17, cex=1.5)
# MCC(0.05)
cond.prop1 = tab3.gbf[3:4,3]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_gb, pch=18, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_gb, pch=18, cex=1.5)
# MCC(0.1)
cond.prop1 = tab3.gbf[3:4,4]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_gb, pch=15, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_gb, pch=15, cex=1.5)
############################################
legend("bottomright", legend=c("ROC, logistic regression ", "ROC, random forest ", 
                           "ROC, gradient boosting ", "RP, logistic regression ",  
                           "RP, random forest ", "RP, gradient boosting "),
       col=rep( c("black",col_rf, col_gb), 2), lty=c(1,1,1,3,3,3), lwd=3, cex=0.7)

################################
## Plot ROC curves for subset ##
################################
pr = 0.014 # proportion in subset of test set
cond.prop = pi.test.lrs 
tpr.fpr = cbind( 1-cond.prop[,2], cond.prop[,1])
plot( tpr.fpr, type="n", xlim=c(0,1), ylim=c(0,1), xaxs = "i", yaxs = "i", lwd=3,
      xlab="FPR / 1-precision", ylab="TPR")             
lines( tpr.fpr, lwd=3, col="black")  #ROC curve, i.e. TPR(recall)-FPR-curve
abline(0,1,lty=2)
# precision-recall-curve
prec.rec = cbind(cond.prop[,1], cond.prop[,1]*pr / (cond.prop[,1]*pr + (1-cond.prop[,2])*(1-pr)) ) 
lines(1-prec.rec[,2], prec.rec[,1], lwd=3, lty=3, col="black") #recall--1-precision-curve
# MCC
cond.prop1 = tab1.lrs[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=16, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=16, cex=1.5)
# MCC(0.01)
cond.prop1 = tab3.lrs[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=17, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=17, cex=1.5)
# MCC(0.05)
cond.prop1 = tab3.lrs[3:4,3]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col="black", pch=18, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col="black", pch=18, cex=1.5)
# MCC(0.1)
cond.prop1 = tab3.lrs[3:4,4]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1 + c(0.02,0), col="black", pch=15, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2]+0.02, prec.rec1[1], col="black", pch=15, cex=1.5)
############################################
cond.prop = pi.test.sub 
tpr.fpr = cbind( 1-cond.prop[,2], cond.prop[,1])
lines( tpr.fpr, lwd=3, col=col_rf)  #ROC curve, i.e. TPR(recall)-FPR-curve
# precision-recall-curve
prec.rec = cbind(cond.prop[,1], cond.prop[,1]*pr / (cond.prop[,1]*pr + (1-cond.prop[,2])*(1-pr)) ) 
lines(1-prec.rec[,2], prec.rec[,1], lwd=3, lty=3, col=col_rf) #recall--1-precision-curve
# MCC
cond.prop1 = tab1.s[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_rf, pch=16, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_rf, pch=16, cex=1.5)
# MCC(0.01)
cond.prop1 = tab3.s[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_rf, pch=17, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_rf, pch=17, cex=1.5)
# MCC(0.05)
cond.prop1 = tab3.s[3:4,3]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_rf, pch=18, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_rf, pch=18, cex=1.5)
# MCC(0.1)
cond.prop1 = tab3.s[3:4,4]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1 + c(0.02,0), col=col_rf, pch=15, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2]+0.02, prec.rec1[1], col=col_rf, pch=15, cex=1.5)
############################################
cond.prop = pi.test.gbs 
tpr.fpr = cbind( 1-cond.prop[,2], cond.prop[,1])
lines( tpr.fpr, lwd=3, col=col_gb)  #ROC curve, i.e. TPR(recall)-FPR-curve
# precision-recall-curve
prec.rec = cbind(cond.prop[,1], cond.prop[,1]*pr / (cond.prop[,1]*pr + (1-cond.prop[,2])*(1-pr)) ) 
lines(1-prec.rec[,2], prec.rec[,1], lwd=3, lty=3, col=col_gb) #recall--1-precision-curve
# MCC
cond.prop1 = tab1.gbs[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_gb, pch=16, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_gb, pch=16, cex=1.5)
# MCC(0.01)
cond.prop1 = tab3.gbs[3:4,2]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_gb, pch=17, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_gb, pch=17, cex=1.5)
# MCC(0.05)
cond.prop1 = tab3.gbs[3:4,3]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_gb, pch=18, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_gb, pch=18, cex=1.5)
# MCC(0.1)
cond.prop1 = tab3.gbs[3:4,4]
tpr.fpr1 = cbind( 1-cond.prop1[2], cond.prop1[1])
points(tpr.fpr1, col=col_gb, pch=15, cex=1.5)
prec.rec1 = cbind(cond.prop1[1], cond.prop1[1]*pr / (cond.prop1[1]*pr + (1-cond.prop1[2])*(1-pr)) ) 
points(1-prec.rec1[2], prec.rec1[1], col=col_gb, pch=15, cex=1.5)
############################################

legend("bottom", legend=c("ROC, logistic regression ", "ROC, random forest ", 
                           "ROC, gradient boosting ", "RP, logistic regression ",  
                           "RP, random forest ", "RP, gradient boosting "),
       col=rep( c("black",col_rf, col_gb), 2), lty=c(1,1,1,3,3,3), lwd=3, cex=0.7)

