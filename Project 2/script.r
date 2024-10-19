A = read.csv("tabella.csv", row.names = 1)
B = A[,-1]
B[,2] = as.numeric(B[,2])

# Comandi per generare la curva ROC presi dal file s2_cmroc.r usato a lezione
s2_roc <- function(response,predictor){
  num=1000
  rpredictor<-round(predictor,log10(num))
  roc=matrix(rep(3*(num+1),0),num+1,3)
  for(i in 0:num){
    p=i/num
    TrueNegative=sum((rpredictor<p)&(response==0))
    TruePositive=sum((rpredictor>p)&(response==1))
    FalsePositive=sum((rpredictor>p)&(response==0))
    FalseNegative=sum((rpredictor<p)&(response==1))
    sens=if(TruePositive==0)0 else TruePositive/(TruePositive+FalseNegative)
    spec=if(TrueNegative==0)0 else TrueNegative/(TrueNegative+FalsePositive)
    roc[i+1,]=c(spec,sens,p)
  }
  roc<-data.frame(roc)
  colnames(roc)<-c("Specificity ","Sensitivity ","Soglia")
  roc
}

s2_roc.plot<-function(roc,...){
  plot(roc[,1:2], type="l",main="ROC curve",
       ylab="Sensitivity",xlab="Specificity",
       xlim=c(1,0),ylim=c(0,1),...)
  segments(0,1,1,0,col="red",lwd=2,lty="dashed")
}

s2_roc.lines<-function(roc,...){
  lines(roc[,1:2], type="l",main="ROC curve",
        ylab="Sensitivity",xlab="Specificity",
        xlim=c(1,0),ylim=c(0,1),...)
}

# Comandi per generare la matrice di confusione presi dal file s2_cmroc.r usato a lezione
s2_confusion <- function(response,predictor,p=0.5){
  cm=matrix(nrow=2,ncol=2)
  
  cm[1,1]=sum((predictor>=p)&(response==1)) # true positive
  cm[1,2]=sum((predictor>=p)&(response==0)) # false positive
  cm[2,1]=sum((predictor<=p)&(response==1)) # false negative
  cm[2,2]=sum((predictor<=p)&(response==0)) # true negative
  
  cm<-data.frame(cm)
  rownames(cm)<-c("predicted 1","predicted 0")
  colnames(cm)<-c("actual 1","actual 0")
  cm
}

# Curva ROC regressione logistica
B.glm = glm(Classe~., data = B, family = binomial)
B.glm.p = predict(B.glm, type = "response")
B.glm.roc = s2_roc(B$Classe, B.glm.p)
s2_roc.plot(B.glm.roc, col = "green3")

# Curva ROC analisi discriminante lineare
library(MASS)
B.lda = lda(Classe~., data = B, CV = F)
B.lda.p = predict(B.lda)
B.lda.post = B.lda.p$posterior[,2]
B.lda.roc = s2_roc(B$Classe, B.lda.post)
s2_roc.lines(B.lda.roc, col = "red")

# Curva ROC analisi discriminante quadratica
B.qda = qda(Classe~., data = B, CV = F)
B.qda.p = predict(B.qda)
B.qda.post = B.qda.p$posterior[,2]
B.qda.roc = s2_roc(B$Classe, B.qda.post)
s2_roc.lines(B.qda.roc, col = "blue")

legend("bottomright",legend = c("Log", "lda", "qda"), col = c("green3", "red", "blue"), lty = c(1,1,1))

# Analisi regressione logistica
library(pROC)
plot.roc(roc(B$Classe, B.glm.p), print.auc = TRUE, print.thres = "best")
s2_confusion(B$Classe, B.glm.p, 0.734)
sum((B.glm.p > 0.734) == (B$Classe > 0.734))/length(B$Classe)

library(ggplot2)
library(WVPlots)
b.glm <- rbind(data.frame(x =  B.glm.p , y = B$Classe))
ThresholdPlot(b.glm, 'x' , 'y' , title = "ROC 'unrolled'", truth_target = TRUE,
              metrics = c("false_negative_rate"), linecolor = "green3") +
  geom_vline(xintercept = 0.49, color = "#d95f02")
ThresholdPlot(b.glm, 'x' , 'y' , title = "ROC 'unrolled'", truth_target = TRUE,
              metrics = c("false_positive_rate"), linecolor = "green3") +
  geom_vline(xintercept = 0.49, color = "#d95f02")

s2_confusion(B$Classe, B.glm.p, 0.49)
sum((B.glm.p > 0.49) == (B$Classe > 0.49))/length(B$Classe)

l = length(B$Classe)
acc.glm = rep(0,75)
for(i in 1:75){
  idx = sample(l,75)
  Bcv = B[-idx,]
  Bcv.glm = glm(Classe~., family = binomial, data = Bcv)
  Bcv.glm.p = predict(Bcv.glm, B[idx,], type="response")
  acc.glm[i] = sum((Bcv.glm.p > 0.49) == (B$Classe[idx] > 0.49))/75
}
mean(acc.glm)
sd(acc.glm)
hist(acc.glm, col = "green3", xlab = "Accuratezza", ylab = "Frequenza")

idx = sample(500,500)
acc1 = rep(0,500)
for(i in 1:500){
  bnf = B
  bnf$Classe[idx[1:i]] = -bnf$Classe[idx[1:i]]+rep(1,i)
  bnf.glm = glm(Classe~., data=bnf, family = binomial)
  bnf.glm.p = predict(bnf.glm, type="response")
  acc1[i] = sum((bnf.glm.p > 0.49) == (B$Classe > 0.49))/l
}
plot(acc1,type="l", xlab = "Errori", ylab = "Accuratezza")
segments(-50,0.8,600,0.8,lty=3,col="red")

# Analisi lda 
plot.roc(roc(B$Classe, B.lda.post), print.auc = TRUE, print.thres = "best")
s2_confusion(B$Classe, B.lda.post, 0.617)
sum((B.lda.post > 0.617) == (B$Classe > 0.617))/length(B$Classe)

b.lda <- rbind(data.frame(x =  B.lda.post , y = B$Classe))
ThresholdPlot(b.lda, 'x' , 'y' , title = "ROC 'unrolled'", truth_target = TRUE,
              metrics = c("false_negative_rate", "false_positive_rate"), linecolor = "red") +
  geom_vline(xintercept = 0.49, color = "#d95f02")

s2_confusion(B$Classe, B.lda.post, 0.49)
sum((B.lda.post > 0.49) == (B$Classe > 0.49))/length(B$Classe)

acc.lda = rep(0,75)
for(i in 1:75){
  idx = sample(l,75)
  Bcv = B[-idx,]
  Bcv.lda = lda(Classe~., data = Bcv)
  Bcv.lda.p = predict(Bcv.lda, B[idx,])$posterior[,2]
  acc.lda[i] = sum((Bcv.lda.p > 0.49) == (B$Classe[idx] > 0.49))/75
}
mean(acc.lda)
sd(acc.lda)
hist(acc.lda, col = "green3")

idx = sample(500,500)
acc2 = rep(0,500)
for(i in 1:500){
  bnf = B
  bnf$Classe[idx[1:i]] = -bnf$Classe[idx[1:i]]+rep(1,i)
  bnf.lda = lda(Classe~., data=bnf)
  bnf.lda.p = predict(bnf.lda)$posterior[,2]
  acc2[i] = sum((bnf.lda.p > 0.49) == (B$Classe > 0.49))/l
}
plot(acc2,type="l")
segments(-50,0.8,600,0.8,lty=3,col="red")

# Analisi qda
plot.roc(roc(B$Classe, B.qda.post), print.auc = TRUE, print.thres = "best")
s2_confusion(B$Classe, B.qda.post, 0.7)
sum((B.qda.post > 0.426) == (B$Classe > 0.426))/length(B$Classe)

library(ROCR)
pred <- prediction(B.qda.post, B$Classe)
perf<-ROCR::performance(pred,"prec","rec")
plot(perf,
     avg= "threshold",
     colorize=TRUE,
     lwd= 3,
     xlab="Sensitivity",
     ylab = "Precision"
     )

acc.qda = rep(0,75)
for(i in 1:75){
  idx = sample(l,75)
  Bcv = B[-idx,]
  Bcv.qda = qda(Classe~., data = Bcv)
  Bcv.qda.p = predict(Bcv.qda, B[idx,])$posterior[,2]
  acc.qda[i] = sum((Bcv.qda.p > 0.426) == (B$Classe[idx] > 0.426))/75
}
mean(acc.qda)
sd(acc.qda)
hist(acc.qda, col = "green3", xlab = "Accuratezza", ylab = "Frequenza")

idx = sample(500,500)
acc3 = rep(0,500)
for(i in 1:500){
  bnf = B
  bnf$Classe[idx[1:i]] = -bnf$Classe[idx[1:i]]+rep(1,i)
  bnf.qda = qda(Classe~., data=bnf)
  bnf.qda.p = predict(bnf.qda)$posterior[,2]
  acc3[i] = sum((bnf.qda.p > 0.426) == (B$Classe > 0.426))/l
}
plot(acc3,type="l", xlab= "Errori", ylab="Accuratezza")
segments(-50,0.7,600,0.7,lty=3,col="red")
segments(-50,0.8,600,0.8,lty=3,col="red")