A = read.csv("tabella.csv", row.names = 1)
B = A[,-13]
colnames(B) = rep(1:12)
cor(B)
library(corrplot)
corrplot(cor(B), method = "circle", tl.col = "black")
pca = princomp(scale(B))
summary(pca)
library(factoextra)
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 35))
plot(cumsum(pca$sdev^2)/sum(pca$sdev^2), type = "b", ylim = c(0.2,1), xlab = "Dimensioni", ylab = "Prop. varianza cumulata" )
segments(0,0.7,14,0.7,lty=3,col="red")
segments(0,0.8,14,0.8,lty=3,col="red")
plot(pca, col = "green3")
segments(0,1,12,1,lty=3,col="red")
loadings(pca)
ld = round(loadings(pca)[,1:4],3)
ld
corrplot(t(ld), method = "circle", tl.col = "black")
varimax(ld)
corrplot(t(varimax(ld)$loadings), method = "circle", tl.col = "black")
fviz_pca_var(pca, axes=c(1,2), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), xlim=c(-1,1), ylim=c(-1,1))
fviz_pca_var(pca, axes=c(1,3), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), xlim=c(-1,1), ylim=c(-1,1))
fviz_pca_var(pca, axes=c(1,4), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), xlim=c(-1,1), ylim=c(-1,1))
fviz_pca_var(pca, axes=c(2,3), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), xlim=c(-1,1), ylim=c(-1,1))
fviz_pca_ind(pca, axes = c(1,2), col.ind = A$EU, palette = c( "orange", "green3", "red", "blue", "black"), legend.title = "Legenda", xlim = c(-5,5), ylim = c(-3,3))
fviz_pca_ind(pca, axes = c(1,3), col.ind = A$EU, palette = c( "orange", "green3", "red", "blue", "black"), legend.title = "Legenda", xlim = c(-5,5), ylim = c(-3,4.5))
B_s = data.frame(scale(B))
Residui = rep(0,24)
for(i in 1:24){
  B_r = B_s[-i,]
  B_r.pca = princomp(B_r)
  B_p = predict(B_r.pca, newdata = B_s[i,])[1:2]
  Residui[i] = mean((B_p-predict(pca)[i,1:2])^2)
}
sqrt(mean(Residui))
hist(Residui, col = "green3")
round(Residui,2)
B[6,]
B[23,]
summary(B)

