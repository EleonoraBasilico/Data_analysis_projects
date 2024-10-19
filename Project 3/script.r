s = read.csv("tabella.csv")
s
s1 = s$HSN1FNSA[325:696]
t = ts(s1, start = c(1990,1), end = c(2020,12), frequency = 12)
plot(t)

#Analisi preliminare
acf(t,30)
acf(diff(t),30)

#Confrontiamo gli andamenti dei diversi anni per vedere se la stagionalità è stazionaria
m_s = matrix(t,12,31)
par(bg = "black")
ts.plot(scale(m_s, scale = F), col = heat.colors(31))
lines(rowMeans(scale(m_s, scale = F)), lwd = 3, col = "white")
par(bg = "white")

boxplot(t(m_s), pch = "*")

#Decomposizione stazionaria: caso additivo e caso moltiplicativo
t.a = decompose(t)
plot(t.a)
t.m = decompose(t, type = "multiplicative")
plot(t.m)

res.a = as.vector(window(t.a$random, c(1990,7), c(2020,6)))
plot(res.a, pch = 20)
res.m = as.vector(window(t.m$random, c(1990,7), c(2020,6)))
plot(res.m, pch = 20)
res.mlog = log(res.m)
plot(res.mlog, pch = 20)

#Calcoliamo la percentuale di residui in valore assoluto minori di 5 e maggiori di 10 nel caso additivo
max(t.a$seasonal)
min(t.a$seasonal)
j = 0
for(i in 1:length(res.a)){
  if(abs(res.a[i])<5){
    j = j+1
  }
}
perc = j/length(res.a)
perc

j = 0
for(i in 1:length(res.a)){
  if(abs(res.a[i])>10){
    j = j+1
  }
}
perc = j/length(res.a)
perc

#Calcoliamo la percentuale di residui in valore assoluto minori di 0.075 e maggiori di 0.15 nel caso moltiplicativo.
#Per semplificare il conto e per facilitare un confronto con il caso additivo passiamo ai logaritmi.
max(log(t.m$seasonal))
min(log(t.m$seasonal))

j = 0
for(i in 1:length(res.mlog)){
  if(abs(res.mlog[i])<0.075){
    j = j+1
  }
}
perc = j/length(res.mlog)
perc

j = 0
for(i in 1:length(res.mlog)){
  if(abs(res.mlog[i])>0.15){
    j = j+1
  }
}
perc = j/length(res.mlog)
perc

#Confrontiamo più in dettaglio i residui dei modelli additivo e moltiplicativo.
acf(res.a, 30)
acf(res.mlog,30)
sd(acf(res.a, plot = F)$acf)
sd(acf(res.m, plot = F)$acf)
var(res.a)/var(window(t, c(1990,7), c(2020,6)))
var(res.mlog)/var(window(log(t), c(1990,7), c(2020,6)))

layout(t(1:2))
hist(res.a, 20 , freq = F)
lines(density(res.a), col="blue")
lines(sort(res.a),dnorm(sort(res.a), mean(res.a), sd(res.a)), col="red")
hist(res.mlog, 20 , freq=F)
lines(density(res.mlog), col="blue")
lines(sort(res.mlog),dnorm(sort(res.mlog), mean(res.mlog), sd(res.mlog)), col="red")
qqnorm(res.a)
qqline(res.a)
qqnorm(res.mlog)
qqline(res.mlog)
layout(1)
shapiro.test(res.a)
shapiro.test(res.mlog)

#Decomposizione non stazionaria: caso additivo e caso moltiplicativo
plot(stl(t, 7))
res.stl = stl(t,7)$time.series[,3]
acf(res.stl,30)
sd(acf(res.stl, plot = F)$acf)

plot(stl(t, 15))
res.stl = stl(t, 15)$time.series[,3]
acf(res.stl,30)
sd(acf(res.stl, plot = F)$acf)

plot(stl(log(t),7))
res.stl2 = stl(log(t),7)$time.series[,3]
acf(res.stl2,30)
sd(acf(res.stl2, plot = F)$acf)

plot(stl(log(t),15))
res.stl2 = stl(log(t),15)$time.series[,3]
acf(res.stl2,30)
sd(acf(res.stl2, plot = F)$acf)

#Metodo di Holt-Winters
t.hw = HoltWinters(t, seasonal = "multiplicative")
plot(t.hw)
t.hw$alpha
t.hw$beta
t.hw$gamma

#Pur cambiando intercetta e coefficiente angolare iniziali vediamo che il grafico rimane pressochè invariato
x = 1:20
coefficients(lm(t[1:20]~x))
plot(HoltWinters(t, seasonal = "multiplicative", l.start = 48.18, b.start = -0.36))

#Vogliamo scegliere dei valori per i parametri alpha, beta e gamma tali per cui si abbiano meno 
#struttura nei residui e maggior capacità di predizione. A tal fine costruiamo la matrice A, le
#cui colonne indicano rispettivamente i valori di alpha, beta, gamma, deviazione standard della
#acf dei residui del metodo HW corrispondente a quei parametri, capacità di predizione del modello,
#percentuale di varianza non spiegata e p-value del test di Shapiro.
#Abbiamo considerato per ogni parametro 5 diversi valori che non si discostano troppo dal valore
#ottimale.
A = matrix(0,125,7)
i = 1
train = window(t, end = c(2019,12))
test = window(t, 2020)
for( a in c(45,50,55,60,65)){
  for(b in c(0.45,0.50,0.55,0.60,0.65)){
    for(c in c(40,45,50,55,60)){
      A[i,1] = a/100
      A[i,2] = b/1000
      A[i,3] = c/100
      residui = resid(HoltWinters(t, seasonal = "multiplicative", alpha = a/100, beta = b/1000, gamma = c/100))
      A[i,4] = sd(acf(residui, plot = F)$acf)
      T.hw = HoltWinters(train, seasonal = "multiplicative", alpha = a/100, beta = b/1000, gamma = c/100)
      T.hw.p = predict(T.hw,12)
      A[i,5] = sqrt(mean((T.hw.p-test)^2))
      A[i,6] = var(residui)/var(window(t,1991))
      A[i,7] = shapiro.test(residui)$p.value
      i = i+1
    }
  }
}
A

#Riportiamo i grafici relativi ai valori dei parametri alpha, beta e gamma in corrispondenza dei quali
#si hanno meno struttura nei residui e maggior capacità di predizione. Abbiamo preso in considerazione 
#questi due parametri tralasciando la percentuale di varianza non spiegata e il p-value del test di Shapiro
#in quanto le differenze tra i valori sono meno significative (nel primo caso parliamo della terza cifra decimale).
plot(HoltWinters(t, seasonal = "multiplicative", alpha = 0.5, beta = 0.00065, gamma = 0.55))
plot(HoltWinters(t, seasonal = "multiplicative", alpha = 0.45, beta = 0.00065, gamma = 0.45))

#Per scegliere la terna di parametri migliore tra le due appena usate ne studiamo in dettaglio
#residui e capacità di predizione.
hw1 = HoltWinters(t, seasonal = "multiplicative", alpha = 0.5, beta = 0.00065, gamma = 0.55)
hw2 = HoltWinters(t, seasonal = "multiplicative", alpha = 0.45, beta = 0.00065, gamma = 0.45)

hw1.r = resid(hw1)
hw2.r = resid(hw2)
var(hw1.r)/var(window(t,1991))
var(hw2.r)/var(window(t,1991))
layout(t(1:2))
plot(hw1.r, type = "p", pch = 20)
plot(hw2.r, type = "p", pch = 20)
plot(as.numeric(hw1$fitted[,1]), as.numeric(hw1.r), type = "p", pch = 20)
plot(as.numeric(hw2$fitted[,1]), as.numeric(hw2.r), type = "p", pch = 20)
acf(hw1.r,30)
acf(hw2.r,30)
sd(acf(hw1.r, plot = F)$acf)
sd(acf(hw2.r, plot = F)$acf)
hist(hw1.r, 20, freq = F)
lines(density(hw1.r), col = "blue")
lines(sort(hw1.r), dnorm(sort(hw1.r), mean(hw1.r), sd(hw1.r)), col = "red")
hist(hw2.r, 20, freq = F)
lines(density(hw2.r), col = "blue")
lines(sort(hw2.r), dnorm(sort(hw2.r), mean(hw2.r), sd(hw2.r)), col = "red")
qqnorm(hw1.r, pch = 20)
qqline(hw1.r)
qqnorm(hw2.r, pch = 20)
qqline(hw2.r)
shapiro.test(hw1.r)
shapiro.test(hw2.r)
layout(1)

train = window(t, end = c(2019,12))
test = window(t, 2020)
hw1.2 = HoltWinters(train, seasonal = "multiplicative", alpha = 0.5, beta = 0.00065, gamma = 0.55)
hw2.2 = HoltWinters(train, seasonal = "multiplicative", alpha = 0.45, beta = 0.00065, gamma = 0.45)
hw1.p = predict(hw1.2,12)
hw2.p = predict(hw2.2,12)
sqrt(mean((hw1.p-test)^2))
sqrt(mean((hw2.p-test)^2))
ts.plot(test, hw1.p, col = c("black", "red"))
ts.plot(test, hw2.p, col = c("black", "red"))

l = length(t)
res.hw1 = rep(0,24)
res.hw2 = rep(0,24)
predizione.hw1 = rep(0, 24)
predizione.hw2 = rep(0, 24)
j=1
for(i in (l-24):(l-1)){
  t_cv = ts(t[1:i], frequency = 12, start = c(1990,1))
  hw1.3 = HoltWinters(t_cv, seasonal="multiplicative", alpha = 0.5, beta = 0.00065, gamma = 0.55)
  hw2.3 = HoltWinters(t_cv, seasonal = "multiplicative", alpha = 0.45, beta = 0.00065, gamma = 0.45)
  hw1.3.p = predict(hw1.3,1)
  hw2.3.p = predict(hw2.3,1)
  res.hw1[j] = hw1.3.p - t[i+1]
  res.hw2[j] = hw2.3.p - t[i+1]
  predizione.hw1[j] = hw1.3.p
  predizione.hw2[j] = hw2.3.p
  j=j+1
}
sqrt(mean(res.hw1^2))
sqrt(mean(res.hw2^2))
plot(res.hw1, type = "b", pch = 20, col = "blue")
lines(res.hw2, type = "b", pch = 20, col = "green3")
ts.plot(window(t, c(2019,1), c(2020,12)), predizione.hw1, col = c("black", "red"))
ts.plot(window(t, c(2019,1), c(2020,12)), predizione.hw2, col = c("black", "red"))

#Scegliamo la prima terna di valori e facciamo una previsione per l'anno 2021.
plot(hw1, predict(hw1, 12))
lines(predict(hw1,12)+quantile(hw1.r, 0.05), col = "blue")
lines(predict(hw1,12)+quantile(hw1.r, 0.95), col = "blue")

#La scarsa capacità di predizione del modello (predizione a 12 mesi) è dovuta al fatto che l'andamento delle vendite nel 2020
#è completamente diverso rispetto all'andamento degli anni precedenti. Se infatti applichiamo il 
#metodo di HW sulla serie troncata al 2018 e prevediamo l'andamento delle vendite del 2019 otteniamo
# un risultato ben migliore in termini di autovalidazione.
train = window(t, end = c(2018,12))
test = window(t, c(2019,1), c(2019,12))
hw1.3 = HoltWinters(train, seasonal = "multiplicative", alpha = 0.5, beta = 0.00065, gamma = 0.55)
hw3.p = predict(hw1.3,12)
sqrt(mean((hw3.p-test)^2))
ts.plot(test, hw3.p, col = c("black", "red"))

#Metodi autoregressivi: modello diretto
pacf(t)

l = length(t)
L = 17  # numero di lag in ingresso
mt = matrix(nrow = l - L, ncol = L + 1)
for (i in 1:(L + 1)) {
  mt[, i] = t[i:(l - L - 1 + i)]
}
mt = data.frame(mt)
t.lm = lm(X18 ~ ., data = mt)
summary(t.lm)

t.lm <- lm(X18 ~ . -X15-X9-X8-X4-X2-X16-X3-X10-X11-X12-X14-X7-X13-X1, data = mt)
summary(t.lm)

#Predizione e analisi del modello diretto ridotto
anni = 1
pt = rep(0, l + 12 * anni)
pt[1:l] = t
for (i in 1:(12 * anni)) {
  pt[l + i] = coef(t.lm) %*% c(1, pt[l + i - 13], pt[l+i-12], pt[l + i - 1])
}
t.lm.pt = ts(pt, frequency = 12, start = c(1990, 1))
t.lm.a = window(t, c(1991, 6)) - resid(t.lm)
ts.plot(t, t.lm.a, window(t.lm.pt, c(2020, 12)), col = c("black","blue","red"))

#Analisi dei residui
t.lm.r = resid(t.lm)
var(t.lm.r)/var(window(t, c(1991,6)))
plot(t.lm.r, type = "p", pch = 20)
plot(t.lm.a, t.lm.r, type = "p", pch = 20)
acf(t.lm.r,30)
sd(acf(t.lm.r, plot = F)$acf)
pacf(t.lm.r)
hist(t.lm.r, 20, freq = F)
lines(density(t.lm.r), col = "blue")
lines(sort(t.lm.r), dnorm(sort(t.lm.r), mean(t.lm.r), sd(t.lm.r)), col = "red")
qqnorm(t.lm.r, pch = 20)
qqline(t.lm.r)
shapiro.test(t.lm.r)

#Autovalidazione
train = window(t, end = c(2019,12))
test = window(t, 2020)

l1 = length(train)
L = 17
mt2 = matrix(nrow = l1 - L, ncol = L + 1)
for (i in 1:(L + 1)) {
  mt2[, i] = train[i:(l1 - L - 1 + i)]
}
mt2 = data.frame(mt2)
t.lm2 = lm(X18 ~ . -X15-X9-X8-X4-X2-X16-X3-X10-X11-X12-X14-X7-X13-X1, data = mt2)
summary(t.lm2)

pt2 = rep(0, l1 + 12)
pt2[1:l1] = train
for (i in 1:12) {
  pt2[l1 + i] = coef(t.lm2) %*% c(1, pt2[l1 + i - 13], pt2[l1 + i - 12], pt2[l1 + i - 1])
}
t.lm.pt2 = ts(pt2, frequency = 12, start = c(1990, 1))
sqrt(mean((window(t.lm.pt2, 2020) - test)^2))
ts.plot(test, t.lm.pt2, col = c("black", "red"))

#Analisi e previsione del modello diretto non ridotto
t.lm1 = lm(X18 ~ ., data = mt)

anni = 1
pt = rep(0, l + 12 * anni)
pt[1:l] = t
for (i in 1:(12 * anni)) {
  pt[l + i] = coef(t.lm1) %*% c(1, rev(pt[l+i-1:L]))
}
t.lm1.pt = ts(pt, frequency = 12, start = c(1990, 1))
t.lm1.a = window(t, c(1991, 6)) - resid(t.lm1)
ts.plot(t, t.lm1.a, window(t.lm1.pt, c(2020, 12)), col = c("black","blue","red"))

#Analisi dei residui del modello diretto non ridotto
t.lm1.r = resid(t.lm1)
var(t.lm1.r)/var(window(t, c(1991,6)))
plot(t.lm1.r, type = "p", pch = 20)
plot(t.lm1.a, t.lm1.r, type = "p", pch = 20)
acf(t.lm1.r,30)
sd(acf(t.lm1.r, plot = F)$acf)
pacf(t.lm1.r)
hist(t.lm1.r, 20, freq = F)
lines(density(t.lm1.r), col = "blue")
lines(sort(t.lm1.r), dnorm(sort(t.lm1.r), mean(t.lm1.r), sd(t.lm1.r)), col = "red")
qqnorm(t.lm1.r, pch = 20)
qqline(t.lm1.r)
shapiro.test(t.lm1.r)

#Autovalidazione
train = window(t, end = c(2019,12))
test = window(t, 2020)

l1 = length(train)
L = 17
mt2 = matrix(nrow = l1 - L, ncol = L + 1)
for (i in 1:(L + 1)) {
  mt2[, i] = train[i:(l1 - L - 1 + i)]
}
mt2 = data.frame(mt2)
t.lm2 = lm(X18 ~ ., data = mt2)
summary(t.lm2)

pt2 = rep(0, l1 + 12)
pt2[1:l1] = train
for (i in 1:12) {
  pt2[l1 + i] = coef(t.lm2) %*% c(1, rev(pt2[l1+i-1:L]))
}
t.lm1.pt2 = ts(pt2, frequency = 12, start = c(1990, 1))
sqrt(mean((window(t.lm1.pt2, 2020) - test)^2))
ts.plot(test, t.lm1.pt2, col = c("black", "red"))

#Metodo Yule-Walker
t.ar = ar(t)
t.ar
ts.plot(t, t - t.ar$resid, col = c("black", "red"))
t.ar.pt = predict(t.ar, n.ahead = 24, se.fit = FALSE)
plot(t.ar.pt)

#Analisi dei residui
Y.W.res = na.omit(t.ar$resid)
var(Y.W.res)/var(window(t, c(1991,7)))
plot(Y.W.res, type = "p", pch = 20)
acf(Y.W.res)
sd(acf(Y.W.res, plot = F)$acf)
pacf(Y.W.res)
hist(Y.W.res, 20, freq = F)
lines(density(Y.W.res), col = "blue")
lines(sort(Y.W.res), dnorm(sort(Y.W.res), mean(Y.W.res), sd(Y.W.res)), col = "red")
qqnorm(Y.W.res, pch = 20)
qqline(Y.W.res)
shapiro.test(Y.W.res)

#Autovalidazione
train = window(t, end = c(2019, 12))
test = window(t, 2020)
YW.p = predict(ar(train), n.ahead = 12, se.fit = FALSE)
ts.plot(test, YW.p, col = c("black", "red"))
sqrt(mean((test - YW.p)^2))

res.ar = rep(0, 11)
for (i in 1:11) {
  train = window(t, end = c(2020, i))
  test = window(t, start = c(2020, i + 1))
  res.ar[i] = predict(ar(train), n.ahead = 1, se.fit = F)
}
test = window(t, start = c(2020, 2))
sqrt(mean((test - res.ar)^2))
ts.plot(test, res.ar, col = c("black", "blue"))

#Metodo dei minimi quadrati
t.ls = ar(t, method = "ols")
t.ls$order
ts.plot(t, t - t.ls$resid, col = c("black", "blue"))

#Analisi dei residui
t.ls.r = as.double(na.omit(t.ls$resid))
t.ls.fitted = as.double(na.omit(t - t.ls$resid))
plot(t.ls.r, pch = 20)
plot(t.ls.fitted, t.ls.r, pch = 20)
var(t.ls.r)/var(t.ls.r + t.ls.fitted)
acf(t.ls.r)
sd(acf(t.ls.r, plot = F)$acf)
pacf(t.ls.r)
hist(t.ls.r, 20, freq = F)
lines(density(t.ls.r), col = "red")
lines(sort(t.ls.r), dnorm(sort(t.ls.r), mean(t.ls.r), sd(t.ls.r)), col = "blue")
qqnorm(t.ls.r, pch = 20)
qqline(t.ls.r)
shapiro.test(t.ls.r)

#Previsione
t.ls.pt = predict(t.ls, n.ahead = 24, se.fit = FALSE)
plot(t.ls.pt)
ts.plot(window(t, 1990), window(t - t.ls$resid, 1990), t.ls.pt, col = c("black", "blue", "red"), lwd = c(1,1,1))
# stima empirica dell'incertezza
lines(t.ls.pt + quantile(t.ls.r, 0.975), col = "green3")
lines(t.ls.pt + quantile(t.ls.r, 0.025), col = "green3")
segments(2021,10,2021,150,lty=3,col="black")

#Autovalidazione
train = window(t, end = c(2019, 12))
test = window(t, 2020)
tcv.ls.p = predict(ar(train, method = "ols"), n.ahead = 12, se.fit = FALSE)
ts.plot(test, tcv.ls.p, col = c("black", "red"))
sqrt(mean((tcv.ls.p - test)^2))

res.ar2 = rep(0, 11)
for (i in 1:11) {
  train = window(t, end = c(2020, i))
  test = window(t, start = c(2020, i + 1))
  res.ar2[i] = predict(ar(train, method = "ols"), n.ahead = 1, se.fit = F)
}
test = window(t, start = c(2020, 2))
sqrt(mean((test - res.ar2)^2))
ts.plot(test, res.ar2, col = c("black", "blue"))

#Il metodo dei minimi quadrati sembra essere migliore degli altri, seppur di poco.
#Confrontiamolo con il metodo HW in termini di analisi e previsione.
ts.plot(t, t - t.ls$resid, t.hw$fitted[, 1], col = c("black", "red", "green3"))
ts.plot(t, t.ls.pt, predict(t.hw, 12), col = c("black", "red", "green3"))

#Confrontiamo la capacità di predizione dei due metodi sugli undici mesi
res.ar2 = rep(0, 11)
res.hw = rep(0, 11)
for (i in 1:11) {
  train = window(t, end = c(2020, i))
  test = window(t, start = c(2020, i + 1))
  res.ar2[i] = predict(ar(train, method = "ols"), n.ahead = 1, se.fit = F)
  res.hw[i] = predict(HoltWinters(train, seasonal = "multiplicative"),1)
}
test = window(t, start = c(2020, 2))
sqrt(mean((test - res.ar2)^2))
sqrt(mean((test - res.hw)^2))
ts.plot(test, res.ar2, res.hw, col = c("black", "blue", "red"))

#Confrontiamo le previsioni dei 5 medoti per il 2021 con i valori noti (fino a novembre).
plot(s$HSN1FNSA[697:707], col = "green3", type = "l", ylim = c(10,100), xlab = "Mesi", ylab = "Vendite (unità di migliaia)")
lines(as.vector(predict(hw1, 11)), col = "red", type = "l")
lines(as.vector(window(t.lm.pt, 2021))[1:11], col = "blue", type = "l")
lines(as.vector(window(t.lm1.pt, 2021))[1:11], col = "black", type = "l")
lines(as.vector(window(t.ar.pt, c(2021,1), c(2021,12)))[1:11], col = "purple", type = "l")
lines(as.vector(window(t.ls.pt, c(2021,1), c(2021,12)))[1:11], col = "orange", type = "l")

legend("bottomright",legend = c("2021", "H-W", "Mod Dir r.", "Mod Dir", "Y-W", "Min Quad"), col = c("green3", "red", "blue", "black", "purple", "orange"), lty = c(1,1,1,1,1,1))