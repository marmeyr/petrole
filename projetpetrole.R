############################################################################### 
### Projet Petrole
############################################################################### 
# Auteur : Alexandra MILLOT
# UE : Econométrie des séries temporelles

############################################################################### 
### Importation des données
############################################################################### 
rm(list=ls())
tab<-read.csv("petrole.csv",header=TRUE,dec=".",sep=",")
str(tab)
tab2=subset(tab,select=c(LOCATION,TIME,Value))
canada=tab2[tab2$LOCATION=="CAN",]
canada=na.omit(canada)
canada
Petrole<-ts(canada$Value,start=1971,freq=1)

# Créer un dataframe pour ggplot
df <- data.frame(Year = time(Petrole), Value = as.numeric(Petrole))

############################################################################### 
### Importation des libraries
############################################################################### 
library(forecast);library(caschrono)
library(lmtest);library(urca)
library(CADFtest);library(ggplot2)
library(doSNOW);library(parallel)
library(TSA);library(FinTS)

############################################################################### 
### Chronogramme de la série d'origine
############################################################################### 
# Tracer le graphique avec ggplot2
ggplot(df, aes(x = Year, y = Value)) +
  geom_line() +
  labs(x = 'Année', y = 'Production de pétrole non raffiné', title = 'Production de pétrole non raffiné de 1970 à 2021 au Canada') +
  # Ajouter la droite de régression
  geom_smooth(method = 'lm', formula = y ~ x, color = 'red', se = FALSE) 

df <- data.frame(Year = time(Petrole), Value = as.numeric(Petrole))

############################################################################### 
### ACF et PACF
############################################################################### 
opp <- par(mfrow = c(1, 2))
Acf(Petrole, lag = 30)
Pacf (Petrole, lag = 30)
par(opp)

######################################################################## 
# Test de racine unitaire : Dickey-Fuller
######################################################################## 
summary(ur.df(Petrole,type="trend",lag=0))
summary(ur.df(Petrole,type="drift",lag=0))
summary(ur.df(Petrole,type="none",lag=0)) 
# DS d'après Dickey-Fuller : "none"

plot(ur.df(Petrole, type="none",lag=0))
# aléa autocorrélé ? Non, pas besoin de faire le Dickey-Fuller Augmenté

######################################################################## 
# Test de racine unitaire : Zivot-Andrews
######################################################################## 
# Formule de Schwert pour trouver le pmax
TPetrole <- length(Petrole)
pmax<-as.integer(12*(TPetrole/100)^(0.25))
pmax # 10

summary(ur.za(Petrole, model="both",lag=pmax))
summary(ur.za(Petrole, model="both",lag=9))
summary(ur.za(Petrole, model="both",lag=8))
# Potential break point at position: 41
# -4.11 > la valeur critique à 5% = -5.08 donc on accepte H0
# On en conclue DS sans changement structurel

######################################################################## 
# Test de racine unitaire : Lee-Strazicich
######################################################################## 
y <- Petrole
source("LeeStrazicichUnitRoot-master/LeeStrazicichUnitRootTest.R")
myBreaks <- 1
myModel <- "break"
myLags <- 4 # car 5 renvoie une erreur

myLS_test <-ur.ls(y=y , model = myModel, breaks = myBreaks, lags = myLags, method = "GTOS",pn = 0.1, print.results = "print" )
# First possible structural break at position: 32 (2002)
# -8.917015 < -4.51 donc vous rejetez H0 donc le PGD qui a généré le Petrole est 
# TS avec un changement structurel dans la constante et la pente de la partie 
# déterminsite de la série
# on prend la valeur de 5% pour le lambda 0.5 car je n'ai pas accès au lambda 0.6 

myBreaks <- 2
myLS_test <-ur.ls(y=y , model = myModel, breaks = myBreaks, lags = myLags, method = "GTOS",pn = 0.1, print.results = "print" )
# First possible structural break at position: 20 (TB1 = 20) (1990)
# Second possible structural break at position: 30 (TB2 = 30) (2000)
# -11.57808 < -5.67 (La valeur à 5% de Break 1 - 0.4 et Break 2 - 0.6)
# Donc vous rejetez H0 donc le PGD qui a généré le Pétrole est TS avec changements structurels.

# Comme le test avec 2 dates est moins puissant que celui avec 1 seule date de 
# rupture vous garderez la conclusion de LS avec 1 seule date

######################################################################## 
# Test de racine unitaire : Lee-Strazicich avec Bootstrap (tirage avec remise)
########################################################################
source("LeeStrazicichUnitRoot-master/LeeStrazicichUnitRootTestParallelization.R")
#Define number of cores to use. By default the maximum available number minus one core is used
cl <- makeCluster(max(1, detectCores() - 1))
registerDoSNOW(cl)
myBreaks <- 1
# Lancer le code du fichier LeeStrazicichUnitRootTestParallelization.R
myParallel_LS <- ur.ls.bootstrap(y=Petrole , model = myModel, breaks = myBreaks, lags = myLags, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "print")
# First possible structural break at position: 35
# -3.288709 > -4.51 donc on accepte H0 et le Petrole qui a généré notre série est DS

# Comme on travaille sur des échantillons assez petits nous allons privilégier les 
# résultats du test LS avec bootstrap donc c’est cette conclusion que je garde.

myBreaks <- 2
myParallel_LS <- ur.ls.bootstrap(y=Petrole , model = myModel, breaks = myBreaks, lags = myLags, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "print")
# First possible structural break at position: 35
# -4.512272 > -5.72

######################################################################## 
# Différenciation à l'ordre 1
########################################################################
dPetrole=diff(Petrole)
eacf(dPetrole)

ddf <- data.frame(Year = time(dPetrole), Value = as.numeric(dPetrole))

############################################################################### 
### Chronogramme de la série d'origine
############################################################################### 
# Tracer le graphique avec ggplot2
ggplot(ddf, aes(x = Year, y = Value)) +
  geom_line() +
  labs(x = 'Année', y = 'Production de pétrole non raffiné', title = 'Production de pétrole non raffiné de 1970 à 2021 au Canada') +
  # Ajouter la droite de régression
  geom_smooth(method = 'lm', formula = y ~ x, color = 'red', se = FALSE) 

######################################################################## 
# ACF et PACF
######################################################################## 
opp <- par(mfrow = c(1, 2))
Acf(dPetrole, lag = 80)
Pacf (dPetrole)
par(opp)

######################################################################## 
# Test de racine unitaire : Dickey-Fuller
######################################################################## 
summary(ur.df(dPetrole,type="trend",lag=0))
# -6.3535 < -3.50 
# Avec Dickey-Fuller le PGD est TS

plot(ur.df(dPetrole, type="trend",lag=0))
# aléa autocorrélé ? Non, pas besoin de faire le Dickey-Fuller Augmenté

######################################################################## 
# Test de racine unitaire : Zivot-Andrews
######################################################################## 
# Formule de Schwert pour trouver le pmax
TdPetrole <- length(dPetrole)
pmax<-as.integer(12*(TdPetrole/100)^(0.25))
pmax # 10

summary(ur.za(dPetrole, model="both",lag=pmax))
summary(ur.za(dPetrole, model="both",lag=9))
summary(ur.za(dPetrole, model="both",lag=8))
summary(ur.za(dPetrole, model="both",lag=7))
# Potential break point at position: 35 (2005)
# -5.5161 < la valeur critique à 5% = -5.08 donc on rejette H0
# On en conclue TS avec un changement structurel

plot(ur.za(dPetrole, model="both",lag=7))

######################################################################## 
# Test de racine unitaire : Lee-Strazicich
######################################################################## 
y <- dPetrole
myBreaks <- 1
myModel <- "break"
myLags <- 4 # car 5 renvoie une erreur

myLS_test <-ur.ls(y=y , model = myModel, breaks = myBreaks, lags = myLags, method = "GTOS",pn = 0.1, print.results = "print" )
# First possible structural break at position: 12 (1982)
# -6.521028 < -4.47 donc vous rejetez H0 donc le PGD qui a généré le dPetrole est 
# TS avec un changement structurel dans la constante et la pente de la partie 
# déterminsite de la série

myBreaks <- 2
myLS_test <-ur.ls(y=y , model = myModel, breaks = myBreaks, lags = myLags, method = "GTOS",pn = 0.1, print.results = "print" )
# First possible structural break at position: 12 (TB1 = 20) (1982)
# Second possible structural break at position: 40 (TB2 = 30) (2010)
# -6.778718 < -5.71 (La valeur à 5% de Break 1 - 0.4 et Break 2 - 0.6)
# Donc vous rejetez H0 donc le PGD qui a généré le Pétrole est TS avec changements structurels.

# Comme le test avec 2 dates est moins puissant que celui avec 1 seule date de 
# rupture vous garderez la conclusion de LS avec 1 seule date

######################################################################## 
# Test de racine unitaire : Lee-Strazicich avec Bootstrap (tirage avec remise)
########################################################################
#Define number of cores to use. By default the maximum available number minus one core is used
cl <- makeCluster(max(1, detectCores() - 1))
registerDoSNOW(cl)
myBreaks <- 1
# Lancer le code du fichier LeeStrazicichUnitRootTestParallelization.R
myParallel_LS <- ur.ls.bootstrap(y=dPetrole , model = myModel, breaks = myBreaks, lags = myLags, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "print")
# First possible structural break at position: 7 (1977)
# -3.991365 > -4.50 donc on accepte H0 et le Petrole qui a généré notre série est DS

# Comme on travaille sur des échantillons assez petits nous allons privilégier les 
# résultats du test LS avec bootstrap donc c’est cette conclusion que je garde.

myBreaks <- 2
myParallel_LS <- ur.ls.bootstrap(y=dPetrole , model = myModel, breaks = myBreaks, lags = myLags, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "print")
# First possible structural break at position: 7 (1977)
# Second possible structural break at position: 41 (2011)
# -5.491944 < -5.71

############################################################################### 
### Différenciation à l'ordre 2
###############################################################################
ddPetrole=diff(Petrole, 2)
ddf <- data.frame(Year = time(ddPetrole), Value = as.numeric(ddPetrole))

############################################################################### 
### Chronogramme de la série différenciée
############################################################################### 
# Tracer le graphique avec ggplot2
ggplot(ddf, aes(x = Year, y = Value)) +
  geom_line() +
  labs(x = 'Année', y = 'Production de pétrole non raffiné', title = ("Production de pétrole non raffiné de 1970 à 2021 au Canada - série différenciée à l'ordre 2")) +
  # Ajouter la droite de régression
  geom_smooth(method = 'lm', formula = y ~ x, color = 'red', se = FALSE) 

############################################################################### 
### ACF et PACF
############################################################################### 
opp <- par(mfrow = c(1, 2))
Acf(ddPetrole, lag = 30)
Pacf (ddPetrole, lag = 30)
par(opp)

############################################################################### 
### Test de racine unitaire : Dickey-Fuller
############################################################################### 
summary(ur.df(ddPetrole,type="trend",lag=0)) 
summary(ur.df(ddPetrole,type="drift",lag=0))
summary(ur.df(ddPetrole,type="none",lag=0))
# -9.569 < -1.95 
# Avec Dickey-Fuller le PGD est AR(1)

plot(ur.df(ddPetrole, type="none",lag=0))
# aléa autocorrélé ? Oui, besoin de faire le Dickey-Fuller Augmenté

############################################################################### 
### Test de racine unitaire : Dickey-Fuller Augmenté
############################################################################### 
# Formule de Schwert pour trouver le pmax
TddPetrole <- length(ddPetrole)
pmax<-as.integer(12*(TddPetrole/100)^(0.25))
pmax # 10

summary(CADFtest(ddPetrole,criterion="MAIC",type="none",max.lag.y=pmax))
summary(ur.df(ddPetrole,type="none",lag=pmax,selectlag="BIC"))

# On rejet H0,

# Ce processus est stationnaire d'après Dickey-Fuller Augmenté. Nous allons 
# maintenant procéder à la recherche du meilleur modèle afin d'effectuer les 
# prédiction.

############################################################################### 
### L'EACF
############################################################################### 
eacf(ddPetrole) # p = 0 et q = 1 sachant que d = 2

############################################################################### 
### Modèle testé 1
############################################################################### 
reg1 = Arima(Petrole, order = c(0, 2, 1), include.constant = T)
coeftest(reg1) # Le coefficient ma1 est significatif.
residu<-reg1$res
jarque.bera.test(residu) # p-value = 0.05943
# p-value > 0.05 : On accepte H0, les aléas sont gaussiens
t.test(residu) # p-value = 0.7686
# p-value > 0.05 : On accepte H0, les aléas du ARIMA(0,2,1) ont une espérance nulle.
residu=(residu-mean(residu))/sd(residu)
K<-12
tmp<-rep(0,K)
temp<-rep(0,K)
for(i in 1:K){
  tmp[i]<-Box.test(residu,lag=i,type="Ljung-Box")$p.value
  temp[i]<-ArchTest(residu, lag=i)$p.value
}
tmp
temp
# Toutes les p-values sont > à 0.05 donc on ne rejette aucune des deux hypothèses 
# nulles pour tout ordre de retard inférieur ou égal à 12. 
# Donc les aléas sont des BB.

reg1$bic # 1011.889

prev<-forecast(reg1,h=4,level=0.95)
prev
plot(prev)
lines(Petrole,col=2)

############################################################################### 
### Modèle testé 2
############################################################################### 
# modèle de la série Petrole
reg2 = Arima(Petrole, order = c(1, 2, 1), include.constant = T)
coeftest(reg2) 
# Le coefficient ar1 n'est pas significatif. Si on estime pas ce coefficient,
# cela revient au modèle 1.

# Nous retrouvons la même conclusion pour les ARIMA(2, 2, 1), ARIMA(3, 2, 1),
# ARIMA(4, 2, 1), ARIMA(5, 2, 1) et ARIMA(6, 2, 1)

############################################################################### 
### Modèle testé 3
############################################################################### 
# modèle de la série Petrole
reg3 = Arima(Petrole, order = c(1, 2, 2), include.constant = T)
coeftest(reg3) # Le coefficient ma1 n'est pas significatif.
reg3 = Arima(Petrole, order = c(1, 2, 2), include.constant = T,
             fixed = c(NA, 0, NA))
coeftest(reg3) # Tous les coefficients sont significatifs
residu<-reg3$res
jarque.bera.test(residu) # p-value = 0.07229
# p-value > 0.05 : On accepte H0, les aléas sont gaussiens
t.test(residu) # p-value = 0.7454
# p-value > 0.05 : On accepte H0, les aléas du ARIMA(1,2,2) ont une espérance nulle.
residu=(residu-mean(residu))/sd(residu)
K<-12
tmp<-rep(0,K)
temp<-rep(0,K)
for(i in 1:K){
  tmp[i]<-Box.test(residu,lag=i,type="Ljung-Box")$p.value
  temp[i]<-ArchTest(residu, lag=i)$p.value
}
tmp
temp
# Toutes les p-values sont > à 0.05 donc on ne rejette aucune des deux hypothèses 
# nulles pour tout ordre de retard inférieur ou égal à 12. 
# Donc les aléas sont des BB.

reg3$bic # 1015.639
# Le bic est supérieur à celui du premier modèle donc nous conservons le modèle 1

############################################################################### 
### Modèle testé 4
############################################################################### 
# modèle de la série Petrole
reg4 = Arima(Petrole, order = c(1, 2, 3), include.constant = T)
coeftest(reg4) # Les coefficients ma1 et ma3 ne sont pas significatifs.
# Si on estime pas ces coefficients, cela revient au modèle 3.

# Nous retrouvons la même conclusion pour les ARIMA(1, 2, 4), ARIMA(1, 2, 5) et
# ARIMA(1,2,6)

############################################################################### 
### Modèle testé 5
############################################################################### 
# modèle de la série Petrole
reg5 = Arima(Petrole, order = c(2, 2, 3), include.constant = T)
coeftest(reg5) # Tous les coefficients sont significatifs.
residu<-reg5$res
jarque.bera.test(residu) # p-value = 0.2964
# p-value > 0.05 : On accepte H0, les aléas sont gaussiens
t.test(residu) # p-value = 0.6282
# p-value > 0.05 : On accepte H0, les aléas du ARIMA(2,2,3) ont une espérance nulle.
residu=(residu-mean(residu))/sd(residu)
K<-12
tmp<-rep(0,K)
temp<-rep(0,K)
for(i in 1:K){
  tmp[i]<-Box.test(residu,lag=i,type="Ljung-Box")$p.value
  temp[i]<-ArchTest(residu, lag=i)$p.value
}
tmp
temp
# Toutes les p-values sont > à 0.05 donc on ne rejette aucune des deux hypothèses 
# nulles pour tout ordre de retard inférieur ou égal à 12. 
# Donc les aléas sont des BB.

reg5$bic # 1024.17
# Le bic est supérieur à celui du premier modèle donc nous conservons le modèle 1


reg5$bic

############################################################################### 
### Modèle testé 6
############################################################################### 
# modèle de la série Petrole
reg6 = Arima(Petrole, order = c(3, 2, 6), include.constant = T)
coeftest(reg6) # Seulement les coefficents associés à ar3 et ma1 sont
# significatifs
reg6 = Arima(Petrole1, order = c(3, 2, 6), include.constant = T,
             fixed = c(0, 0, NA, NA, 0, 0, 0, 0, 0))
coeftest(reg6) # Seulement ma1 est significatif on retombe sur le modèle 1

############################################################################### 
### Modèle testé 7
############################################################################### 
# modèle de la série Petrole
reg7 = Arima(Petrole, order = c(3, 2, 2), include.constant = T)
coeftest(reg7) # Aucun coefficient n'est significatif, on abandonne ce modèle.

############################################################################### 
### Modèle testé 8
############################################################################### 
# modèle de la série Petrole
reg8 = Arima(Petrole, order = c(1, 2, 0), include.constant = T)
coeftest(reg8) # Le coefficient ar1 est significatif.
residu<-reg8$res
jarque.bera.test(residu) # p-value = 0.1235
# p-value > 0.05 : On accepte H0, les aléas sont gaussiens
t.test(residu) # p-value = 0.9093
# p-value > 0.05 : On accepte H0, les aléas du ARIMA(1,2,0) ont une espérance nulle.
residu=(residu-mean(residu))/sd(residu)
K<-12
tmp<-rep(0,K)
temp<-rep(0,K)
for(i in 1:K){
  tmp[i]<-Box.test(residu,lag=i,type="Ljung-Box")$p.value
  temp[i]<-ArchTest(residu, lag=i)$p.value
}
tmp
temp
# Pas toutes les p-values sont > à 0.05 donc il y a des effets ARCH, ce modèle
# n'est donc pas adéquat.

reg8$bic # 1024.83
# Le bic est supérieur à celui du premier modèle donc nous conservons le modèle 1

###
#### CONCLUSION : le meilleur modèle est le ARIMA(0,2,1)
###

############################################################################### 
## Comparer les valeurs réelles au valeurs prédites pur les années 2022 et 2023
###############################################################################
# Source : https://www.statcan.gc.ca/o1/fr/plus/5781-la-production-de-petrole-brut-atteint-des-sommets-surtout-en-raison-des-sables-bitumineux
# La variation d’une année à l’autre est de +1,9 % pour 2022 et de +2,1 % pour 2023.
# 2022 = 198050,802
# 2023 = 202209,8688

#      Point Forecast    Lo 95    Hi 95
# 2022       199679.1 186400.7 212957.6
# 2023       205000.2 185200.5 224799.9

# Vecteurs de valeurs réelles et prédites
valeurs_reelles <- c(198050.802, 202209.8688)
valeurs_predites <- c(199679.1, 205000.2)

# Calcul des erreurs
erreurs <- valeurs_reelles - valeurs_predites

# Affichage des erreurs
print(erreurs) # -1628.298 -2790.331

# Calcul de la marge d'erreur
marge_erreur <- mean(abs(erreurs))
print(marge_erreur) # 2209.315

# Calcul des erreurs en pourcentage
erreurs_en_pourcentage <- ((valeurs_predites - valeurs_reelles) / valeurs_reelles) * 100

# Affichage des erreurs en pourcentage
print(erreurs_en_pourcentage) # 0.8221618 1.3799184

############################################################################### 
### L'estimation et la prévision avec l'année 2022 en plus
###############################################################################
# Créer une série chronologique avec les valeurs jusqu'en 2021
Petrole1 <- window(Petrole, start = c(1971), end = c(2021))
# Ajouter une valeur pour l'année 2022
nouvelles_valeurs <- c(coredata(Petrole1), 198050.802)
annees <- c(time(Petrole1), 2022)
Petrole1 <- ts(nouvelles_valeurs, start = min(annees), end = max(annees), frequency = 1)

df1 <- data.frame(Year = time(Petrole1), Value = as.numeric(Petrole1))

############################################################################### 
### Chronogramme
############################################################################### 
# Tracer le graphique avec ggplot2
ggplot(df1, aes(x = Year, y = Value)) +
  geom_line() +
  labs(x = 'Année', y = 'Production de pétrole non raffiné', title = 'Production de pétrole non raffiné de 1970 à 2022 au Canada') +
  # Ajouter la droite de régression
  geom_smooth(method = 'lm', formula = y ~ x, color = 'red', se = FALSE) 

############################################################################### 
### ACF et PACF
############################################################################### 
opp <- par(mfrow = c(1, 2))
Acf(Petrole1, lag = 30)
Pacf (Petrole1, lag = 30)
par(opp)

############################################################################### 
### Différentiation à l'ordre 2
############################################################################### 
ddPetrole1=diff(Petrole1, 2)

ddf1 <- data.frame(Year = time(ddPetrole1), Value = as.numeric(ddPetrole1))

############################################################################### 
### Chronogramme
############################################################################### 
# Tracer le graphique avec ggplot2
ggplot(ddf1, aes(x = Year, y = Value)) +
  geom_line() +
  labs(x = 'Année', y = 'Production de pétrole non raffiné', title = "Chronogramme de la série différenciée à l'ordre 2") +
  # Ajouter la droite de régression
  geom_smooth(method = 'lm', formula = y ~ x, color = 'red', se = FALSE) 

############################################################################### 
### ACF et PACF
############################################################################### 
opp <- par(mfrow = c(1, 2))
Acf(ddPetrole1, lag = 30)
Pacf (ddPetrole1, lag = 30)
par(opp)

############################################################################### 
### L'EACF
############################################################################### 
eacf(ddPetrole1) # p = O et q = 1 sachant que d = 2

############################################################################### 
### Modèle testé 1
############################################################################### 
reg1 = Arima(Petrole1, order = c(0, 2, 1), include.constant = T)
coeftest(reg1) # Le coefficient ma1 est significatif.
residu<-reg1$res
jarque.bera.test(residu) # p-value = 0.06052
# p-value > 0.05 : On accepte H0, les aléas sont gaussiens
t.test(residu) # p-value = 0.7814
# p-value > 0.05 : On accepte H0, les aléas du ARIMA(0,2,1) ont une espérance nulle.
residu=(residu-mean(residu))/sd(residu)
K<-12
tmp<-rep(0,K)
temp<-rep(0,K)
for(i in 1:K){
  tmp[i]<-Box.test(residu,lag=i,type="Ljung-Box")$p.value
  temp[i]<-ArchTest(residu, lag=i)$p.value
}
tmp
temp
# Toutes les p-values sont > à 0.05 donc on ne rejette aucune des deux hypothèses 
# nulles pour tout ordre de retard inférieur ou égal à 12. 
# Donc les aléas sont des BB.

reg1$bic # 1031.437

prev<-forecast(reg1,h=3,level=0.95)
prev
plot(prev, main = "Prévisions de la production annuelle de pétrole pour 2023 à 2025")
lines(Petrole1,col=2)

############################################################################### 
### Modèle testé 2
############################################################################### 
# modèle de la série Petrole
reg2 = Arima(Petrole1, order = c(1, 2, 1), include.constant = T)
coeftest(reg2) 
# Le coefficient ar1 n'est pas significatif. Si on estime pas ce coefficient,
# cela revient au modèle 1.

# Nous retrouvons la même conclusion pour les ARIMA(2, 2, 1), ARIMA(3, 2, 1),
# ARIMA(4, 2, 1), ARIMA(5, 2, 1) et ARIMA(6, 2, 1)

############################################################################### 
### Modèle testé 3
############################################################################### 
# modèle de la série Petrole
reg3 = Arima(Petrole1, order = c(1, 2, 2), include.constant = T)
coeftest(reg3) # Le coefficient ma1 n'est pas significatif.
reg3 = Arima(Petrole1, order = c(1, 2, 2), include.constant = T,
             fixed = c(NA, 0, NA))
coeftest(reg3) # Tous les coefficients sont significatifs
residu<-reg3$res
jarque.bera.test(residu) # p-value = 0.0757
# p-value > 0.05 : On accepte H0, les aléas sont gaussiens
t.test(residu) # p-value = 0.7539
# p-value > 0.05 : On accepte H0, les aléas du ARIMA(1,2,2) ont une espérance nulle.
residu=(residu-mean(residu))/sd(residu)
K<-12
tmp<-rep(0,K)
temp<-rep(0,K)
for(i in 1:K){
  tmp[i]<-Box.test(residu,lag=i,type="Ljung-Box")$p.value
  temp[i]<-ArchTest(residu, lag=i)$p.value
}
tmp
temp
# Toutes les p-values sont > à 0.05 donc on ne rejette aucune des deux hypothèses 
# nulles pour tout ordre de retard inférieur ou égal à 12. 
# Donc les aléas sont des BB.

reg3$bic # 1035.21
# Le bic est supérieur à celui du premier modèle donc nous conservons le modèle 1

############################################################################### 
### Modèle testé 4
############################################################################### 
# modèle de la série Petrole
reg4 = Arima(Petrole1, order = c(1, 2, 3), include.constant = T)
coeftest(reg4) # Les coefficients ma1 et ma3 ne sont pas significatifs.
# Si on estime pas ces coefficients, cela revient au modèle 3.

# Nous retrouvons la même conclusion pour les ARIMA(1, 2, 4) et ARIMA(1, 2, 5)

############################################################################### 
### Modèle testé 5
############################################################################### 
# modèle de la série Petrole
reg5 = Arima(Petrole1, order = c(2, 2, 3), include.constant = T)
coeftest(reg5) # Tous les coefficients sont significatifs.
residu<-reg5$res
jarque.bera.test(residu) # p-value = 0.2473
# p-value > 0.05 : On accepte H0, les aléas sont gaussiens
t.test(residu) # p-value = 0.7149
# p-value > 0.05 : On accepte H0, les aléas du ARIMA(2,2,3) ont une espérance nulle.
residu=(residu-mean(residu))/sd(residu)
K<-12
tmp<-rep(0,K)
temp<-rep(0,K)
for(i in 1:K){
  tmp[i]<-Box.test(residu,lag=i,type="Ljung-Box")$p.value
  temp[i]<-ArchTest(residu, lag=i)$p.value
}
tmp
temp
# Toutes les p-values sont > à 0.05 donc on ne rejette aucune des deux hypothèses 
# nulles pour tout ordre de retard inférieur ou égal à 12. 
# Donc les aléas sont des BB.

reg5$bic # 1043.14
# Le bic est supérieur à celui du premier modèle donc nous conservons le modèle 1

############################################################################### 
### Modèle testé 6
############################################################################### 
# modèle de la série Petrole
reg6 = Arima(Petrole1, order = c(3, 2, 6), include.constant = T)
coeftest(reg6) # Seulement les coefficents associés à ar2, ar3 et ma4 sont
# significatifs
reg6 = Arima(Petrole1, order = c(3, 2, 6), include.constant = T,
             fixed = c(0, NA, NA, 0, 0, 0, NA, 0, 0))
coeftest(reg6) # Aucun coefficient n'est significatif, on abandonne ce modèle.

############################################################################### 
### Modèle testé 7
############################################################################### 
# modèle de la série Petrole
reg7 = Arima(Petrole1, order = c(3, 2, 2), include.constant = T)
coeftest(reg7) # Aucun coefficient n'est significatif, on abandonne ce modèle.

############################################################################### 
### Modèle testé 8
############################################################################### 
# modèle de la série Petrole
reg8 = Arima(Petrole1, order = c(1, 2, 0), include.constant = T)
coeftest(reg8) # Le coefficient ar1 est significatif.
residu<-reg8$res
jarque.bera.test(residu) # p-value = 0.09825
# p-value > 0.05 : On accepte H0, les aléas sont gaussiens
t.test(residu) # p-value = 0.9197
# p-value > 0.05 : On accepte H0, les aléas du ARIMA(1,2,0) ont une espérance nulle.
residu=(residu-mean(residu))/sd(residu)
K<-12
tmp<-rep(0,K)
temp<-rep(0,K)
for(i in 1:K){
  tmp[i]<-Box.test(residu,lag=i,type="Ljung-Box")$p.value
  temp[i]<-ArchTest(residu, lag=i)$p.value
}
tmp
temp
# Pas toutes les p-values sont > à 0.05 donc il y a des effets ARCH, ce modèle
# n'est donc pas adéquat.

reg8$bic # 1044.62
# Le bic est supérieur à celui du premier modèle donc nous conservons le modèle 1

############################################################################### 
### Modèle testé 9
############################################################################### 
# modèle de la série Petrole
reg9 = Arima(Petrole1, order = c(1, 2, 6), include.constant = T)
coeftest(reg9) # Seul le coefficient ma2 est significatif.
reg9 = Arima(Petrole1, order = c(1, 2, 6), include.constant = T,
             fixed = c(0, 0, NA, 0, 0, 0, 0))
coeftest(reg9) # Le coefficient ma2 n'est plus significatif, on abandonne ce modèle.

###
#### CONCLUSION : le meilleur modèle est le ARIMA(0,2,1)
###