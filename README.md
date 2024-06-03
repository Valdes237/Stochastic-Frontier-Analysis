# Stochastic-Frontier-Analysis
Enhält Code zum SFA

#####################################################################################################
 Thema: Absch?tzung der Produktionsfunktion mittels stochastischer Grenzanalyse.           
                                                                                          
#################################### Pakage laden und lesen ##########################################

rm(list = ls())

# Pakages laden
if(!require(HDInterval)){install.packages("HDInterval")}
if(!require(sfa)){install.packages("sfa")}
if(!require(MASS)){install.packages(MASS)}
if(!require(invgamma)) {install.packages("invgamma")}
if(!require(coda)) {isntall.packages(coda)}
if(!require(BayesFactor)){install.packages("BayesFactor")}
if(!require(BayesFactor)){install.packages("statmod")}
library(HDInterval)
library(ggplot2)
library(MASS)
library(invgamma)
library(coda)
library(sfa)
library(mvtnorm)
library(BayesFactor)
library(dplyr)
library(statmod)

##################################### Vorbereitung der Daten ########################################

# Daten laden
Dairy_Farm <- read.csv2("H:/Semester 5/Masterarbeit/Daten/Spanish Dairy Prod/Spanish Dairy Farm Production.csv",
                        sep = ",", quote = "", dec = ".")

# Daten Visualisieren
View(Dairy_Farm)

# Dimension der Daten
dim(Dairy_Farm)

# Extraktion der Daten aus dem Jahr 93
Dairy_Farm_new <- filter(Dairy_Farm, Dairy_Farm$YEAR == 93)
dim(Dairy_Farm_new)

# Visualisierung des neuen Datensatz
View(Dairy_Farm_new)
attach(Dairy_Farm_new)

## Daten logarithmieren
# wir logarithmieren die relevante Daten f?r unsere Simulation.
logMilk <-log(MILK)
logLabor <- log(LABOR)
logCows <- log(COWS)
logFeed <- log(FEED)
logLand <- log(LAND)

# Neuer Datensatz mit nur die relevanten Daten
Dairy_Farm_new <- data.frame(Dairy_Farm_new, logMilk, logLabor, logCows, logFeed, logLand)

# Anzahl der Beobachtungen
n <- nrow(Dairy_Farm_new)

# Visiualisierung der Daten
View(Dairy_Farm_new)
attach(Dairy_Farm_new)



####################### Sch?tzung mittels Maximum Likelihood Estimation Methode ######################

# Durchf?hrung der SFA Funktion 
Modell_mle <- sfa( logMilk ~ logCows + logLabor + logFeed + logLand, data = Dairy_Farm_new, 
            fun = "exp", form = "production", pars = NULL, pars_mu = NULL, method = "BFGS")

summary(Modell_mle)

# Koeffizienten
coef(Modell_mle) # oder auch
print(Modell_mle)

# Betas Koeffizienten
beta_0 <- coef(Modell_mle)[[1]]
beta_1 <- coef(Modell_mle)[[2]]
beta_2 <- coef(Modell_mle)[[3]]
beta_3 <- coef(Modell_mle)[[4]]
beta_4 <- coef(Modell_mle)[[5]]


# varianz der Ineffizienz
sigma2u <- coef(Modell_mle)[[6]]

# varianz des statistischen random noise.
sigma2v <- coef(Modell_mle)[[7]]

# Anteil  der Ineffizienz Varianz zu der gesamte Varianz
I <- sqrt(coef(Modell_mle)[[6]])/sqrt(coef(Modell_mle)[[7]]) # I n?hrt sich von 0, die Varianz von 
                                                            # v dominiert die Varianz von u|e

# Total Residuum (epsilon) aus dem Modell
e <- residuals(Modell_mle) # epsilon = v - u

# Total Varianz
sigma2 <- sigma2u + sigma2v

# mu tilde   
mu_tild <- -e - (sigma2v*sqrt(sigma2u))

A <- -(-e/sqrt(sigma2v) - sqrt(sigma2v)*sqrt(sigma2u))


# Sch?tzung von u_i bedingt auf epsilon (u_hat oder E(u|e) bzw. Ineffizienz der Firmen
u_hat <- mu_tild + sqrt(sigma2v)*(dnorm(A)/pnorm(-A))

# Mittelwert der Ineffizienz
mean(u_hat)

# Produzentenspezifisch technische Effizienz
TE1_i <- exp(-u_hat) 

# Technische Effizienz in Industrie Bereich
TE1_i_strich <- mean(TE1_i) 

## Sch?tzung von u_i bedingt auf epsilon mit der Mode
mode.u <- exp(mu_tild)
mean(mode.u)

## Visualisierung der technische Effizienz
# Histogramm technische Effizienz
windows()
G <- ggplot(Dairy_Farm_new, aes(x = TE1_i)) +
  geom_histogram(aes(y= ..density..),bins = 30, color = "black", fill = "grey") +
  geom_density(alpha= 0.4, fill = "tomato") + ggtitle("Plot efficiency") + xlab("technical efficiency") + ylab("count")
G 
dev.off()

# Plot  technische Effizienz gegen Farmen
windows()
P <- ggplot(Dairy_Farm_new, aes(x = FARM, y = sort(TE1_i))) +
  geom_line() + ggtitle("Efficiency") + ylab("technical efficiency") + xlab("Farm")
P
dev.off()

## Konfidenzintervall der Ineffizienz
# Untere Grenze von u|e
LB_i <- mu_tild + pinvgauss(1 - (1 - 0.05/2)*(1-pnorm(A))) * sqrt(sigma2v) 
range(LB_i)

# Obere Grenze von u|e
UB_i <- mu_tild + pinvgauss(1 - (0.05/2) * (1-pnorm(A))) * sqrt(sigma2v)
range(UB_i)

# Untere und obere Grenze von exp(-u|e).
# Wegen der monotone Umwandlung von u_i die unter und obere Grenze von u_i|e wird zu obere und untere
# Grenze von TE1_i.
Li <- exp(-UB_i)
range(Li)

Ui <- exp(-LB_i)
range(Ui)


## Plot TE1_i mit Konfindenzintervall
Predframe <- with(Dairy_Farm_new, data.frame(TE1_i, Li, Ui))
windows()
W <- ggplot(Dairy_Farm_new, aes(x = FARM, y = TE1_i)) +
  geom_point() + geom_line(data = Predframe) +
  ggtitle("technical efficiency") + ylab("technical efficiency") + xlab("FARM")
W <- W + geom_ribbon(data = Predframe, aes(ymin = Li, ymax = Ui), 
                     alpha = 0.2, color ="red", fill = "blue") + scale_y_continuous(breaks=c(seq(0, max(Ui), 0.3)))


W
dev.off()

## Rangkorrelationskoffizient nach Spearman.
cor(Dairy_Farm_new[,c("logMilk","logCows","logFeed","logLabor","logLand")], method = "spearman")

## Test zur Heteroskedaszitit?t von Breusch-Pagan
v <- e + u_hat
v2 <- v^2
y <- as.matrix(logMilk)
X <- cbind(rep(1,length(n)), logCows, logLabor, logFeed, logLand)

# Hilfregression
HR <- lm(v2 ~ X, data = Dairy_Farm_new)
summary(HR)

LM <- n*0.00772
# LM ist chi Quadrat verteilt a 4 DF.
# LM < chi quadrat(4, 1-5%) --> H0 nicht abgelehnt: es liegt Homoskadastizit?t vor



########################## Sch?tzung mittels Bayesschen Sch?tzung Methode ############################

set.seed(1234)


b <- 5000 # burn-in
R <- 100000 # Iterationen nach burn-in

## Speicher Vektoren
Beta <- matrix(0, nrow = b+R, ncol = 5) # Sample Vektor f?r die Betas Koeffizienten
sigma2V <- matrix(0, nrow = b+R, ncol = 1) # Sample Vektor f?r Varianz des statistischen random Noise
sigmaU <- matrix(0, nrow = b+R, ncol = 1) # Sample Vektor f?r Ineffizienz Parameter
Ineff <- matrix(0, nrow = b+R, ncol = n) # Sample Vektor von Ineffizienz
inef <- matrix(0,n, 1) # Ineffizienz Vektor

## priors parameter
# hyperparameter f?r Beta 
Beta_0 <- c(0, 0, 0, 0, 0)
inv_B_0 <- 1 * diag(5) 

# hyperparametr f?r Sigmav quadrat
alpha_v_0 <- 1  # phi_v
beta_v_0 <- 1   # p_v

# hyperparameter f?r die Ineffizienz
to <- 0.875
lambda_0 <- 1/log(to)

# hyperparameter f?r sigmaU
alpha_u_0 <- 1
beta_u_0 <- -log(to)

## Initiale Werte
# Initiale Werte f?r Betas
beta <-  c(beta_0, beta_1, beta_2, beta_3, beta_4)

# Initiale Werte f?r die Ineffizienz 
u <- rexp(n, -lambda_0) 


## Berechnungen
y <- as.matrix(logMilk)
X <- cbind(rep(1,length(n)), logCows, logLabor, logFeed, logLand) 

# Quadrat Summe der Residuen.
SSE <- sum((y - X%*%beta + u)^2)

# Initiale Werte f?r sigmav quadrat (Prr?zision Parameter h = 1/sigma2v0)
sigma2v0 <- SSE/(n-5) 

# alpha_v_n <- alpha_v_0 + n/2  # hyperparameter for 1/sigma2v
alpha_v_n <- alpha_v_0 + n/2

# Vektor Anzahl der Farmen
r <- rep(1,n)


for (i in 1:(b+R)) {

  # sample beta/y, h, u, lambda
  B_n <- solve((1/(sigma2v0))*t(X)%*%X + solve(inv_B_0))
  beta_n <- solve((1/(sigma2v0))*t(X)%*%X + solve(inv_B_0))%*%((1/(sigma2v0))*t(X)%*%(y+u) + solve(inv_B_0)%*%Beta_0)
  beta_i <- mvrnorm(1, beta_n, B_n)
  
  # sample h/y, u, lambda, beta
  beta_v_n <- SSE/2 + beta_v_0 
  h <- rgamma(1, alpha_v_n, beta_v_n)
  
  # sample sigmau(1/lambda)/y, h, beta, u
  beta_u_n <- u - log(to)
  inv_lambda <- rinvgamma(1, n + 1, beta_u_n) 
  
  # sample u_i/y,h, beta,lambda mit rejection sampling  
  u_i <- - y + X%*%beta_i - (1/h)*(inv_lambda)*r
  
  # mean von u_i
  mu <- mean(u_i)
  
  for (j in 1:n) {
     a = (u_i[j]) # truncated point
     c = Inf
     lambda = a
     repeat{
     z1 <- a - 1/lambda*log(1 - runif(1)) # inverse Umwandlung der translated exponential.
      if(z1 >= a) {
        z <- lambda*exp(-lambda*(z1 - a))
      } 
      r <- exp(-0.5*(z-lambda)^2)
      w <- runif(1)
      if (w <= r) 
        break
    } 
    inef[j] <- z
    
  }
  
  # gezogene Parameter speichern
  Beta[i,] <-beta_i
  sigma2V[i] <- 1/h
  sigmaU[i] <-inv_lambda
  Ineff[i,] <- inef
  
}


# final beta, sigma und u|e 

BI_CHAIN <- as.mcmc(data.frame(beta0 = Beta[1:b,1],
                               beta1 = Beta[1:b,2],
                               beta2 = Beta[1:b,3],
                               beta3 = Beta[1:b,4],
                               beta4 = Beta[1:b,5],
                               sigm_v = sigma2V[1:b],
                               sigm_U = sigmaU[1:b],
                               Ineffi1 = Ineff[1:b, 1],
                               Ineffi2 = Ineff[1:b, 2],
                               Ineffi3 = Ineff[1:b, 3],
                               Ineffi4 = Ineff[1:b, 4],
                               Ineffi5 = Ineff[1:b, 5],
                               Ineffi242 = Ineff[1:b, 247],
                               Ineffi243 = Ineff[1:b, 246],
                               Ineffi244 = Ineff[1:b, 245],
                               Ineffi245 = Ineff[1:b, 244],
                               Ineffi246 = Ineff[1:b, 243],
                               Ineffi247 = Ineff[1:b, 242]))

CHAIN <- as.mcmc(data.frame(beta0 = Beta[(b+1):(b+R), 1],
                            beta1 = Beta[(b+1):(b+R), 2],
                            beta2 = Beta[(b+1):(b+R), 3],
                            beta3 = Beta[(b+1):(b+R), 4],
                            beta4 = Beta[(b+1):(b+R), 5],
                            sigma2_v = sigma2V[(b+1):(b+R)],
                            sigma_U = sigmaU[(b+1):(b+R)],
                            Ineffi1 = Ineff[(b+1):(b+R), 1],
                            Ineffi2 = Ineff[(b+1):(b+R), 2],
                            Ineffi3 = Ineff[(b+1):(b+R), 3],
                            Ineffi4 = Ineff[(b+1):(b+R), 4],
                            Ineffi5 = Ineff[(b+1):(b+R), 5],
                            Ineffi242 = Ineff[(b+1):(b+R), 242],
                            Ineffi243 = Ineff[(b+1):(b+R), 243],
                            Ineffi244 = Ineff[(b+1):(b+R), 244],
                            Ineffi245 = Ineff[(b+1):(b+R), 245],
                            Ineffi246 = Ineff[(b+1):(b+R), 246],
                            Ineffi247 = Ineff[(b+1):(b+R), 247]))

summary(CHAIN)

# Sch?tzer von beta, sigma2v, sigmau und u|e
beta0 = Beta[(b+1):(b+R),1]
beta1 = Beta[(b+1):(b+R),2]
beta2 = Beta[(b+1):(b+R),3]
beta3 = Beta[(b+1):(b+R),4]
beta4 = Beta[(b+1):(b+R),5]
sigma2_v = sigma2V[(b+1):(b+R)]
sigma_u = sigmaU[(b+1):(b+R)]
Ineff <- Ineff[(b+1):(b+R),]
Ineffi1 <- Ineff[,1]
Ineffi2 <- Ineff[,2]
Ineffi3 <- Ineff[,3]
Ineffi4 <- Ineff[,4]
Ineffi5 <- Ineff[,5]
Ineffi243 <- Ineff[,243]
Ineffi244 <- Ineff[,244]
Ineffi245 <- Ineff[,245]
Ineffi246 <- Ineff[,246]
Ineffi247 <- Ineff[,247]


# Analyse der Konvergenz
windows()
plot(CHAIN)
dev.off()

# Gesch?tzte Parameter
Beta0 <- mean(beta0)
Beta1 <- mean(beta1)
Beta2 <- mean(beta2)
Beta3 <- mean(beta3)
Beta4 <- mean(beta4)
Sigma2_v = mean(sigma2V[(b+1):(b+R)])
Sigma_u = mean(sigmaU[(b+1):(b+R)])
sigma2 <- Sigma_u^2 + Sigma2_v

# Relative Beitrag der Ineffizienz zu der gesamte Varianz
I <- Sigma_u^2 / Sigma2_v   # I n?hrt sich zu null : v dominiert u|e im Modell
                           

# Produzentspezifische Ineffizienz (u_i|e)
Ineffizienz <- matrix(0, n, 1)
for (i in 1:n) {
  Ineffizienz[i] <- mean(Ineff[, i])
  
}

# High Posterior Density f?r u|e
hdi(Ineffizienz)

# Mittelwert der Ineffizienz (E(u|e))
u <- mean(Ineffizienz)

# Produzent spezifische technische Effizienz (TE2_i) 
TE2_i <- exp(- Ineffizienz)

# Technische Effizienz in Industrie Bereich
TE2_i_strich <- mean(TE2_i)

## 95% High Posterior Density f?r TE2_i 
# Berechnung der untere und obere Grenze von u_i|e
LB <- matrix(0, n, 1)
UB <- matrix(0, n, 1)

for (i in 1:n) {
  LB[i] <- sort(Ineff[, i])[0.025*length(Ineff[,i])]
  UB[i] <- sort(Ineff[,i])[0.975*length(Ineff[,i])]
  
 }

# Berechnung der unter und obere Grenze von TE2_i
# untereGrenze von TE2_i
Li <- exp(-UB)
range(Li)

# obere Grenze von TE2_i
Ui <- exp(-LB)
range(Ui)

# Visualisierung der TE2_i mit Konfidenz Interval
windows()
Pred_datfram <- with(Dairy_Farm_new, data.frame(Ineffizienz, Li, Ui))
g <- ggplot(Dairy_Farm_new, aes(x = FARM, y = TE2_i)) +
  geom_point() + geom_line(data = Pred_datfram) +
  ggtitle("Plot technical efficiency") + ylab("technical efficiency") + xlab("FARM")
g <- g + geom_ribbon(data = Pred_datfram, aes(ymin = Li, ymax = Ui), 
                     alpha = 0.2, color ="red", fill = "blue") + scale_y_continuous(breaks=c(seq(0, max(Ui), 0.3)))


g
dev.off()


## Test zur Hypothese mit Bayes faktor
# Test zur Hypothese f?r H0: P(sigma2u) = 0 VS H1: P(sigma2u) > 0 
BETAS <- c(Beta0, Beta1, Beta2, Beta3, Beta4)
y1 <- X%*%BETAS + rnorm(n,0,sqrt(Sigma2_v)) - Ineffizienz # F?r das erste Modell
y0 <- X%*%BETAS + rnorm(n,0,sqrt(Sigma2_v))              # F?r das null Modell

all.equal(y1, y0)

# Daten f?r das erste und null Modell
DatFr <- data.frame(y1,y0)

# wir gehen davon aus, dass die Voraussetzung f?r den t-Test erf?llt sind
with(DatFr, t.test(DatFr$y0, DatFr$y1, alternative = "greater", paired = TRUE))

# t-Wert
t_wert <- 3.8384

# Anzahl der gesamten Farmen
n <- nrow(Dairy_Farm_new)

# hyperparameter Cauchy-Prior
scale <-  sqrt(2)/2 # default value in Bayes Faktor
location <- 0 # default value
Sig2u <- Sigma_u^2

# Funktion der Likelihood f?r das erste Modell
funk <- function(Sig2u){
  dt(t_wert, df = n -1, ncp = Sig2u*sqrt(n))*dcauchy(Sig2u, location = 0, scale = sqrt(2)/2)
}

# Likelihood f?r das erste und null Modell
L1 <- integrate(funk, 0, Inf)[[1]]/0.5
L0 <- dt(t_wert, df = n-1)

# Bayes Faktor
BF10 <- L1/L0
BF10

# Alternative
# Bayes Faktor mit Package
B <- ttestBF(x = y0, y = y1, paired = TRUE, rscale = sqrt(2)/2, nullInterval = c(0, Inf))
summary(B)

