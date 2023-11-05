
#########################################################################

# R code to analyse mate choice in zebrafish
# Smith & Spence

#########################################################################

#Import the txt file into a dataframe
dr <- read.table(file = "DRdata.txt", 
               header = TRUE, 
                  dec = ".", 
     stringsAsFactors = TRUE)

str(dr, vec.len=3)

#########################################
# Load packages
library(arm)
library(car)
library(ggplot2)
library(lattice)
library(lawstat)
library(outliers)
library(tidyverse)
library(lme4)
library(car)
library(lmerTest)
library(MuMIn)
library(performance)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(gridExtra)
library(GGally)
library(plyr)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(lubridate)
library(aweek)
library(zoo)
library(rio)
library(lmtest)
library(AICcmodavg)

# Install the latest stable version of INLA:
install.packages("INLA", repos=c(getOption("repos"), 
                                 INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# Also install brinla
devtools::install_github("julianfaraway/brinla")

# And install inlatools:
ip <- rownames(installed.packages())
if (!"remotes" %in% ip) {
  install.packages("remotes")}
if (!"INLA" %in% ip) {
  install.packages(
    "INLA", 
    repos = c(getOption("repos"), "https://inla.r-inla-download.org/R/stable"))}
remotes::install_github("inbo/inlatools")

library(INLA)
library(brinla)
library(inlatools)

# For more details on installing and running the INLA package, see: https://www.r-inla.org/home

#=======================================
#Data exploration

#Any NAs
colSums(is.na(dr))
# missing length data for some males
# not a problem as these data are not used
# since they are a component of male identity

# Make factors
dr$fTrial <- as.factor(dr$trial)
dr$fDay   <- as.factor(dr$day)
dr$fRep   <- as.factor(dr$rep)
dr$fMale  <- as.factor(dr$male)
dr$fFem   <- as.factor(dr$female)
dr$fCross <- as.factor(dr$cross)
dr$fTank  <- as.factor(dr$tank)

# Define preferred figure format
My_theme <- theme(panel.background = element_blank(),  
                  panel.border = element_rect(fill = NA, size = 1),  
                  strip.background = element_rect(fill = "white",  
                                                 color = "white"), 
                  text = element_text(size = 12),  
                  panel.grid.major = element_line(colour = "white"),  
                  panel.grid.minor = element_line(colour = "white")) 

# A function for dotplots
multi_dotplot <- function(filename, Xvar, Yvar){
  filename %>%
    ggplot(aes(x = {{Xvar}})) +
    geom_point(aes(y = {{Yvar}})) +
    theme_bw() +
    coord_flip() +
    labs(x = "Order of Data")}

#=======================================

# 1. Outliers in response and independent variables

#Order data
dr <- dr %>%
  mutate(order = seq(1:nrow(dr)))

#Select continuous variables to plot
p1 <- multi_dotplot(dr, order, eggs)
p2 <- multi_dotplot(dr, order, femlen)

#Plot as a grid
grid.arrange(p1, p2, nrow = 1)
# Nothing problematic

#=======================================

# 2. Distribution of the response variable

# Frequency polygon plot for egg number
dr %>% ggplot(aes(eggs)) +
  geom_freqpoly(bins = 15) +
  labs(x = "Eggs laid", y = "Frequency") +
  My_theme +
  theme(panel.border = element_rect(colour = "black", 
                                    fill=NA, size = 1))
#Looks like zero inflation and overdispersion

#=======================================

# 3. Balance of categorical variables

table(dr$fRep)
table(dr$fTrial)
table(dr$fDay)
table(dr$fMale)
table(dr$fFem)
table(dr$fCross)
length(unique(dr$fCross))
# 96 unique crosses

table(dr$fTank)
table(dr$fMale, dr$fTank)
table(dr$fFem, dr$fTank)
# Can't include Tank as a term - it is a random 'male' effect

######################################

# 4. An excess of zeros in the response variable

# What is the percentage of zeros in the response variable

round(sum(dr$eggs == 0) * 100 / nrow(dr),0)
# 38% - a high proportion of zeros

######################################

# 5. Multicollinearity among covariates

######################################

# 6. Relationships among response and independent variables

ggplot(dr, aes(x = rank, y = (eggs))) +
  geom_jitter(shape = 16, size = 2.5, alpha = 0.3, height = 0.25, width = 0.25) +
  geom_smooth(formula = y ~ x, method = 'lm', 
              colour = 'red', se = FALSE, size = 1.5) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, size = 1)) +
  theme(strip.background = element_rect(fill = "white", 
                                        color = "white", size = 1)) +
  theme(text = element_text(size=13)) +
  xlab("Male rank") + ylab("Eggs")

ggplot(dr, aes(x = femlen, y = (eggs))) +
  geom_jitter(shape = 16, size = 2.5, alpha = 0.3, height = 1, width = 0.1) +
  geom_smooth(formula = y ~ x, method = 'lm', 
              colour = 'red', se = FALSE, size = 1.5) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, size = 1)) +
  theme(strip.background = element_rect(fill = "white", 
                                        color = "white", size = 1)) +
  theme(text = element_text(size=13)) +
  xlab("Female length") + ylab("Eggs")

# Female
ggplot(dr, aes(x = fFem, y = (eggs))) +
  geom_boxplot() +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, size = 1)) +
  theme(strip.background = element_rect(fill = "white", 
                                        color = "white", size = 1)) +
  theme(text = element_text(size=13)) +
  xlab("Female number") + ylab("Eggs")

# Male
ggplot(dr, aes(x = fMale, y = (eggs))) +
  geom_boxplot() +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, size = 1)) +
  theme(strip.background = element_rect(fill = "white", 
                                        color = "white", size = 1)) +
  theme(text = element_text(size=13)) +
  xlab("Male number") + ylab("Eggs")

# Cross
ggplot(dr, aes(x = fCross, y = (eggs))) +
  geom_boxplot() +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, size = 1)) +
  theme(strip.background = element_rect(fill = "white", 
                                        color = "white", size = 1)) +
  theme(text = element_text(size=13)) +
  xlab("Cross number") + ylab("Eggs")

#=======================================

# 7. Independence of the response variable?

dim(dr)
# 288 rows of data

length(unique(dr$fCross))
# 96 crosses

length(unique(dr$fTank))
# 24 tanks

#####################################

# FIT MODEL

# Nested design (not fully crossed)
poisson1 <- glmmTMB(eggs ~ day + 
                           (1|fFem) +
                           (1|fMale) +
                           (1|fCross),
                           family = "poisson"(link = "log"),
                           ziformula=~0,
                           data = dr)

check_overdispersion(poisson1)   # 90 - overdispersed
Res1 <- simulateResiduals(fittedModel = poisson1, plot = F)
plotResiduals(Res1)
plotQQunif(Res1)
testZeroInflation(Res1) # zero inflation
# zero inflation and overdispersion

# Deal with overdispersion with a negative binomial (linear) model
nbinom1 <- glmmTMB(eggs ~ day +
                          (1|fCross) + 
                          (1|fFem) + 
                          (1|fMale),
                          family = nbinom1(link = "log"),
                          ziformula=~0,
                          data = dr)

# Simulate data using model parameters
Res2 <- simulateResiduals(fittedModel = nbinom1, plot = F)

# Use simulated data to test zero-inflation in both (overdispersed) models
par(mfrow=c(1,1), mar=c(5,5,4,4), cex.lab = 1)
testZeroInflation(Res2)
# Negative binomial can handle this number of zeros

plotQQunif(Res2) # Still problem with dispersion

# Alternatively, deal with zeros with a zero-inflated negative binomial (linear) model
zinb <- glmmTMB(eggs ~ day + 
                       (1|fCross) +
                       (1|fFem) +
                       (1|fMale),
                       family = nbinom1(link = "log"),
                       ziformula=~ day,
                       data = dr)
summary(zinb)

# Simulate data using model parameters
SimZINB <- simulateResiduals(fittedModel = zinb, plot = F)
testZeroInflation(SimZINB)# No problem with zero-inflation
plotQQunif(SimZINB)# Still problem with dispersion

# Compare  models with AIC
models <- list(poisson1, nbinom1, zinb)
aictab(cand.set = models)
# zinb most probable
# But variance explained by random terms is unusually low
# An option is to regularise the model fit with priors

#===================================================
# Fit Bayesian model
I1.def <- inla(eggs ~ day +
                      f(fFem, model = "iid") +
                      f(fMale, model = "iid") +
                      f(fCross, model = "iid"),
                      control.compute = list(dic = TRUE, waic = TRUE),    
                      control.predictor = list(compute = TRUE),  
                      family = "poisson", 
                      data = dr)

#Posterior mean values and 95% CI
Beta1.def <- I1.def$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta1.def, digits = 2)  
# Identical output to frequentist model

bri.fixed.plot(I1.def)
bri.random.plot(I1.def)
bri.hyperpar.plot(I1.def)

# Hyper parameters
summary(I1.def)
HyperPar.I1.def <- bri.hyperpar.summary(I1.def)
round(HyperPar.I1.def, digits = 3)
# Random terms for male and female look suspicious (like the frequentist model)

# Add PC priors to Bayesian model
# Derive PC prior
sdres <- sd(dr$eggs)
round(sdres,0)
pcprior <- list(prec = list(prior="pc.prec", param = c(3*sdres,0.01)))

# Run model with PC priors
I1.pc <- inla(eggs ~ day +
                     f(fFem, model = "iid", hyper = pcprior) +
                     f(fMale, model = "iid", hyper = pcprior) +
                     f(fCross, model = "iid", hyper = pcprior),
                     control.compute = list(dic = TRUE, waic = TRUE),    
                     control.predictor = list(compute = TRUE),  
                     family = "poisson", 
                     data = dr)

#Posterior mean values and 95% CI
Beta1.pc <- I1.pc$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(Beta1.pc, digits = 2)  

I1.pc$summary.fixed
bri.hyperpar.summary(I1.pc)
bri.hyperpar.plot(I1.pc)
# Random terms look sensible

overdisp_post <- inla.tmarginal(fun=function(x) 1/x, marg=I1.pc$marginals.hyperpar[[1]])
round(inla.emarginal(fun=function(x) x, marg=overdisp_post), 4)
round(inla.qmarginal(c(0.025, 0.975), overdisp_post), 4)
# Model may be underdispersed

bri.random.plot(I1.pc)

#=============================

# Check for overdispersion 
I1.sim <- inla(eggs ~ day +
                f(fFem, model = "iid", hyper = pcprior) +
                f(fMale, model = "iid", hyper = pcprior) +
                f(fCross, model = "iid", hyper = pcprior),
                control.compute = list(config = TRUE),     #Allow for simulation    
                family= "poisson",
                data = dr)


# Simulate 1000 sets of regression parameters
NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = I1.sim)

#1.
nd <- length(rownames(SimData[[1]]$latent))
LastBeta <- rownames(SimData[[1]]$latent)[nd]

#2.
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }
Last2Character <- substrRight(LastBeta, 2)

#3. 
BetasInModel <- rownames(I1.sim$summary.fixed) #<--Change the name 'I1.sim' for your model
MyID         <- function(x){ which(rownames(SimData[[1]]$latent) == x) }
BetasInINLA <- BetasInModel
if (Last2Character == ":1") { BetasInINLA  <- paste(BetasInModel, ":1", sep ="") }

#4.
BetaRows <- lapply(BetasInINLA, MyID)
BetaRows <- as.numeric(BetaRows)

# Random effects
MyGrep <- function(x, SimulatedData){ 
  names.SimData <- attributes(SimulatedData[[1]]$latent)$dimnames[[1]]
  names.SimData[grep(x, names.SimData)] 
}

# Female
names.SimData <- attributes(SimData[[1]]$latent)$dimnames[[1]] 
FemNames   <- names.SimData[grep("fFem", names.SimData)]
FemNames

# Determine on which rows the random effects are
FemRows <- lapply(FemNames, MyID)
FemRows <- as.numeric(FemRows)
FemRows

# Male
MaleNames <- names.SimData[grep("fMale", names.SimData)]

# Determine on which rows the random effects are
MaleRows <- lapply(MaleNames, MyID)
MaleRows <- as.numeric(MaleRows)

# Cross
CrossNames <- names.SimData[grep("fCross", names.SimData)]

# Determine on which rows the random effects are
CrossRows <- lapply(CrossNames, MyID)
CrossRows <- as.numeric(CrossRows)

# Extract betas and random effects and calculate
# the fitted values and simulate count data
N    <- nrow(dr)
Ysim <- matrix(nrow = N, ncol = NSim)
mu.i <- matrix(nrow = N, ncol = NSim)
X    <- model.matrix(~ day,
                     data = dr)
X   <- as.matrix(X)

FemID <- as.numeric(as.factor(dr$fFem))
MaleID   <- as.numeric(as.factor(dr$fMale))
CrossID    <- as.numeric(as.factor(dr$fCross))

for (i in 1: NSim){
  Betas     <- SimData[[i]]$latent[BetaRows] 
  Fem_i     <- SimData[[i]]$latent[FemRows] 
  Male_ij   <- SimData[[i]]$latent[MaleRows]   
  Cross_ijk <- SimData[[i]]$latent[CrossRows] 
  
  eta <- X %*% Betas + Fem_i[FemID] + Male_ij[MaleID] + Cross_ijk[CrossID] 
  mu.i[,i] <- exp(eta)                        
  Ysim[,i] <- rpois(n = N, lambda = mu.i[,i]) 
}

# First 6 rows of the first 10 simulated data sets
head(Ysim[,1:10])

# Calculate the dispersion statistic in each of the 1,000
# data sets.
Disp <- vector(length = NSim)
N    <- nrow(dr)
Np   <- length(rownames(I1.sim$summary.fixed)) + 1 + 1 + 1 #(3 sigmas)
for(i in 1:NSim){
  ei <- (Ysim[,i] - mu.i[,i]) / sqrt(mu.i[,i])
  Disp[i] <- sum(ei^2) / (N - Np)
}

#Plot this as a table
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(Disp, 
     xlab = "Dispersion statistic",
     ylab = "Frequency",
     xlim = c(0, 100),
     main = "Simulation results")
points(x = Dispersion, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
# The red dot is the dispersion in the original data set.
# Poisson model does not handle overdispersion in the data

## Plot Pearson residuals vs fitted values
mu1 <- I1.pc$summary.fitted.values[,"mean"]
E1  <- (dr$eggs - mu1) / sqrt(mu1)

par(mfrow = c(1,1), mar = c(5,5,3,3), cex.lab = 1)
plot(x = mu1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline (h = 0, lty = 2)
# Not bad

par(mfrow = c(1,1), mar = c(5,5,3,3), cex.lab = 1)
plot(x = dr$day,
     y = E1,
     xlab = "Day",
     ylab = "Pearson residuals")
abline (h = 0, lty = 2)
# Not bad

# Zero inflation
N <- nrow(dr)
ZerosInData <- 100 * sum(dr$eggs == 0) / N
ZerosInData

# In the simulated data:
Zeros <- vector(length = NSim)
for(i in 1:NSim){
  Zeros[i] <- 100 * sum(Ysim[,i] == 0) / N
}

#Plot these as a table
par(mar = c(5,5,2,2), cex.lab = 1.5)
hist(Zeros, 
     xlim = c(0, 40),
     xlab = "Percentage of zeros",
     ylab = "Frequency",
     main = "Simulation results")
points(x = ZerosInData, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
# The red dot is the percentage of zeros in the original data set.
# Problem with zero inflation


# ============================================
# Try zero-adjusted Poisson GLMM
# ============================================

# Create binomial variable - eggs.01
dr$eggs.01 <- ifelse(test = dr$eggs >0, yes = 1, no = 0)

# Fit a Bernoulli model to binomial data
bern <- inla(eggs.01 ~ day +
                         f(fFem, model = "iid", hyper = pcprior) +
                         f(fMale, model = "iid", hyper = pcprior) +
                         f(fCross, model = "iid", hyper = pcprior),
                         control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),    
                         control.predictor = list(compute = TRUE),  
                         family = "binomial", 
                         data = dr)

#Posterior mean values and 95% CI
Betabern <- bern$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Betabern, digits = 2)  
bri.fixed.plot(bern, together = FALSE)

# Random terms
round(bri.hyperpar.summary(bern),2)
bri.hyperpar.plot(bern, together = FALSE)

# Remove random terms sequentially and compare WAIC
bern.F <- inla(eggs.01 ~ day +
                         #f(fFem, model = "iid", hyper = pcprior) +
                         f(fMale, model = "iid", hyper = pcprior) +
                         f(fCross, model = "iid", hyper = pcprior),
                         control.compute = list(waic = TRUE),    
                         control.predictor = list(compute = TRUE),  
                         family = "binomial", 
                         data = dr)

bern.M <- inla(eggs.01 ~ day +
                         f(fFem, model = "iid", hyper = pcprior) +
                         #f(fMale, model = "iid", hyper = pcprior) +
                         f(fCross, model = "iid", hyper = pcprior),
                         control.compute = list(waic = TRUE),    
                         control.predictor = list(compute = TRUE),  
                         family = "binomial", 
                         data = dr)

bern.C <- inla(eggs.01 ~ day +
                         f(fFem, model = "iid", hyper = pcprior) +
                         f(fMale, model = "iid", hyper = pcprior),
                         #f(fCross, model = "iid", hyper = pcprior),
                         control.compute = list(waic = TRUE),    
                         control.predictor = list(compute = TRUE),  
                         family = "binomial", 
                         data = dr)

bern.FM <- inla(eggs.01 ~ day +
                          #f(fFem, model = "iid", hyper = pcprior) +
                          #f(fMale, model = "iid", hyper = pcprior) +
                          f(fCross, model = "iid", hyper = pcprior),
                          control.compute = list(waic = TRUE),    
                          control.predictor = list(compute = TRUE),  
                          family = "binomial", 
                          data = dr)

FemRE     <- c(bern$waic$waic - bern.F$waic$waic)
MaleRE    <- c(bern$waic$waic - bern.M$waic$waic) 
CrossRE   <- c(bern$waic$waic - bern.C$waic$waic)
FemMaleRE <- c(bern$waic$waic - bern.FM$waic$waic)

round(FemRE,     digits = 1) #additive female effects not important
round(MaleRE,    digits = 1) #additive male effects not important
round(CrossRE,   digits = 1) #non-additive effects important
round(FemMaleRE, digits = 1) #additive effects not important

# within 2 units is considered about the same, 
# within 7 units similar model fit, 
# >7 different

# Create plot

MyData <- expand.grid(
  day = seq(from = min(dr$day), 
              to = max(dr$day), length = 50))
Xpred <- model.matrix(~ 1 + day,
                      data = MyData)
head(Xpred)
Xpred <- as.data.frame(Xpred)
lcb <- inla.make.lincombs(Xpred)
I1.Pred <- inla(eggs.01 ~ day +
                          f(fFem, model = "iid", hyper = pcprior) +
                          f(fMale, model = "iid", hyper = pcprior) +
                          f(fCross, model = "iid", hyper = pcprior),
                          control.predictor = list(link = 1,
                                                compute = TRUE),
                          data = dr, 
                          lincomb = lcb,
                          family = "binomial")

# Get the marginal distributions:
Pred.marg <- I1.Pred$marginals.lincomb.derived

# The output above is on the logit-scale. 
# This function converts x into exp(x) / (1 + exp(x))
MyLogit <- function(x) {exp(x)/(1+exp(x))}

MyData$mu <- unlist( 
  lapply(
    Pred.marg,
    function(x) inla.emarginal(MyLogit,x)))

MyData$selo <- unlist( 
  lapply(
    Pred.marg,
    function(x) 
      inla.qmarginal(c(0.025), 
                     inla.tmarginal(MyLogit, x))))

MyData$seup <- unlist( 
  lapply(
    Pred.marg,
    function(x) 
      inla.qmarginal(c(0.975), 
                     inla.tmarginal(MyLogit, x))))

plotA <- ggplot() +
  geom_jitter(data = dr, aes(y = eggs.01, x = day),
              shape = 16, size = 4, height = 0.015,
              width = 0.45, alpha = 0.6) +
  xlab("Day of experiment") + ylab("Probability of spawning") + 
  My_theme + xlim(-0.5,12.5) +
  geom_line(data = MyData, aes(x = day, y = mu)) +
  geom_ribbon(data = MyData, aes(x = day, ymax = seup, ymin = selo),
            alpha = 0.4)
plotA

# ======================================================
# Model zero-truncated data with gpoisson distribution

dr1 <- dr
dr1$eggs.pos <- ifelse(test = dr$eggs ==0, yes = NA, no = dr$eggs)
dr2 <- na.omit(dr1)

gp <- inla(eggs.pos ~ day +
                     f(fFem, model = "iid", hyper = pcprior) +
                     f(fMale, model = "iid", hyper = pcprior) +
                     f(fCross, model = "iid", hyper = pcprior),
                     control.compute = list(dic = TRUE, waic = TRUE),    
                     control.predictor = list(compute = TRUE),  
                     family = "gpoisson", 
                     data = dr2)

#Posterior mean values and 95% CI
Betagp <- gp$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")] 
print(Betagp, digits = 2)  
bri.fixed.plot(gp, together = FALSE)

bri.hyperpar.summary(gp)
bri.hyperpar.plot(gp, together = FALSE)
bri.random.plot(gp)

gp.F <- inla(eggs.pos ~ day +
                       #f(fFem, model = "iid", hyper = pcprior) +
                        f(fMale, model = "iid", hyper = pcprior) +
                        f(fCross, model = "iid", hyper = pcprior),
                        control.compute = list(waic = TRUE),    
                        control.predictor = list(compute = TRUE),  
                        family = "gpoisson", 
                        data = dr2)

gp.M <- inla(eggs.pos ~ day +
                        f(fFem, model = "iid", hyper = pcprior) +
                       #f(fMale, model = "iid", hyper = pcprior) +
                        f(fCross, model = "iid", hyper = pcprior),
                        control.compute = list(waic = TRUE),    
                        control.predictor = list(compute = TRUE),  
                        family = "gpoisson", 
                        data = dr2)

gp.C <- inla(eggs.pos ~ day +
                        f(fFem, model = "iid", hyper = pcprior) +
                        f(fMale, model = "iid", hyper = pcprior),
                       #f(fCross, model = "iid", hyper = pcprior),
                        control.compute = list(waic = TRUE),    
                        control.predictor = list(compute = TRUE),  
                        family = "gpoisson", 
                        data = dr2)

gp.FM <- inla(eggs.pos ~ day +
                        #f(fFem, model = "iid", hyper = pcprior) +
                        #f(fMale, model = "iid", hyper = pcprior) +
                         f(fCross, model = "iid", hyper = pcprior),
                         control.compute = list(waic = TRUE),    
                         control.predictor = list(compute = TRUE),  
                         family = "gpoisson", 
                         data = dr2)

# Compare WAICs
FemREgp     <- c(gp$waic$waic-gp.F$waic$waic)
MaleREgp    <- c(gp$waic$waic-gp.M$waic$waic) 
CrossREgp   <- c(gp$waic$waic-gp.C$waic$waic)
FemMaleREgp <- c(gp$waic$waic-gp.FM$waic$waic)

round(FemREgp,     digits = 1) #additive female effects not important
round(MaleREgp,    digits = 1) #additive male effects not important
round(CrossREgp,   digits = 2) #non-additive effects highly important
round(FemMaleREgp, digits = 1) #additive effects not important

# Create plot
MyData <- expand.grid(
  day = seq(from = min(dr$day), 
            to = max(dr$day), length = 50))
Xpred <- model.matrix(~ 1 + day, data = MyData)
Xpred <- as.data.frame(Xpred)
lcb <- inla.make.lincombs(Xpred)

I1.Predgp <- inla(eggs.pos ~ day +
                  f(fFem, model = "iid", hyper = pcprior) +
                  f(fMale, model = "iid", hyper = pcprior) +
                  f(fCross, model = "iid", hyper = pcprior),
                  data = dr2, 
                  lincomb = lcb,
                  control.predictor = list(link = 1,
                                         compute = TRUE),
                  family = "gpoisson")

# Get the marginal distributions:
Pred.marg <- I1.Predgp$marginals.lincomb.derived

MyData$mu <- unlist( 
  lapply(
    Pred.marg,
    function(x) inla.emarginal(exp,x)))

MyData$selo <- unlist( 
  lapply(
    Pred.marg,
    function(x) 
      inla.qmarginal(c(0.025), 
                     inla.tmarginal(exp, x))))

MyData$seup <- unlist( 
  lapply(
    Pred.marg,
    function(x) 
      inla.qmarginal(c(0.975), 
                     inla.tmarginal(exp, x))))

plotB <- ggplot() +
  geom_jitter(data = dr2, aes(y = eggs.pos, x = day),
                     shape = 16, size = 4, height = 0.01,
                     width = 0.45, alpha = 0.6) +
  xlab("Day of experiment") + ylab("Number of eggs spawned") + 
  My_theme + xlim(-0.5,12.5) +
  geom_line(data = MyData, aes(x = day, y = mu)) +
  geom_ribbon(data = MyData, aes(x = day, ymax = seup, ymin = selo),
                     alpha = 0.4)
plotB

# Combine plots
ggarrange(plotA, plotB,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

#=======================================

# Plot residuals - Bernoulli
Fit <- I1.Pred$summary.fitted.values[, "mean"]
Res <- dr$eggs.01 - Fit
ResPlot <- cbind.data.frame(Fit,Res,dr$eggs.01,dr$day)

FigA <- ggplot(ResPlot, aes(x=Fit, y=Res)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Bayesian residuals") + xlab("Fitted values") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

FigB <- ggplot(ResPlot, aes(x=dr$day, y=Res)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("day") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

ggarrange(FigA, FigB,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# Plot residuals - generalised Poisson
Fit2 <- I1.Predgp$summary.fitted.values[, "mean"]
Res2 <- dr2$eggs.pos - Fit
ResPlot2 <- cbind.data.frame(Fit,Res,dr2$eggs.pos,dr2$day)

# Plot residuals against fitted
FigA2 <- ggplot(ResPlot2, aes(x=Fit2, y=Res2)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Bayesian residuals") + xlab("Fitted values") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# And plot residuals against variables in the model
FigB2 <- ggplot(ResPlot2, aes(x=dr2$day, y=Res2)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("day") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Combine plots
ggarrange(FigA2, FigB2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# ==END==
