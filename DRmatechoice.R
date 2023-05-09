
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

#'data.frame':	288 obs. of  9 variables:
# $ trial  : int  1 1 ...
# $ day    : int  1 1 ...
# $ rep    : int  1 1 ...
# $ male   : int  1 2 ...
# $ malelen: int  35 36 ...
# $ female : int  1 3 ...
# $ femlen : int  39 40 ...
# $ eggs   : int  0 0 ...
# $ tank   : int  5 3 ...

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
My_theme <- theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.ticks.x=element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, size = 1),
                  strip.background = element_rect(fill = "white", 
                                                  color = "white", size = 1),
                  text = element_text(size = 14),
                  panel.grid.major = element_line(colour = "white", size = 0.1),
                  panel.grid.minor = element_line(colour = "white", size = 0.1))

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

# 2. Normality and homogeneity of the response variable

# Frequency polygon plot for catch
dr %>% ggplot(aes(eggs)) +
  geom_freqpoly(bins = 15) +
  labs(x = "Eggs laid", y = "Frequency") +
  My_theme +
  theme(panel.border = element_rect(colour = "black", 
                                    fill=NA, size = 1))

#Shapiro-Wilk test for deviation from normality
shapiro.test(dr$eggs)

# Shapiro-Wilk normality test
# 
# data:  dr$eggs
# W = 0.78262, p-value <0.001

# Data positively skewed, with a lot of zeros

#=======================================

# 3. Balance of categorical variables

# Among bricks
table(dr$fRep)
#  1  2  3  4  5  6 
# 48 48 48 48 48 48

table(dr$fTrial)
# 1  2  3 
# 96 96 96 

table(dr$fDay)
# 1  2  3  4  5  6  7  8  9 10 11 12 
# 24 24 24 24 24 24 24 24 24 24 24 24

table(dr$fMale)
# mA mB mC mD mE mF mG mH mI mJ mK mL mM mN mO mP mQ mR mS mT mU mV mW mX 
# 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 

table(dr$fFem)
# fA fB fC fD fE fF fG fH fI fJ fK fL fM fN fO fP fQ fR fS fT fU fV fW fX 
# 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12  

table(dr$fCross)
length(unique(dr$fCross))
# 96 unique crosses

table(dr$fTank)
# 1  2  3  5  6  8  9 10 11 12 14 15 16 18 19 20 21 22 23 25 26 27 29 30 
# 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 

table(dr$fMale, dr$fTank)
# Male and tank collinear

table(dr$fFem, dr$fTank)
# And female

# Can't include Tank as a term - it is a random 'male' effect

######################################

# 4. An excess of zeros in the response variable

# What is the percentage of zeros in the response variable

round(sum(dr$eggs == 0) * 100 / nrow(dr),0)
# 38% - a high proportion of zeros in the data....

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

# Each row of data cannot be assumed independent

#####################################

# Nested design (not fully crossed - coding of male, females, crosses important)
poisson1 <- glmmTMB(eggs ~ day + 
                           (1|fFem) +
                           (1|fMale) +
                           (1|fCross),
                           family = "poisson"(link = "log"),
                           ziformula=~0,
                           data = dr)

check_overdispersion(poisson1)   # 90 - overdispersed

# Simulate data using model parameters
SimPois <- simulateResiduals(fittedModel = poisson1, plot = F)

# Use simulated data to test zero-inflation in both (overdispersed) models
par(mfrow=c(1,1), mar=c(5,5,4,4), cex.lab = 1)
testZeroInflation(SimPois)
# Too many zeros for a Poisson distribution

# Deal with overdispersion with a negative binomial (linear) model
nbinom1 <- glmmTMB(eggs ~ day +
                          (1|fCross) + 
                          (1|fFem) + 
                          (1|fMale),
                          family = nbinom1(link = "log"),
                          ziformula=~0,
                          data = dr)

check_overdispersion(nbinom1)  #0.4 Underdispersed

# Simulate data using model parameters
SimPois <- simulateResiduals(fittedModel = nbinom1, plot = F)

# Use simulated data to test zero-inflation in both (overdispersed) models
par(mfrow=c(1,1), mar=c(5,5,4,4), cex.lab = 1)
testZeroInflation(SimPois)
# Negative binomial can handle this number of zeros

# Alternatively, deal with zeros with a zero-inflated negative binomial (linear) model
zinb <- glmmTMB(eggs ~ day + 
                       (1|fCross) + 
                       (1|fFem) + 
                       (1|fMale),
                       family = nbinom1(link = "log"),
                       ziformula=~ day,
                       data = dr)

check_overdispersion(zinb)  #0.55 Still some underdispersion

# Simulate data using model parameters
SimZINB <- simulateResiduals(fittedModel = zinb, plot = F)

# Examine zero-inflation
par(mfrow=c(1,1), mar=c(5,5,4,4), cex.lab = 1)
testZeroInflation(SimZINB)
# No problem with zero-inflation

# Compare  models with AIC
models <- list(poisson1, nbinom1, zinb)
aictab(cand.set = models)
# zinb most probable

# Which random term is significant?
# Likelihood Ratio Test
zinbA <- glmmTMB(eggs ~ day +
                       (1|fCross) + 
                       (1|fMale),
                       family = nbinom1(link = "log"),
                       ziformula=~day,
                       data = dr)

zinbB <- glmmTMB(eggs ~ day +
                       (1|fCross) + 
                       (1|fFem),
                       family = nbinom1(link = "log"),
                       ziformula=~day,
                       data = dr)

zinbC <- glmmTMB(eggs ~ day +
                      (1|fFem) + 
                      (1|fMale),
                       family = nbinom1(link = "log"),
                       ziformula=~day,
                       data = dr)

zinbD <- glmmTMB(eggs ~ day +
                       (1|fCross),
                       family = nbinom1(link = "log"),
                       ziformula=~day,
                       data = dr)

# Likelihood ratio test
lrtest(zinb,zinbA) #additive female effects not important
lrtest(zinb,zinbB) #additive male effects not important
lrtest(zinb,zinbC) #non-additive effects highly important
lrtest(zinb,zinbD) #additive effects not important

###############################

# MODEL VALIDATION

#Obtain residuals and fitted values
Res <- resid(zinb)
Fit <- fitted(zinb)

# Plot residuals against fitted values
par(mfrow = c(1,2), mar = c(5,5,2,2))
plot(x = Fit,
     y = Res,
     xlab = "Fitted values",
     ylab = "Residuals",
     pch = 16, cex = 1.5)
abline(h = 0, lty = 2)

# Day
plot(x = dr$day,
     y = Res,
     xlab = "Fitted values",
     ylab = "Residuals",
     pch = 16, cex = 1.5)
abline(h = 0, lty = 2)


# Examine random effects
set_theme(base = theme_bw(),
 axis.textsize = 1)

plot_model(zinb,
           vline.color = "black",
                 title = "",
                  type = "re")

# Plot fixed effect day
My_theme <- theme(panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, size = 1),
                  strip.background = element_rect(fill = "white", 
                                                  color = "white", size = 1),
                  text = element_text(size = 14),
                  panel.grid.major = element_line(colour = "white", size = 0.1),
                  panel.grid.minor = element_line(colour = "white", size = 0.1))


plot_model(zinb,
                type = "pred",
               terms = c("day [all]"),
        show.zeroinf = T,
           show.data = T,
           pred.type = c("fe"),
               title = "",
         show.legend = F,
              jitter = 0.5,
          axis.title = c("Day of experiment","Eggs spawned")) + 
          My_theme +
         scale_x_continuous(limits = c(0, 12), 
                            breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0, 700), 
                     breaks = c(0, 100, 200, 300, 400, 500, 600, 700))

# Tabulate results
tab_model(zinbD,
          show.zeroinf = T,
             dv.labels = c("ZINB GLMM (D. rerio mate choice)"),
           string.pred = "Coeffcient",
             string.ci = "Conf. Int (95%)",
              string.p = "P-value",
               p.style = c("numeric"),
                emph.p = FALSE,
             transform = NULL
          )
