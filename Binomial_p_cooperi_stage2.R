#' Highland Statistics Ltd.
#' www.highstat.com

# Section 1: Data description----

#' continued with the code from stage 1


#'   Germination success was determined as the number of seeds germinated per 
#'   plate. 

#' The response variable is `per`, which is the success of germination
#' in a plate in a proportional data. 
#
#' There are three covariates `fun`, `temp`, and `pop`.


#* Subsection 1.3: The response variable and covariates----

#' Each row in the data represents data from a plate.

#' Response Variable:
#'   Germination success: The success proportion of seed germination.


#' Covariate:
#'  Fungi: The four fungi used in experiment. They are named pev_cerato, pev_tul,
#'  FF_cerato, and FF_tul
#'  Temperature: The two tempeatures of experiments taken placed. They are 
#'  15C and 23C
#'  Populations: The two populations where the seeds collected from. They are
#'  FF and MX

#' Dependency:
#'     We have 200 observations (plates)


#  Section 2: Import the data and load the packages----


#' Import the data.(s3 and s4)
s3 <- read.table(file = "d124_count_stage2.txt", 
                 header = TRUE,
                 na.strings = "NA",
                 stringsAsFactors = TRUE,
                 dec = ".")

#' Check the imported data.
dim(s3)
names(s3)
str(s3)


#' We need the following packages:
#install.packages(Pkgs2Load)    #' Run the first time

Pkgs2Load <- c("lattice", "ggplot2", "mgcv", "plyr", "MASS", 
               "pscl", "glmmTMB", "DHARMa", "cowplot", "rgl", 
               "GGally", "stringi", "performance", "emmeans",
               "scales", "multcomp")
invisible(lapply(Pkgs2Load, library, character.only = TRUE))

#' Load our support file.
source("HighstatLibV13.R")




# Section 3: Prepare the data----

#' Define success as:  germination / number of seeds
s3$success <- s3$count / s3$num 

#' remove samples that are contaminated.
s4<-subset(s3, num!=0)
dim(s4)



# Section 4: Data exploration----


#* Subsection 4.1: Missing values----

colSums(is.na(s4))
#' No missing value


#* Subsection 4.2: Zeros----
#' How many observations are equal to 0?
100 * sum(s4$per == 0) / nrow(s4)
#' 13% of zeros
plot(table(s3$success))
#' plenty of zeros


#* Subsection 4.3: Outliers----

#' Make a Cleveland dotplot for the response variable (per) and the
#' covariate fun, temp, and pop. There are three categorical covariates 
#' in the data set.
#' We rename them.

s4$ffun <- factor(s4$fun) 
table(s4$ffun)

s4$ftemp <- factor(s4$temp) 
table(s4$ftemp)

s4$fpop <- factor(s4$pop) 
table(s4$fpop)

table(s4$ffun, s4$ftemp)
table(s4$ftemp, s4$fpop)
table(s4$ffun, s4$fpop)
#' balanced for an interaction term


#* Subsection 4.4: Collinearity----
#' Make a pairplot of the covariates.
ToPlot <- c("ffun", "ftemp", "fpop")
ggpairs(s4[,ToPlot])
#' balanced. non collinearity

#* Subsection 4.5: Relationships----
#' The figure below shows a graph of the response variable
#' The success of seed germination versus each of the covariates (made with MyMultipanel.ggp2
#'  which is in our support file).
#'  Patterns in fungal species in different temp
p1 <- ggplot(s4, aes(x = ffun, y = per)) + 
  geom_boxplot() + xlab("Fungi")

p2 <- ggplot(s4, aes(x = ffun, y = per)) +
  geom_boxplot() + facet_wrap(~ftemp, ncol = 2)  + xlab("Fungi")

p3 <- ggplot(s4, aes(x = ffun, y = per)) +
  geom_boxplot() + facet_wrap(~fpop, ncol = 2)  + xlab("Fungi")

plot_grid(p1, p2, p3, ncol = 1, labels = c("A", "B", "C"))
#' variable in GC and FF

#' We also show a scatterplot of Depth versus the number of dolphin 
#' sightings for each year. There seems to be a weak non-linear relationship. 

ggplot(s4, aes(x = ffun, y = per)) +
  geom_boxplot() +
  facet_wrap(~ftemp, ncol = 2, scales = "free_y")
#' in GC, higher seed germination using cerato > tul
#' in RT, higher when there is off_cerato

ggplot(s4, aes(x = ffun, y = per)) +
  geom_boxplot() +
  facet_wrap(~fpop, ncol = 2, scales = "free_y")
#' in FF, higher seed germination is obtained with 2 cerato and pev_cerato
#' in MX, higher seed germination is obtained in pff_cerato


#* Subsection 4.6: Dependency----

#' We do not have any dependency until we add in random effect
#' in the experiment. 

# Section 5: Model formulation----

#' We will execute the following binomial GLMM

#' count_i ~ Binomial(Pi_i, num_i)
#' P_i us the number of success
#' E[success_i]   = num_i * Pi_i  
#' var[success_i] = num_i * Pi_i  * (1 - Pi_i)

#' count_i is the number of germinated seed on in plate i.


#'              exp(eta_i)
#' P_i = --------------------
#'           1 + exp(eta_i)

#' where
#'   eta_i = Fun_i * pop_i * temp_i 

#' no random intercept.


# Section 6: Apply the binomial GLMM----

#' Calculate Failure (needed in glmmTMB):
s4$Failure <- s4$num - s4$count

#' Execute the binomial GLMM using glmmTMB
M3 <- glmmTMB(cbind(count, Failure) ~ ftemp*fpop*ffun,
              data = s4,
              family = binomial)

#' Numerical outpout
summary(M3)
#' No odd numbers.

# Section 7: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M1, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M3)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' Overdispersion!

#' We get the same message from DHARMa.
testDispersion(M3)


# Section 8 Model Validation binomial GLMM----


#* Subsection 8.1: Plot residuals vs fitted values----

#' Get the Pearson residuals and the fitted values.
E1 <- resid(M3, type = "pearson")
F1 <- fitted(M3)


#' Plot residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
#' Ok pattern.

#* Subsection 8.1: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' Ok


#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' not sure

#' Residuals vs pop.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' OK?

#' What does DHARMa tell us?
#' Get scaled quantile residuals.
E1Binqr <- simulateResiduals(fittedModel = M1, plot = FALSE)
#' warning message that we have a deviation from uniformity


#' Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(E1Binqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' That is a 'no'.


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E1Binqr, quantreg = TRUE, smoothScatter = FALSE) 
#' Problems for the larger residuals? And also for the smaller residuals
#' Those are the two clouds of residuals in our graph?
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)

plotResiduals(E1Binqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E1Binqr, form = s4$fpop, xlab = "Rank-transformed temp") 
plotResiduals(E1Binqr, form = s4$ftemp, xlab = "Rank-transformed pop") 
#' Trouble for the first and last? There is certainly a non-linear pattern
#' in the residuals with fungi, population, and temp.
#' some within-group residual distributions that are not uniform

#* Subsection 8.2: Simulation ---

#' Simulate 1000 data sets from the binomial GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M1, NSim, size = s4$num, seed = 12345)


#' This is the first simulated data set:
YBin[,1]
#' The first column is the number of successes, and the second column is failure.
#' We will extract the first column, and store it in a matrix Ysim.

Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s4$num)
}


#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0, 1),
     main = "Simulation results")
points(x = sum(s4$count == s4$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)

#' Let's plot the number of zeros in a histogram.
#' table(Zeros)
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(0, 950))
points(x = sum(s4$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)

#' What does DHARMa tell us?
testZeroInflation(M1)
#' we cannot cope with the zeros


#' Expressed as some sort of p-value:
ZerosInData <- sum(s4$count == 0)
sum(Zeros > ZerosInData) / 1000



# Section 9: Visualise the results of the binomial GLMM----


#* Subsection 9.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus distance to the
#' floor, for average values of all other covariates.



#' We will not use the predict function, because we need to do
#' the predictions manually once we reach the zero-inflation
#' model.

#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(fpop = levels(s2$fpop), 
                      ffun = levels(s2$ffun), 
                      ftemp = levels(s2$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s2, 
                .(fpop, ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun * fpop * ftemp, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M1)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M1)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' Back-standardise fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + facet_grid(ftemp~fpop, scales = "fixed")
p1
#' The red dots (including 95% interval) are the predicted
#' values. They are not fit with the observed value. So,
#' this model does not fit


# Section 10: ZIB GLMM---- 


#* Subsection 10.1: Model formulation----
#' A zero-inflated binomial (ZIB)  model is given by the following expression.

#'  Success_i ~ ZIB(P_i, N_i, Pi)
#'  E(success_i)   = (1 - Pi) * mu_i
#'  var(sucess_i) = (1 - Pi) * (V_i + mu_i^2) - (1 - Pi)^2 * mu_i^2

#' Pi_i is the probability of success.
#' N_i is number of seeds in a plate 
#' Pi is the probability of a false zero.

#' And:
#'   mu_i = N_i * P_i  
#'   V_i  = N_i * P_i * (1 - P_i)



#' We will use covariates to model the P_ij
#'  logit(P_ij) = Covariate stuff


#' For the probability of a false zero, we use:
#'  logit(Pi) = Intercept

#' It is possible in glmmTMB to model Pi as a function of covariates.
#' Note: If Pi = 0, then we obtain the ordinary binomial GLMM.


#' One more time: 
#'  -Pi is for the zero-inflation stuff.
#'  -Pi is constant (at least that is what we will do)
#'  -P is for the binomial part. 
#'  -P is a function of covariates and the random effects.
#'  -The logistic link function is used for both terms.



#* Subsection 10.2: Execute the ZIB GLMM----

M4 <- glmmTMB(cbind(count,Failure) ~ fpop*ftemp*ffun,
              family = "binomial",
              ziformula =~ 1,
              data = s4)

#' The numerical output of the NB GLM is as follows.
summary(M4)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M4)


drop1(M4, test = "Chi")
#' The three interactions (ffun * ftemp * fpop) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

#* Subsection 10.3 Assess overdispersion----  

#' Get the fitted values
F1 <- fitted(M4)

#' Get the Pearson residuals
E1 <- resid(M4, type = "pearson")

#' Let´s calculate pearson residuals ourselves
#' Get an X matrix
X <- model.matrix(~ fpop*ftemp*ffun,
                  data = s4)
#' Or: X <- model.matrix(M4)

#' Get the betas for the binomial part.
beta.count <- fixef(M4)$cond

#' Calculate eta = X * beta
eta.count <- X %*% beta.count

#' Calculate P for the binomial part.
P <- exp(eta.count) / (1 + exp(eta.count))

#' Calculate Pi (probability of a false zero).
#' Get the estimated regression part for the binary part.
gamma  <- summary(M4)$coefficients$zi[1,"Estimate"] 

#' And this is Pi.
Pi <- exp(gamma) / (1 + exp(gamma))
Pi
#' Pi is 0.10

#' Recall that these are the mean and variance of a ZIB GLM(M):
#'  E(Y)   <- (1 - Pi) * mu
#'  var(Y) <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2 * mu^2

#'Where:
#'    mu = N * P
#'    V  = N * P * (1 - P) 
#'    N is the  number of seeds
#'    p is the probility of success of the binomial process.


#' We first calculate E(Y).
Ntrials <- s4$num
mu      <- Ntrials * P
ExpY    <- (1 - Pi) * mu

#' Now calculate V and var(Y).
V    <- Ntrials * P * (1 - P)
VarY <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2  * mu^2 


#' And here are the Pearson residuals.
PRes2 <- (s4$count - ExpY) / sqrt(VarY)


#' Assess overdispersion.
Npar <- length(fixef(M4)$cond) + length(fixef(M4)$zi) +1

sum(PRes2^2) / (nrow(s4) - Npar)
# Minor overdispersion


#' What does DHARMa tell us?
testDispersion(M4)
#' Minor overdispersion. 

#' alternative- broken formular
#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M4)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' Overdispersion!


# Section 11: Model Validation ZIB GLMM----     



#* Subsection 11.1: Plot residuals vs fitted values and covariates----

# Residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = ExpY,
     y = PRes2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)            
#' better than the binomial GLMM!

#' What does DHARMa tell us?
E2BINqr <- simulateResiduals(fittedModel = M4, plot = FALSE)
plotQQunif(E2BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' better than the binomial GLMM!


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E2BINqr, quantreg = TRUE, smoothScatter = FALSE) 



# Residuals vs covariates in the model.
MyVar <- c("ffun", "fpop",
           "ftemp")
s4$PRes2 <- PRes2
MyMultipanel.ggp2(Z = s4, 
                  varx = MyVar, 
                  vary = "PRes2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)


plotResiduals(E2BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E2BINqr, form = s4$fpop, xlab = "Rank-transformed temp") 
plotResiduals(E2BINqr, form = s4$ftemp, xlab = "Rank-transformed pop") 
#still acceptable?

#* Subsection 11.2: Simulation for zero-inflation and one-inflation----

testZeroInflation(M4)
#' Pretty accurate and we can skip the simulate 1000 data onwards

#' Alternative if we want to simulate 1000 data sets
#' Simulate 1000 data sets from the ZIB GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M4, NSim, size = s4$num, seed = 12345)


#' We will extract the first column, and store it in a matrix Ysim.
Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Now we have 1,000 simulated data sets from the model. What shall we do 
#' with these simulated data sets? We can calculate the number of zeros in 
#' each of the 1,000 data sets. And we can also count how often we predict
#' a number of success equal to the total (which is a proportion of 1).

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s4$num)
}




#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(-0.5, 0.25),
     main = "Simulation results")
points(x = sum(s4$count == s4$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#' The red dot is the % of ones in the original data set.
#' The model predicts data sets that have just similar number of ones  



#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(1, 50))
points(x = sum(s4$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)
#' The red dot is the number of ones in the original data set.
#' The model predicts data sets that have similar number of zeros
#' as the observed data set.

#' Summarising: The ZIB GLMM can cope with the excessive number of zeros,
#' and also the excessive number of ones.

# Section 12: Visualise the results of the binomial GLMM----


#* Subsection 12.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun to the
#' floor.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(fpop = levels(s4$fpop), 
                      ffun = levels(s4$ffun), 
                      ftemp = levels(s4$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(fpop, ffun, ftemp), summarize,
                stage = seq(from = min(s4$stage), 
                            to   = max(s4$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun * fpop * ftemp, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M4)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M4)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' Back-standardise fun:

p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + facet_grid(ftemp~fpop, scales = "fixed")
p1
#' It is the fitted value for the interaction of ffun, fpop, and ftemp + 95% 
#' interval. It is a bit off in RT

# Section 13: Beta glmm interpretation----

#* Subsection 13.1: Prepare the data----

#' Define the response variable.
#' We will analyse the following response variable: 
#'   Proportion of time that a caribou feeds
#'   This is a number between 0 and 100.
#'   We need to convert it into a number between 0 and 1

s4$PER <- s4$per / 100

#' The beta GLMM does not want to see a 0 or a 1. See the
#' short video on 'What is a beta distribution' on how to
#' deal with this is a quick-and-dirty way. This approach is
#' not recommended if there are plenty of zeros and/or ones.

N <- nrow(s4)
s4$PERT <- (s4$PER * (N - 1) + 0.5 ) / N

#' The transformation changes the numbers in digits number 3 and
#' converts 0 into something small. Same for 1; it becomes 
#' 0.99-something.
cbind(s4$PER, s4$PERT)
#' As you can see, these numbers are very similar.

# Section 14: Model formulation of the beta GLMM----

#' We will execute the following model:
#'  P_PER_j is the proportion of germination for observation j
#'  J is observed plate

#'  P_PER_j ~ beta(Pi_j, phi)
#'  E[P_PER_j]   = Pi_j
#'  var[P_PER_j] = Pi_j * (1 - Pi_j) / (1 + phi)

#'  logit(Pi_i) = Intercept + Covariates


# Section 15: Execute the beta GLMM----


#' Execute the beta GLMM
M5 <- glmmTMB(PERT ~ fpop*ffun*ftemp,
              family = beta_family(link = "logit"), 
              data = s4)

#' Here is the numerical output
summary(M5)

#' The fungi, temp, pop, 2-way (pop/fungi or temp/fungi), and 3-way 
#' interactions are significant at the 5% level.
#' The effect of fungi and temp alone is positive and significant at 
#' the 5% level. Also, the two-way (pop/fungi or temp/fungi) have significantly 
#' negative effects
#' Fungal treatment and seed germination differed at different
#' temperstures

#' Model selection
# Use classical backwards model selection using the AIC:
step(M5)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M5, test = "Chi")
#' The three interactions (ffun * ftemp * fpop) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

add1(M5, scope = ~fpop*ffun*ftemp, test = "Chisq")

glmmTMB:::Anova.glmmTMB(M5, test.statistic = c("Chisq"))
allMeans <- emmeans(M5, specs = c("ffun", "ftemp", "fpop"))
allMeans
pairs(allMeans, comparisons = TRUE)

allMeans <- emmeans(M5, specs = c("ffun", "ftemp", "fpop"))
allMeans
allmeans_pairs <- pairs(allMeans, comparisons = TRUE)
allmeans_pairs_tb <- as.data.frame(allmeans_pairs)
View(allmeans_pairs_tb)
allmeans_pairs_tb <- write.csv(allmeans_pairs_tb, file = "/Volumes/pcooperi_seed/allmeanpairs_stage2_M5.csv")


tempMeans <- emmeans(M5, "ftemp", data = s4)
tempMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M5, "ffun", data = s4)
funMeans
pairs(funMeans, comparisons = TRUE)
popMeans <- emmeans(M5, "fpop", data = s4)
popMeans
pairs(popMeans, comparisons = TRUE)


# Section 16: Model validation----


#' Get the residuals and the fitted values.
E1 <- resid(M5, type = "pearson")
F1 <- fitted(M5)


#* Subsection 16.1: Plot residuals versus fitted values----
par(mfrow = c(1, 1))
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2, col = 1)


#* Subsection 16.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok?


#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' not so okay

#' Residuals vs pop.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' OK?

#* Subsection 16.4 Assess overdispersion----  

# Section 16.3: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M5, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M5)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' underdispersion!

#' We get the same message from DHARMa.
testDispersion(M5)
#'0.82. closed to 1

#' What does DHARMa tell us?
E3BINqr <- simulateResiduals(fittedModel = M5, plot = FALSE)
plotQQunif(E3BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' passed the dispersion test

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E3BINqr, quantreg = TRUE, smoothScatter = FALSE) 


#' Standard model validation steps:
plotResiduals(E3BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E3BINqr, form = s4$fpop, xlab = "Rank-transformed temp") 
plotResiduals(E3BINqr, form = s4$ftemp, xlab = "Rank-transformed pop") 
#' temp has a lot of within- variations 

# Section 17: Simulation----
testZeroInflation(M5)
#' not working with M5. Go stimulate


# The simulation function does also work for the beta GLMM:
Ysim <- simulate(M5, 1000)

#' Now we have 1000 simulated data set. Here is the first one:
Ysim[,1]

#' We can make a histogram of this one.
hist(Ysim[,1])

Interval <- seq(0, 1, by = 0.1)
SimData.int <- cut(Ysim[,1], breaks = Interval)
Tab <- table(SimData.int)
Tab

#' Here is a visualisation of this:
plot(Tab, type = "h", xlab = "Simulated values for coverage", ylab = "Frequency")

#' How is this for the observed data?
ObsData.int <- cut(s4$PERT, breaks = Interval)
TabObs      <- table(ObsData.int)
plot(Tab, type = "h", col = 2, xlab = "Values for coverage", ylab = "Frequency")

# What about doing this 1000 times?
Tab1000 <- matrix(NA, nrow = length(Interval)-1, ncol = 1000)
rownames(Tab1000) <- names(TabObs)
for (i in 1:1000){
  SimData.int <- cut(Ysim[,i], breaks = Interval)
  Tab1000[,i] <- table(SimData.int)
}

#' This is for the first simulated data set:  
Tab1000[,1]

#' This is for the second simulated data set:  
Tab1000[,2]

#' Etc.

MyData <- data.frame(SimulatedFreq = as.vector(Tab1000),
                     ID = rep(rownames(Tab1000), 1000))
MyData$ID <- factor(MyData$ID)
head(MyData,20)

#' Plot the residuals as a boxlot
boxplot(MyData$SimulatedFreq ~ MyData$ID)
points(1:10, TabObs, col = 2, cex = 2, pch = 16)


#' Instead of a boxplot, we can also plot this as a Cleveland dotplot.
#' First make a data frame for the frequencies for the observed coverage data.
DataOriginal <- data.frame(TabObs = TabObs,
                           ID = 1:10)
DataOriginal

#' And plot the simulated frequencies and the observed ones.
p <- ggplot()
p <- p + geom_point(data = MyData,
                    aes(x = SimulatedFreq,
                        y = ID),
                    col = "purple")
p <- p + geom_point(data = DataOriginal,
                    aes(x = TabObs,
                        y = ID),
                    col = "red",
                    size = 4)
p <- p + xlab("Frequeny") + xlab("p_germinating success (in bins")
p
#' The red dot is for the observed data. The purple dots are for the 
#' simulated data.
#' seems like the model is good to predict numbers


# Section 18: Visualise the results of the binomial GLMM----


#* Subsection 18.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun to the
#' floor.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(fpop = levels(s4$fpop), 
                      ffun = levels(s4$ffun), 
                      ftemp = levels(s4$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(fpop, ffun, ftemp), summarize,
                stage = seq(from = min(s4$stage), 
                            to   = max(s4$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun * fpop * ftemp, 
                   data = MyData)
Xp


#'Extract parameters and parameter covariance matrix
Betas    <- fixef(M5)$cond
CovBetas <- vcov(M5)$cond


#'Calculate the fitted values on the predictor scale.
MyData$eta <- Xp %*% Betas
MyData$mu  <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' Calculate the SEs on the scale of the predictor function.
MyData$se    <- sqrt(diag(Xp %*% vcov(M5)$cond %*% t(Xp)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  - 1.96 *MyData$se))
View(MyData)
#write.csv(MyData, file = "/Volumes/pcooperi_seed/predict_stage2.csv")

#' E. Plot everything
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun, colour = ftemp),
                      position=position_dodge(width=1), alpha = 0.3)
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm formation")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = mu, colour = ftemp), 
                      position=position_dodge(width=1))
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo, colour = ftemp),
                         position=position_dodge(width=1))
p1 <- p1 + facet_grid(ftemp~fpop, scales = "fixed") + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#' It is the fitted value for the interaction of ffun, fpop, and ftemp + 95% 
#' interval. 
#' Increasing levels of protocorm formation

#' The left bottom (RT*FF*FUN) has the lowest value and largest variations
#' In comparison of both temperatures, we found higher probability 
#' of protocorm formation (success) in GC. 
#' In comparision across fungi, we found higher probability of 
#' protocorm formation (success) when co-cultured with pff_tul in both 
#' pop and temperatures
#' And the lowest prabability of protocorm formation (success) when
#' there was no fungi (the control treatment). 
#' With the interaction along with temperature and pop and fun, 
#' we found higher probability of seed germination in GC * MX with 
#' pff_tul"
#' 

#' F. Plot everything in 1 graph
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun, colour = ftemp, 
                          shape = fpop),
                      position=position_dodge(width=0.8), alpha = 0.3)
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm formation")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = mu, colour = ftemp,
                          shape = fpop), 
                      position=position_dodge(width=0.8))
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo, colour = ftemp, 
                             shape = fpop),
                         position=position_dodge(width=0.8))
p1 <- p1 + scale_y_continuous(labels = label_number(accuracy=0.01)) + scale_color_manual(values=c("#4EBAD4", "#E64A35"))
p1

#' G. Plot everything in 1 graph
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = interaction(ffun, fpop), colour = ftemp),
                      position=position_dodge(width=0.70), alpha = 0.3)
p1 <- p1 + xlab("") + ylab("Protocorm formation (%)")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = interaction(ffun, fpop), 
                          y = mu, colour = ftemp), 
                      position=position_dodge(width=0.70))
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = interaction(ffun, fpop), 
                             ymax = SeUp, 
                             ymin = SeLo, colour = ftemp),
                         position=position_dodge(width=0.70))
p1 <- p1 + scale_y_continuous(labels = label_number(accuracy=0.01), expand = c(0, 0), limits = c(0, 0.85)) + 
  scale_color_manual(values=c("black", "black")) +
  theme(legend.position="none")
p1

#' H. Plot everything in 1 graph_bar
#' 
s4$Success = s4$count/s4$num
p1 <- ggplot(data = s4, aes(x = interaction(fpop, ffun), y = Success, fill = ftemp))
p1 <- p1 + geom_bar(position=position_dodge(width=0.89), stat = "identity")
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm development")
p1 <- p1 + theme(text = element_text(size = 15)) + 
  scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#4EBAD4", "#E64A35")) +
  expand_limits(y=c(0, 0.7))
p1

#' I. Plot everything in 1 boxplot
#' 
fun_mean <- function(x){
return(round(data.frame(y=mean(x),label=mean(x,na.rm=T)),digit=2))}

s4$Success = s4$count/s4$num
p1 <- ggplot(data = s4, aes(x = interaction(ftemp, fpop), y = Success, fill = ffun))
p1 <- p1 + geom_boxplot()
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm development")
p1 <- p1 + theme(text = element_text(size = 15)) + 
  scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#ffffff", "#ffffff","#ffffff","#ffffff","#ffffff")) +
  expand_limits(y=c(0, 0.7)) + theme(legend.position="none") + stat_summary(fun.data = fun_mean, geom="text", vjust=-0.5, position = position_dodge(width=0.7))
p1

#* Subsection 12.2: summary table of seed germination----
data_summary(s4, varname="Success", groupnames=c("fpop","ffun", "ftemp"))


# Section 19: Zi-Beta glmm interpretation----

#* Section 19.1: Model formulation of the beta GLMM----

#' We will execute the following model:
#'  P_PER_j is the proportion of germination for observation j
#'  J is observed plate

#'  P_PER_j ~ beta(Pi_j, phi)
#'  E[P_PER_j]   = Pi_j
#'  var[P_PER_j] = Pi_j * (1 - Pi_j) / (1 + phi)

#'  logit(Pi_i) = Intercept + Covariates




# Section 20: Execute the Zi-Beta GLMM----


#' Execute the beta GLMM
M6 <- glmmTMB(PERT ~ fpop*ftemp*ffun,
              family = beta_family(link = "logit"), 
              ziformula =~ 1,
              data = s4)

#' Here is the numerical output
summary(M6)

# Section 21: Model validation----

#' Get the residuals and the fitted values.
E1 <- resid(M6, type = "pearson")
F1 <- fitted(M6)


#* Subsection 21.1: Plot residuals versus fitted values----
par(mfrow = c(1, 1))
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2, col = 1)


#* Subsection 21.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok?


#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' not so okay

#' Residuals vs pop.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok?

#* Subsection 21.3 Assess overdispersion----  

# Section 21.4: 1 Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M6, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M6)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' underdispersion!

#' We get the same message from DHARMa.
testDispersion(M6)
#same as beta glmm

#' What does DHARMa tell us?
E4BINqr <- simulateResiduals(fittedModel = M6, plot = FALSE)
plotQQunif(E4BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' better than the binomial GLMM and beta glmm
#' 
#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E4BINqr, quantreg = TRUE, smoothScatter = FALSE) 


#' Standard model validation steps:
plotResiduals(E4BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E4BINqr, form = s4$fpop, xlab = "Rank-transformed temp") 
plotResiduals(E4BINqr, form = s4$ftemp, xlab = "Rank-transformed pop") 
#variation in ftemp

# Section 22: Simulation----
testZeroInflation(M6)
#' not working with M5. Go stimulate


# The simulation function does also work for the beta GLMM:
Ysim <- simulate(M6, 1000)

#' Now we have 1000 simulated data set. Here is the first one:
Ysim[,1]

#' We can make a histogram of this one.
hist(Ysim[,1])

Interval <- seq(0, 1, by = 0.1)
SimData.int <- cut(Ysim[,1], breaks = Interval)
Tab <- table(SimData.int)
Tab

#' Here is a visualisation of this:
plot(Tab, type = "h", xlab = "Simulated values for coverage", ylab = "Frequency")

#' How is this for the observed data?
ObsData.int <- cut(s4$PERT, breaks = Interval)
TabObs      <- table(ObsData.int)
plot(Tab, type = "h", col = 2, xlab = "Values for coverage", ylab = "Frequency")

# What about doing this 1000 times?
Tab1000 <- matrix(NA, nrow = length(Interval)-1, ncol = 1000)
rownames(Tab1000) <- names(TabObs)
for (i in 1:1000){
  SimData.int <- cut(Ysim[,i], breaks = Interval)
  Tab1000[,i] <- table(SimData.int)
}

#' This is for the first simulated data set:  
Tab1000[,1]

#' This is for the second simulated data set:  
Tab1000[,2]

#' Etc.

MyData <- data.frame(SimulatedFreq = as.vector(Tab1000),
                     ID = rep(rownames(Tab1000), 1000))
MyData$ID <- factor(MyData$ID)
head(MyData,20)

#' Plot the residuals as a boxlot
boxplot(MyData$SimulatedFreq ~ MyData$ID)
points(1:10, TabObs, col = 2, cex = 2, pch = 16)


#' Instead of a boxplot, we can also plot this as a Cleveland dotplot.
#' First make a data frame for the frequencies for the observed coverage data.
DataOriginal <- data.frame(TabObs = TabObs,
                           ID = 1:10)
DataOriginal

#' And plot the simulated frequencies and the observed ones.
p <- ggplot()
p <- p + geom_point(data = MyData,
                    aes(x = SimulatedFreq,
                        y = ID),
                    col = "purple")
p <- p + geom_point(data = DataOriginal,
                    aes(x = TabObs,
                        y = ID),
                    col = "red",
                    size = 4)
p <- p + xlab("Frequeny") + xlab("p_germinating success (in bins")
p
#' The red dot is for the observed data. The purple dots are for the 
#' simulated data.
#' seems like the model is good to predict numbers


# Section 23: Visualise the results of the binomial GLMM----


#* Subsection 23.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(fpop = levels(s4$fpop), 
                      ffun = levels(s4$ffun), 
                      ftemp = levels(s4$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(fpop, ffun, ftemp), summarize,
                stage = seq(from = min(s4$stage), 
                            to   = max(s4$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun * fpop * ftemp, 
                   data = MyData)
Xp


#'Extract parameters and parameter covariance matrix
Betas    <- fixef(M6)$cond
CovBetas <- vcov(M6)$cond


#'Calculate the fitted values on the predictor scale.
MyData$eta <- Xp %*% Betas
MyData$mu  <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' Calculate the SEs on the scale of the predictor function.
MyData$se    <- sqrt(diag(Xp %*% vcov(M6)$cond %*% t(Xp)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  - 1.96 *MyData$se))
head(MyData)

#' E. Plot everything
#' fun:

p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = PERT, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = mu), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + facet_grid(ftemp~fpop, scales = "fixed")
p1
#' It is the fitted value for the interaction of ffun, fpop, and ftemp + 95% 
#' interval. It is a bit off in RT*MX*fun
#' 
AIC(M1, M4, M5, M6)
#'M5 has the lowest AIC value and the visulaise makes more sense
#'

################ REMOVE POP #######################################
#' We would like to remove pop because pop does not have any significant 
#' effect on protocorm formation (stage 2)

# Section 24: Model formulation----

#' We will execute the following binomial GLMM

#' count_i ~ Binomial(Pi_i, num_i)
#' P_i us the number of success
#' E[success_i]   = num_i * Pi_i  
#' var[success_i] = num_i * Pi_i  * (1 - Pi_i)

#' count_i is the number of germinated seed on in plate i.


#'              exp(eta_i)
#' P_i = --------------------
#'           1 + exp(eta_i)

#' where
#'   eta_i = Fun_i * temp_i 

#' no random intercept.


# Section 14: Apply the binomial GLMM----

#' Execute the binomial GLMM using glmmTMB
M3a <- glmmTMB(cbind(count, Failure) ~ ftemp*ffun,
               data = s4,
               family = binomial)

#' Numerical outpout
summary(M3a)
#' No odd numbers.

step(M3a)
drop1(M3a, test = "Chi")

# Section 25: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M3a, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M3a)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' Overdispersion!

#' We get the same message from DHARMa.
testDispersion(M3a)
#1.4

# Section 26 Model Validation binomial GLMM----


#* Subsection 26.1: Plot residuals vs fitted values----

#' Get the Pearson residuals and the fitted values.
E1 <- resid(M3a, type = "pearson")
F1 <- fitted(M3a)


#' Plot residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
#' Ok pattern.

#* Subsection 26.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' Ok?


#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' not so much

#' What does DHARMa tell us?
#' Get scaled quantile residuals.
E1Binqr <- simulateResiduals(fittedModel = M3a, plot = FALSE)
#' warning message that we have a deviation from uniformity


#' Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(E1Binqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' That is a 'no'.


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E1Binqr, quantreg = TRUE, smoothScatter = FALSE) 
#' Problems for the larger residuals? And also for the smaller residuals
#' Those are the two clouds of residuals in our graph?
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)

plotResiduals(E1Binqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E1Binqr, form = s4$ftemp, xlab = "Rank-transformed temp") 
#' Trouble for the first and last? There is certainly a non-linear pattern
#' in the residuals with fungi, population, and temp.
#' some within-group residual distributions that are not uniform

#* Subsection 8.2: Simulation ---
testZeroInflation(M3a)

#' alternative please use the below codes
#' Simulate 1000 data sets from the binomial GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M1, NSim, size = s2$num, seed = 12345)


#' This is the first simulated data set:
YBin[,1]
#' The first column is the number of successes, and the second column is failure.
#' We will extract the first column, and store it in a matrix Ysim.

Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s2$num)
}


#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0, 1),
     main = "Simulation results")
points(x = sum(s2$count == s2$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)

#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(0, 950))
points(x = sum(s2$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)

#' What does DHARMa tell us?
testZeroInflation(M1)
#' we cannot cope with the zeros


#' Expressed as some sort of p-value:
ZerosInData <- sum(s2$count == 0)
sum(Zeros > ZerosInData) / 1000



# Section 27: Visualise the results of the binomial GLMM----


#* Subsection 27.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus distance to the
#' floor, for average values of all other covariates.



#' We will not use the predict function, because we need to do
#' the predictions manually once we reach the zero-inflation
#' model.

#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun), 
                      ftemp = levels(s4$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s2, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun * ftemp, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M3a)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M3a)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' Back-standardise fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + facet_grid(~ftemp, scales = "fixed")
p1
#this result is not accurate due to overdispersion

# Section 28: ZIB GLMM---- 


#* Subsection 28.1: Model formulation----
#' A zero-inflated binomial (ZIB)  model is given by the following expression.

#'  Success_i ~ ZIB(P_i, N_i, Pi)
#'  E(success_i)   = (1 - Pi) * mu_i
#'  var(sucess_i) = (1 - Pi) * (V_i + mu_i^2) - (1 - Pi)^2 * mu_i^2

#' Pi_i is the probability of success.
#' N_i is number of seeds in a plate 
#' Pi is the probability of a false zero.

#' And:
#'   mu_i = N_i * P_i  
#'   V_i  = N_i * P_i * (1 - P_i)



#' We will use covariates to model the P_ij
#'  logit(P_ij) = Covariate stuff


#' For the probability of a false zero, we use:
#'  logit(Pi) = Intercept

#' It is possible in glmmTMB to model Pi as a function of covariates.
#' Note: If Pi = 0, then we obtain the ordinary binomial GLMM.


#' One more time: 
#'  -Pi is for the zero-inflation stuff.
#'  -Pi is constant (at least that is what we will do)
#'  -P is for the binomial part. 
#'  -P is a function of covariates and the random effects.
#'  -The logistic link function is used for both terms.



#* Subsection 28.2: Execute the ZIB GLMM----

M4a <- glmmTMB(cbind(count,Failure) ~ ftemp*ffun,
               family = "binomial",
               ziformula =~ 1,
               data = s4)

#' The numerical output of the NB GLM is as follows.
summary(M4a)

glmmTMB:::Anova.glmmTMB(M4a, test.statistic = c("Chisq"))
#' both results are consistent

tempMeans <- emmeans(M4a, "ftemp", data = s2)
tempMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M4a, "ffun", data = s2)
funMeans
pairs(funMeans, comparisons = TRUE)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M4a)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M4a, test = "Chi")
#' The three interactions (ffun * ftemp ) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

#* Subsection 28.3 Assess overdispersion----  

#' Get the fitted values
F1 <- fitted(M4a)

#' Get the Pearson residuals
E1 <- resid(M4a, type = "pearson")

#' Let´s calculate pearson residuals ourselves
#' Get an X matrix
X <- model.matrix(~ ftemp*ffun,
                  data = s4)
#' Or: X <- model.matrix(M2)

#' Get the betas for the binomial part.
beta.count <- fixef(M4a)$cond

#' Calculate eta = X * beta
eta.count <- X %*% beta.count

#' Calculate P for the binomial part.
P <- exp(eta.count) / (1 + exp(eta.count))

#' Calculate Pi (probability of a false zero).
#' Get the estimated regression part for the binary part.
gamma  <- summary(M4a)$coefficients$zi[1,"Estimate"] 

#' And this is Pi.
Pi <- exp(gamma) / (1 + exp(gamma))
Pi
#' Pi is 0.10

#' Recall that these are the mean and variance of a ZIB GLM(M):
#'  E(Y)   <- (1 - Pi) * mu
#'  var(Y) <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2 * mu^2

#'Where:
#'    mu = N * P
#'    V  = N * P * (1 - P) 
#'    N is the  number of seeds
#'    p is the probility of success of the binomial process.


#' We first calculate E(Y).
Ntrials <- s4$num
mu      <- Ntrials * P
ExpY    <- (1 - Pi) * mu

#' Now calculate V and var(Y).
V    <- Ntrials * P * (1 - P)
VarY <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2  * mu^2 


#' And here are the Pearson residuals.
PRes2 <- (s4$count - ExpY) / sqrt(VarY)


#' Assess overdispersion.
Npar <- length(fixef(M4a)$cond) + length(fixef(M4a)$zi) +1

sum(PRes2^2) / (nrow(s4) - Npar)
# minor overdispersion.

#' What does DHARMa tell us?
testDispersion(M4a)
#1.14

# Section 29: Model Validation ZIB GLMM----     



#* Subsection 29.1: Plot residuals vs fitted values and covariates----

# Residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = ExpY,
     y = PRes2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)            
#' better than the binomial GLMM!

#' What does DHARMa tell us?
E2BINqr <- simulateResiduals(fittedModel = M4a, plot = FALSE)
plotQQunif(E2BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' better than the binomial GLMM!


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E2BINqr, quantreg = TRUE, smoothScatter = FALSE) 



# Residuals vs covariates in the model.
MyVar <- c("ffun",
           "ftemp")
s4$PRes2 <- PRes2
MyMultipanel.ggp2(Z = s4, 
                  varx = MyVar, 
                  vary = "PRes2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)


plotResiduals(E2BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E2BINqr, form = s4$ftemp, xlab = "Rank-transformed temp") 
#still acceptable as other parts of DHARMS still look good!

#* Subsection 29.2: Simulation for zero-inflation and one-inflation----

testZeroInflation(M3a)
#' Pretty accurate and we can skip the simulate 1000 data onwards

#' Alternative if we want to simulate 1000 data sets
#' Simulate 1000 data sets from the ZIB GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M2, NSim, size = s2$num, seed = 12345)


#' We will extract the first column, and store it in a matrix Ysim.
Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Now we have 1,000 simulated data sets from the model. What shall we do 
#' with these simulated data sets? We can calculate the number of zeros in 
#' each of the 1,000 data sets. And we can also count how often we predict
#' a number of success equal to the total (which is a proportion of 1).

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s2$num)
}




#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0.1, 0.25),
     main = "Simulation results")
points(x = sum(s2$count == s2$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#' The red dot is the % of ones in the original data set.
#' The model predicts data sets that have just similar number of ones  



#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(300, 390))
points(x = sum(s2$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)
#' The red dot is the number of ones in the original data set.
#' The model predicts data sets that have similar number of zeros
#' as the observed data set.

#' Summarising: The ZIB GLMM can cope with the excessive number of zeros,
#' and also the excessive number of ones.

# Section 30: Visualise the results of the binomial GLMM----


#* Subsection 30.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun), 
                      ftemp = levels(s4$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun * ftemp, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M4a)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M4a)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + facet_grid(~ftemp, scales = "fixed") + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#' It is the fitted value for the interaction of ffun and ftemp + 95% 
#' interval. 

# Section 31: Beta glmm interpretation----

#' Execute the beta GLMM
M5a <- glmmTMB(PERT ~ ffun*ftemp,
               family = beta_family(link = "logit"), 
               data = s4)

#' Here is the numerical output
summary(M5a)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M5a)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M5a, test = "Chi")
#' The three interactions (ffun * ftemp * fpop) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

glmmTMB:::Anova.glmmTMB(M5a, test.statistic = c("Chisq"))

emmeans(M5a, list(pairwise~ffun*ftemp))

allMeans <- emmeans(M5a, specs = c("ffun", "ftemp"))
allMeans
allmeans_pairs <- pairs(allMeans, adjust = "tukey")
allmeans_pairs_tb <- as.data.frame(allmeans_pairs)
View(allmeans_pairs_tb)
#lsmeans(M5a, pairwise~ffun*ftemp)
allmeans_pairs_tb <- write.csv(allmeans_pairs_tb, file = "/Volumes/pcooperi_seed/allmeanpairs_stage2.csv")

popMeans <- emmeans(M5a, "fpop", data = s2)
tempMeans
tempMeans <- emmeans(M5a, "ftemp", data = s2)
tempMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M5a, "ffun", data = s2)
funMeans
pairs(funMeans, comparisons = TRUE)


g1 <- glht(M5a, linfct = mcp(period = "Tukey"))
summary(g1)


plot(allMeans, comparisons= TRUE)
pwpp(allMeans)

# Section 32: Model validation----


#' Get the residuals and the fitted values.
E1 <- resid(M5a, type = "pearson")
F1 <- fitted(M5a)


#* Subsection 32.1: Plot residuals versus fitted values----
par(mfrow = c(1, 1))
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2, col = 1)


#* Subsection 32.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok?


#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' not so okay

#* Subsection 32.3 Assess overdispersion----  

# Section 33: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M5a, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M5a)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' underdispersion!

#' We get the same message from DHARMa.
testDispersion(M5a)
#'0.80 underdispersion

#' What does DHARMa tell us?
E3BINqr <- simulateResiduals(fittedModel = M5a, plot = FALSE)
plotQQunif(E3BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' passed the dispersion test

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E3BINqr, quantreg = TRUE, smoothScatter = FALSE) 


#' Standard model validation steps:
plotResiduals(E3BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E3BINqr, form = s4$ftemp, xlab = "Rank-transformed temp") 
#' temp has a lot of within- variations 

# Section 34: Simulation----
testZeroInflation(M5a)
#' not working with M5a. Go stimulate


# The simulation function does also work for the beta GLMM:
Ysim <- simulate(M5a, 1000)

#' Now we have 1000 simulated data set. Here is the first one:
Ysim[,1]

#' We can make a histogram of this one.
hist(Ysim[,1])

Interval <- seq(0, 1, by = 0.1)
SimData.int <- cut(Ysim[,1], breaks = Interval)
Tab <- table(SimData.int)
Tab

#' Here is a visualisation of this:
plot(Tab, type = "h", xlab = "Simulated values for coverage", ylab = "Frequency")

#' How is this for the observed data?
ObsData.int <- cut(s4$PERT, breaks = Interval)
TabObs      <- table(ObsData.int)
plot(Tab, type = "h", col = 2, xlab = "Values for coverage", ylab = "Frequency")

# What about doing this 1000 times?
Tab1000 <- matrix(NA, nrow = length(Interval)-1, ncol = 1000)
rownames(Tab1000) <- names(TabObs)
for (i in 1:1000){
  SimData.int <- cut(Ysim[,i], breaks = Interval)
  Tab1000[,i] <- table(SimData.int)
}

#' This is for the first simulated data set:  
Tab1000[,1]

#' This is for the second simulated data set:  
Tab1000[,2]

#' Etc.

MyData <- data.frame(SimulatedFreq = as.vector(Tab1000),
                     ID = rep(rownames(Tab1000), 1000))
MyData$ID <- factor(MyData$ID)
head(MyData,20)

#' Plot the residuals as a boxlot
boxplot(MyData$SimulatedFreq ~ MyData$ID)
points(1:10, TabObs, col = 2, cex = 2, pch = 16)


#' Instead of a boxplot, we can also plot this as a Cleveland dotplot.
#' First make a data frame for the frequencies for the observed coverage data.
DataOriginal <- data.frame(TabObs = TabObs,
                           ID = 1:10)
DataOriginal

#' And plot the simulated frequencies and the observed ones.
p <- ggplot()
p <- p + geom_point(data = MyData,
                    aes(x = SimulatedFreq,
                        y = ID),
                    col = "purple")
p <- p + geom_point(data = DataOriginal,
                    aes(x = TabObs,
                        y = ID),
                    col = "red",
                    size = 4)
p <- p + xlab("Frequeny") + xlab("p_germinating success (in bins")
p
#' The red dot is for the observed data. The purple dots are for the 
#' simulated data.
#' seems like the model is good to predict numbers


# Section 35: Visualise the results of the binomial GLMM----


#* Subsection 35.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun to the
#' floor.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun), 
                      ftemp = levels(s4$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s4$stage), 
                            to   = max(s4$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun * ftemp, 
                   data = MyData)
Xp


#'Extract parameters and parameter covariance matrix
Betas    <- fixef(M5a)$cond
CovBetas <- vcov(M5a)$cond


#'Calculate the fitted values on the predictor scale.
MyData$eta <- Xp %*% Betas
MyData$mu  <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' Calculate the SEs on the scale of the predictor function.
MyData$se    <- sqrt(diag(Xp %*% vcov(M5a)$cond %*% t(Xp)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  - 1.96 *MyData$se))
head(MyData)

#' E. Plot everything
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun),
                      position=position_dodge(width=1), alpha = 0.3)
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm formation")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = mu, colour = "red"), 
                      position=position_dodge(width=1))
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo, colour = "red"),
                         position=position_dodge(width=1))
p1 <- p1 + facet_grid(~ftemp, scales = "fixed") + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#' It is the fitted value for the interaction of ffun, fpop, and ftemp + 95% 
#' interval. 

#' plot them all in a graph
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun, colour = ftemp),
                      position=position_dodge(width=0.7), alpha = 0.3)
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm formation")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = mu, colour = ftemp), 
                      position=position_dodge(width=0.7))
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo, colour = ftemp),
                         position=position_dodge(width=0.7))
p1 <- p1 + scale_y_continuous(labels = label_number(accuracy=0.01)) + scale_color_manual(values=c("#4EBAD4", "#E64A35"))
p1

#' G. Plot everything in 1 graph_bar
#' 
s4$Success = s4$count/s4$num
p1 <- ggplot(data = s4, aes(x = ffun, y = Success, fill = ftemp))
p1 <- p1 + geom_bar(position=position_dodge(width=0.89), stat = "identity")
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm development")
p1 <- p1 + theme(text = element_text(size = 15)) + 
  scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#4EBAD4", "#E64A35")) +
  expand_limits(y=c(0, 0.7))
p1

#' H. Plot everything in 1 box_plot
#' 
s4$Success = s4$count/s4$num
p1 <- ggplot(data = s4, aes(x = ffun, y = Success, fill = ftemp))
p1 <- p1 + geom_boxplot()
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm development")
p1 <- p1 + theme(text = element_text(size = 15)) + 
  scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#4EBAD4", "#E64A35")) +
  expand_limits(y=c(-0, 0.7))
p1



#pop
s4_summary <- data_summary(s4, varname="Success", groupnames=c("ffun", "ftemp"))
d124_count_stage1_pop_summary$pop=as.factor(d124_count_stage1_pop_summary$pop)
View(d124_count_stage1_pop_summary)



# Section 36: Zi-Beta glmm interpretation----

#* Section 36.1: Model formulation of the beta GLMM----

#' We will execute the following model:
#'  P_PER_j is the proportion of germination for observation j
#'  J is observed plate

#'  P_PER_j ~ beta(Pi_j, phi)
#'  E[P_PER_j]   = Pi_j
#'  var[P_PER_j] = Pi_j * (1 - Pi_j) / (1 + phi)

#'  logit(Pi_i) = Intercept + Covariates




# Section 37: Execute the Zi-Beta GLMM----


#' Execute the beta GLMM
M6a <- glmmTMB(PERT ~ ftemp*ffun,
               family = beta_family(link = "logit"), 
               ziformula =~ 1,
               data = s4)

#' Here is the numerical output
summary(M6a)

# Section 38: Model validation----

#' Get the residuals and the fitted values.
E1 <- resid(M6a, type = "pearson")
F1 <- fitted(M6a)


#* Subsection 38.1: Plot residuals versus fitted values----
par(mfrow = c(1, 1))
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2, col = 1)


#* Subsection 38.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok?


#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' not so okay

#* Subsection 38.3 Assess overdispersion----  

# Section 38.4: 1 Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M6a, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M6a)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' underdispersion!

#' We get the same message from DHARMa.
testDispersion(M6a)
#same as beta glmm

#' What does DHARMa tell us?
E4BINqr <- simulateResiduals(fittedModel = M6a, plot = FALSE)
plotQQunif(E4BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' better than the binomial GLMM 
#' 
#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E4BINqr, quantreg = TRUE, smoothScatter = FALSE) 


#' Standard model validation steps:
plotResiduals(E4BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E4BINqr, form = s4$ftemp, xlab = "Rank-transformed temp") 
#variation in ftemp

# Section 39: Simulation----
testZeroInflation(M6a)
#' not working with M5. Go stimulate


# The simulation function does also work for the beta GLMM:
Ysim <- simulate(M6a, 1000)

#' Now we have 1000 simulated data set. Here is the first one:
Ysim[,1]

#' We can make a histogram of this one.
hist(Ysim[,1])

Interval <- seq(0, 1, by = 0.1)
SimData.int <- cut(Ysim[,1], breaks = Interval)
Tab <- table(SimData.int)
Tab

#' Here is a visualisation of this:
plot(Tab, type = "h", xlab = "Simulated values for coverage", ylab = "Frequency")

#' How is this for the observed data?
ObsData.int <- cut(s4$PERT, breaks = Interval)
TabObs      <- table(ObsData.int)
plot(Tab, type = "h", col = 2, xlab = "Values for coverage", ylab = "Frequency")

# What about doing this 1000 times?
Tab1000 <- matrix(NA, nrow = length(Interval)-1, ncol = 1000)
rownames(Tab1000) <- names(TabObs)
for (i in 1:1000){
  SimData.int <- cut(Ysim[,i], breaks = Interval)
  Tab1000[,i] <- table(SimData.int)
}

#' This is for the first simulated data set:  
Tab1000[,1]

#' This is for the second simulated data set:  
Tab1000[,2]

#' Etc.

MyData <- data.frame(SimulatedFreq = as.vector(Tab1000),
                     ID = rep(rownames(Tab1000), 1000))
MyData$ID <- factor(MyData$ID)
head(MyData,20)

#' Plot the residuals as a boxlot
boxplot(MyData$SimulatedFreq ~ MyData$ID)
points(1:10, TabObs, col = 2, cex = 2, pch = 16)


#' Instead of a boxplot, we can also plot this as a Cleveland dotplot.
#' First make a data frame for the frequencies for the observed coverage data.
DataOriginal <- data.frame(TabObs = TabObs,
                           ID = 1:10)
DataOriginal

#' And plot the simulated frequencies and the observed ones.
p <- ggplot()
p <- p + geom_point(data = MyData,
                    aes(x = SimulatedFreq,
                        y = ID),
                    col = "purple")
p <- p + geom_point(data = DataOriginal,
                    aes(x = TabObs,
                        y = ID),
                    col = "red",
                    size = 4)
p <- p + xlab("Frequeny") + xlab("p_germinating success (in bins")
p
#' The red dot is for the observed data. The purple dots are for the 
#' simulated data.
#' seems like the model is good to predict numbers


# Section 40: Visualise the results of the binomial GLMM----


#* Subsection 40.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun), 
                      ftemp = levels(s4$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(fpop, ffun, ftemp), summarize,
                stage = seq(from = min(s4$stage), 
                            to   = max(s4$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun*ftemp, 
                   data = MyData)
Xp


#'Extract parameters and parameter covariance matrix
Betas    <- fixef(M6a)$cond
CovBetas <- vcov(M6a)$cond


#'Calculate the fitted values on the predictor scale.
MyData$eta <- Xp %*% Betas
MyData$mu  <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' Calculate the SEs on the scale of the predictor function.
MyData$se    <- sqrt(diag(Xp %*% vcov(M6a)$cond %*% t(Xp)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  - 1.96 *MyData$se))
head(MyData)

#' E. Plot everything
#' fun:

p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = PERT, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = mu), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + facet_grid(~ftemp, scales = "fixed")
p1
#' It is the fitted value for the interaction of ffun, fpop, and ftemp + 95% 
#' interval. It is a bit off because of the dispersion problem
#' 
AIC(M1a, M4a, M5a, M6a)
#'M5a has the lowest AIC value and the visulaise makes more sense



################ REMOVE TEMP #######################################
#' We would like to remove temp because temp does not have any significant 
#' effect on protocorm formation (stage 2)

# Section 41: Model formulation----

#' We will execute the following binomial GLMM

#' count_i ~ Binomial(Pi_i, num_i)
#' P_i us the number of success
#' E[success_i]   = num_i * Pi_i  
#' var[success_i] = num_i * Pi_i  * (1 - Pi_i)

#' count_i is the number of germinated seed on in plate i.


#'              exp(eta_i)
#' P_i = --------------------
#'           1 + exp(eta_i)

#' where
#'   eta_i = Fun_i * pop_i 

#' no random intercept.


# Section 42: Apply the binomial GLMM----

#' Execute the binomial GLMM using glmmTMB
M3b <- glmmTMB(cbind(count, Failure) ~ fpop*ffun,
               data = s4,
               family = binomial)

#' Numerical outpout
summary(M3b)
#' No odd numbers.

step(M3b)
drop1(M3b, test = "Chi")

# Section 43: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M3b, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M3b)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' Overdispersion!

#' We get the same message from DHARMa.
testDispersion(M3b)


# Section 44 Model Validation binomial GLMM----


#* Subsection 44.1: Plot residuals vs fitted values----

#' Get the Pearson residuals and the fitted values.
E1 <- resid(M3b, type = "pearson")
F1 <- fitted(M3b)


#' Plot residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
#' Ok pattern.

#* Subsection 44.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' Ok?


#' Residuals vs pop.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' not so much

#' What does DHARMa tell us?
#' Get scaled quantile residuals.
E1Binqr <- simulateResiduals(fittedModel = M3b, plot = FALSE)
#' warning message that we have a deviation from uniformity


#' Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(E1Binqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' That is a 'no'.


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E1Binqr, quantreg = TRUE, smoothScatter = FALSE) 
#' Problems for the larger residuals? And also for the smaller residuals
#' Those are the two clouds of residuals in our graph?
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)

plotResiduals(E1Binqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E1Binqr, form = s4$fpop, xlab = "Rank-transformed pop") 
#' Trouble for the first and last? There is certainly a non-linear pattern
#' in the residuals with fungi, population, and temp.
#' some within-group residual distributions that are not uniform

#* Subsection 8.2: Simulation ---
testZeroInflation(M3b)

#' alternative please use the below codes
#' Simulate 1000 data sets from the binomial GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M1, NSim, size = s2$num, seed = 12345)


#' This is the first simulated data set:
YBin[,1]
#' The first column is the number of successes, and the second column is failure.
#' We will extract the first column, and store it in a matrix Ysim.

Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s2$num)
}


#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0, 1),
     main = "Simulation results")
points(x = sum(s2$count == s2$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)

#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(0, 950))
points(x = sum(s2$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)

#' What does DHARMa tell us?
testZeroInflation(M1)
#' we cannot cope with the zeros


#' Expressed as some sort of p-value:
ZerosInData <- sum(s2$count == 0)
sum(Zeros > ZerosInData) / 1000



# Section 45: Visualise the results of the binomial GLMM----


#* Subsection 45.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus distance to the
#' floor, for average values of all other covariates.



#' We will not use the predict function, because we need to do
#' the predictions manually once we reach the zero-inflation
#' model.

#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun), 
                      fpop = levels(s4$fpop))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s2, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun * fpop, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M3b)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M3b)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' Back-standardise fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + facet_grid(~fpop, scales = "fixed")
p1
#this result is not accurate due to overdispersion

# Section 46: ZIB GLMM---- 


#* Subsection 46.1: Model formulation----
#' A zero-inflated binomial (ZIB)  model is given by the following expression.

#'  Success_i ~ ZIB(P_i, N_i, Pi)
#'  E(success_i)   = (1 - Pi) * mu_i
#'  var(sucess_i) = (1 - Pi) * (V_i + mu_i^2) - (1 - Pi)^2 * mu_i^2

#' Pi_i is the probability of success.
#' N_i is number of seeds in a plate 
#' Pi is the probability of a false zero.

#' And:
#'   mu_i = N_i * P_i  
#'   V_i  = N_i * P_i * (1 - P_i)



#' We will use covariates to model the P_ij
#'  logit(P_ij) = Covariate stuff


#' For the probability of a false zero, we use:
#'  logit(Pi) = Intercept

#' It is possible in glmmTMB to model Pi as a function of covariates.
#' Note: If Pi = 0, then we obtain the ordinary binomial GLMM.


#' One more time: 
#'  -Pi is for the zero-inflation stuff.
#'  -Pi is constant (at least that is what we will do)
#'  -P is for the binomial part. 
#'  -P is a function of covariates and the random effects.
#'  -The logistic link function is used for both terms.



#* Subsection 46.2: Execute the ZIB GLMM----

M4b <- glmmTMB(cbind(count,Failure) ~ fpop*ffun,
               family = "binomial",
               ziformula =~ 1,
               data = s4)

#' The numerical output of the NB GLM is as follows.
summary(M4b)

glmmTMB:::Anova.glmmTMB(M4b, test.statistic = c("Chisq"))
#' both results are consistent

tempMeans <- emmeans(M4b, "ftemp", data = s2)
tempMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M4b, "ffun", data = s2)
funMeans
pairs(funMeans, comparisons = TRUE)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M4b)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M4b, test = "Chi")
#' The three interactions (ffun * ftemp ) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

#* Subsection 46.3 Assess overdispersion----  

#' Get the fitted values
F1 <- fitted(M4b)

#' Get the Pearson residuals
E1 <- resid(M4b, type = "pearson")

#' Let´s calculate pearson residuals ourselves
#' Get an X matrix
X <- model.matrix(~ fpop*ffun,
                  data = s4)
#' Or: X <- model.matrix(M2)

#' Get the betas for the binomial part.
beta.count <- fixef(M4b)$cond

#' Calculate eta = X * beta
eta.count <- X %*% beta.count

#' Calculate P for the binomial part.
P <- exp(eta.count) / (1 + exp(eta.count))

#' Calculate Pi (probability of a false zero).
#' Get the estimated regression part for the binary part.
gamma  <- summary(M4b)$coefficients$zi[1,"Estimate"] 

#' And this is Pi.
Pi <- exp(gamma) / (1 + exp(gamma))
Pi
#' Pi is 0.12

#' Recall that these are the mean and variance of a ZIB GLM(M):
#'  E(Y)   <- (1 - Pi) * mu
#'  var(Y) <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2 * mu^2

#'Where:
#'    mu = N * P
#'    V  = N * P * (1 - P) 
#'    N is the  number of seeds
#'    p is the probility of success of the binomial process.


#' We first calculate E(Y).
Ntrials <- s4$num
mu      <- Ntrials * P
ExpY    <- (1 - Pi) * mu

#' Now calculate V and var(Y).
V    <- Ntrials * P * (1 - P)
VarY <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2  * mu^2 


#' And here are the Pearson residuals.
PRes2 <- (s4$count - ExpY) / sqrt(VarY)


#' Assess overdispersion.
Npar <- length(fixef(M4b)$cond) + length(fixef(M4b)$zi) +1

sum(PRes2^2) / (nrow(s4) - Npar)
# minor overdispersion.

#' What does DHARMa tell us?
testDispersion(M4b)

# Section 47: Model Validation ZIB GLMM----     



#* Subsection 47.1: Plot residuals vs fitted values and covariates----

# Residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = ExpY,
     y = PRes2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)            
#' better than the binomial GLMM!

#' What does DHARMa tell us?
E2BINqr <- simulateResiduals(fittedModel = M4b, plot = FALSE)
plotQQunif(E2BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' better than the binomial GLMM!


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E2BINqr, quantreg = TRUE, smoothScatter = FALSE) 



# Residuals vs covariates in the model.
MyVar <- c("ffun",
           "fpop")
s4$PRes2 <- PRes2
MyMultipanel.ggp2(Z = s4, 
                  varx = MyVar, 
                  vary = "PRes2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)


plotResiduals(E2BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E2BINqr, form = s4$fpop, xlab = "Rank-transformed pop") 
#still acceptable as other parts of DHARMS still look good!

#* Subsection 47.2: Simulation for zero-inflation and one-inflation----

testZeroInflation(M4b)
#' Pretty accurate and we can skip the simulate 1000 data onwards

#' Alternative if we want to simulate 1000 data sets
#' Simulate 1000 data sets from the ZIB GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M2, NSim, size = s2$num, seed = 12345)


#' We will extract the first column, and store it in a matrix Ysim.
Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Now we have 1,000 simulated data sets from the model. What shall we do 
#' with these simulated data sets? We can calculate the number of zeros in 
#' each of the 1,000 data sets. And we can also count how often we predict
#' a number of success equal to the total (which is a proportion of 1).

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s2$num)
}




#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0.1, 0.25),
     main = "Simulation results")
points(x = sum(s2$count == s2$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#' The red dot is the % of ones in the original data set.
#' The model predicts data sets that have just similar number of ones  



#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(300, 390))
points(x = sum(s2$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)
#' The red dot is the number of ones in the original data set.
#' The model predicts data sets that have similar number of zeros
#' as the observed data set.

#' Summarising: The ZIB GLMM can cope with the excessive number of zeros,
#' and also the excessive number of ones.

# Section 48: Visualise the results of the binomial GLMM----


#* Subsection 48.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun), 
                      fpop = levels(s4$fpop))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun * fpop, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M4b)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M4b)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + facet_grid(~fpop, scales = "fixed") + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#' It is the fitted value for the interaction of ffun and ftemp + 95% 
#' interval. 

# Section 49: Beta glmm interpretation----

#' Execute the beta GLMM
M5b <- glmmTMB(PERT ~ ffun*fpop,
               family = beta_family(link = "logit"), 
               data = s4)

#' Here is the numerical output
summary(M5b)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M5b)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M5b, test = "Chi")
#' The three interactions (ffun * fpop) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.


glmmTMB:::Anova.glmmTMB(M5b, test.statistic = c("Chisq"))
popMeans <- emmeans(M5b, "fpop", data = s4)
popMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M5b, "ffun", data = s4)
funMeans
pairs(funMeans, comparisons = TRUE)


# Section 50: Model validation----


#' Get the residuals and the fitted values.
E1 <- resid(M5b, type = "pearson")
F1 <- fitted(M5b)


#* Subsection 32.1: Plot residuals versus fitted values----
par(mfrow = c(1, 1))
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2, col = 1)


#* Subsection 50.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok?


#' Residuals vs pop.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' not so okay

#* Subsection 50.3 Assess overdispersion----  

# Section 51: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M5b, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M5b)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' underdispersion!

#' We get the same message from DHARMa.
testDispersion(M5b)
#'1.03 underdispersion

#' What does DHARMa tell us?
E3BINqr <- simulateResiduals(fittedModel = M5b, plot = FALSE)
plotQQunif(E3BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' passed the dispersion test

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E3BINqr, quantreg = TRUE, smoothScatter = FALSE) 


#' Standard model validation steps:
plotResiduals(E3BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E3BINqr, form = s4$fpop, xlab = "Rank-transformed pop") 
#' temp has a lot of within- variations 

# Section 52: Simulation----
testZeroInflation(M5b)
#' not working with M5b. Go stimulate


# The simulation function does also work for the beta GLMM:
Ysim <- simulate(M5b, 1000)

#' Now we have 1000 simulated data set. Here is the first one:
Ysim[,1]

#' We can make a histogram of this one.
hist(Ysim[,1])

Interval <- seq(0, 1, by = 0.1)
SimData.int <- cut(Ysim[,1], breaks = Interval)
Tab <- table(SimData.int)
Tab

#' Here is a visualisation of this:
plot(Tab, type = "h", xlab = "Simulated values for coverage", ylab = "Frequency")

#' How is this for the observed data?
ObsData.int <- cut(s4$PERT, breaks = Interval)
TabObs      <- table(ObsData.int)
plot(Tab, type = "h", col = 2, xlab = "Values for coverage", ylab = "Frequency")

# What about doing this 1000 times?
Tab1000 <- matrix(NA, nrow = length(Interval)-1, ncol = 1000)
rownames(Tab1000) <- names(TabObs)
for (i in 1:1000){
  SimData.int <- cut(Ysim[,i], breaks = Interval)
  Tab1000[,i] <- table(SimData.int)
}

#' This is for the first simulated data set:  
Tab1000[,1]

#' This is for the second simulated data set:  
Tab1000[,2]

#' Etc.

MyData <- data.frame(SimulatedFreq = as.vector(Tab1000),
                     ID = rep(rownames(Tab1000), 1000))
MyData$ID <- factor(MyData$ID)
head(MyData,20)

#' Plot the residuals as a boxlot
boxplot(MyData$SimulatedFreq ~ MyData$ID)
points(1:10, TabObs, col = 2, cex = 2, pch = 16)


#' Instead of a boxplot, we can also plot this as a Cleveland dotplot.
#' First make a data frame for the frequencies for the observed coverage data.
DataOriginal <- data.frame(TabObs = TabObs,
                           ID = 1:10)
DataOriginal

#' And plot the simulated frequencies and the observed ones.
p <- ggplot()
p <- p + geom_point(data = MyData,
                    aes(x = SimulatedFreq,
                        y = ID),
                    col = "purple")
p <- p + geom_point(data = DataOriginal,
                    aes(x = TabObs,
                        y = ID),
                    col = "red",
                    size = 4)
p <- p + xlab("Frequeny") + xlab("p_protocorm development success (in bins")
p
#' The red dot is for the observed data. The purple dots are for the 
#' simulated data.
#' seems like the model is good to predict numbers


# Section 53: Visualise the results of the binomial GLMM----


#* Subsection 53.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun), 
                      fpop = levels(s4$fpop))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun, fpop), summarize,
                stage = seq(from = min(s4$stage), 
                            to   = max(s4$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun * fpop, 
                   data = MyData)
Xp


#'Extract parameters and parameter covariance matrix
Betas    <- fixef(M5b)$cond
CovBetas <- vcov(M5b)$cond


#'Calculate the fitted values on the predictor scale.
MyData$eta <- Xp %*% Betas
MyData$mu  <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' Calculate the SEs on the scale of the predictor function.
MyData$se    <- sqrt(diag(Xp %*% vcov(M5b)$cond %*% t(Xp)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  - 1.96 *MyData$se))
head(MyData)

#' E. Plot everything
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun),
                      position=position_dodge(width=1), alpha = 0.3)
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm formation")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = mu, colour = "red"), 
                      position=position_dodge(width=1))
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo, colour = "red"),
                         position=position_dodge(width=1))
p1 <- p1 + facet_grid(~fpop, scales = "fixed") + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#' It is the fitted value for the interaction of ffun, fpop, and ftemp + 95% 
#' interval. 
#' 
#' 
#' plot them all in a graph
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun, colour = fpop),
                      position=position_dodge(width=1), alpha = 0.3)
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm formation")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = mu, colour = fpop), 
                      position=position_dodge(width=1))
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo, colour = fpop),
                         position=position_dodge(width=1))
p1 <- p1 + scale_y_continuous(labels = label_number(accuracy=0.01)) + scale_color_manual(values=c("#4bcdb3", "#454648"))
p1
#' It is the fitted value for the interaction of ffun and fpop + 95% 
#' interval. 
#' 
################ REMOVE POP and TEMP #######################################
#' We would like to remove temp and pop because both do not have any significant 
#' effect on protocorm formation (stage 2)

# Section 54: Model formulation----

#' We will execute the following binomial GLMM

#' count_i ~ Binomial(Pi_i, num_i)
#' P_i us the number of success
#' E[success_i]   = num_i * Pi_i  
#' var[success_i] = num_i * Pi_i  * (1 - Pi_i)

#' count_i is the number of germinated seed on in plate i.


#'              exp(eta_i)
#' P_i = --------------------
#'           1 + exp(eta_i)

#' where
#'   eta_i = Fun_i * temp_i 

#' no random intercept.


# Section 55: Apply the binomial GLMM----

#' Execute the binomial GLMM using glmmTMB
M3c <- glmmTMB(cbind(count, Failure) ~ ffun,
               data = s4,
               family = binomial)

#' Numerical outpout
summary(M3c)
#' No odd numbers.

step(M3c)
drop1(M3c, test = "Chi")

# Section 56: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M3c, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M3c)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' Overdispersion!

#' We get the same message from DHARMa.
testDispersion(M3c)


# Section 57 Model Validation binomial GLMM----


#* Subsection 57.1: Plot residuals vs fitted values----

#' Get the Pearson residuals and the fitted values.
E1 <- resid(M3c, type = "pearson")
F1 <- fitted(M3c)


#' Plot residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
#' Ok pattern.

#* Subsection 57.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' Ok?

#' What does DHARMa tell us?
#' Get scaled quantile residuals.
E1Binqr <- simulateResiduals(fittedModel = M3c, plot = FALSE)
#' warning message that we have a deviation from uniformity


#' Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(E1Binqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' That is a 'no'.


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E1Binqr, quantreg = TRUE, smoothScatter = FALSE) 
#' Problems for the larger residuals? And also for the smaller residuals
#' Those are the two clouds of residuals in our graph?
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)

plotResiduals(E1Binqr, form = s4$ffun, xlab = "Rank-transformed fun") 
#' Trouble for the first and last? There is certainly a non-linear pattern
#' in the residuals with fungi, population, and temp.
#' some within-group residual distributions that are not uniform

#* Subsection 8.2: Simulation ---
testZeroInflation(M3c)

#' alternative please use the below codes
#' Simulate 1000 data sets from the binomial GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M1, NSim, size = s4$num, seed = 12345)


#' This is the first simulated data set:
YBin[,1]
#' The first column is the number of successes, and the second column is failure.
#' We will extract the first column, and store it in a matrix Ysim.

Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s4$num)
}


#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0, 1),
     main = "Simulation results")
points(x = sum(s4$count == s4$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)

#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(0, 950))
points(x = sum(s2$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)

#' What does DHARMa tell us?
testZeroInflation(M3c)
#' we cannot cope with the zeros


#' Expressed as some sort of p-value:
ZerosInData <- sum(s4$count == 0)
sum(Zeros > ZerosInData) / 1000



# Section 58: Visualise the results of the binomial GLMM----


#* Subsection 58.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus distance to the
#' floor, for average values of all other covariates.



#' We will not use the predict function, because we need to do
#' the predictions manually once we reach the zero-inflation
#' model.

#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M3c)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M3c)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' Back-standardise fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1
#this result is not accurate due to overdispersion

# Section 59: ZIB GLMM---- 


#* Subsection 59.1: Model formulation----
#' A zero-inflated binomial (ZIB)  model is given by the following expression.

#'  Success_i ~ ZIB(P_i, N_i, Pi)
#'  E(success_i)   = (1 - Pi) * mu_i
#'  var(sucess_i) = (1 - Pi) * (V_i + mu_i^2) - (1 - Pi)^2 * mu_i^2

#' Pi_i is the probability of success.
#' N_i is number of seeds in a plate 
#' Pi is the probability of a false zero.

#' And:
#'   mu_i = N_i * P_i  
#'   V_i  = N_i * P_i * (1 - P_i)



#' We will use covariates to model the P_ij
#'  logit(P_ij) = Covariate stuff


#' For the probability of a false zero, we use:
#'  logit(Pi) = Intercept

#' It is possible in glmmTMB to model Pi as a function of covariates.
#' Note: If Pi = 0, then we obtain the ordinary binomial GLMM.


#' One more time: 
#'  -Pi is for the zero-inflation stuff.
#'  -Pi is constant (at least that is what we will do)
#'  -P is for the binomial part. 
#'  -P is a function of covariates and the random effects.
#'  -The logistic link function is used for both terms.



#* Subsection 59.2: Execute the ZIB GLMM----

M4c <- glmmTMB(cbind(count,Failure) ~ ffun,
               family = "binomial",
               ziformula =~ 1,
               data = s4)

#' The numerical output of the NB GLM is as follows.
summary(M4c)

glmmTMB:::Anova.glmmTMB(M4c, test.statistic = c("Chisq"))
#' both results are consistent

tempMeans <- emmeans(M4c, "ftemp", data = s2)
tempMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M4c, "ffun", data = s2)
funMeans
pairs(funMeans, comparisons = TRUE)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M4c)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M4c, test = "Chi")
#' The three interactions (ffun * ftemp ) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

#* Subsection 59.3 Assess overdispersion----  

#' Get the fitted values
F1 <- fitted(M4c)

#' Get the Pearson residuals
E1 <- resid(M4c, type = "pearson")

#' Let´s calculate pearson residuals ourselves
#' Get an X matrix
X <- model.matrix(~ ffun,
                  data = s4)
#' Or: X <- model.matrix(M2)

#' Get the betas for the binomial part.
beta.count <- fixef(M4c)$cond

#' Calculate eta = X * beta
eta.count <- X %*% beta.count

#' Calculate P for the binomial part.
P <- exp(eta.count) / (1 + exp(eta.count))

#' Calculate Pi (probability of a false zero).
#' Get the estimated regression part for the binary part.
gamma  <- summary(M4c)$coefficients$zi[1,"Estimate"] 

#' And this is Pi.
Pi <- exp(gamma) / (1 + exp(gamma))
Pi
#' Pi is 0.12

#' Recall that these are the mean and variance of a ZIB GLM(M):
#'  E(Y)   <- (1 - Pi) * mu
#'  var(Y) <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2 * mu^2

#'Where:
#'    mu = N * P
#'    V  = N * P * (1 - P) 
#'    N is the  number of seeds
#'    p is the probility of success of the binomial process.


#' We first calculate E(Y).
Ntrials <- s4$num
mu      <- Ntrials * P
ExpY    <- (1 - Pi) * mu

#' Now calculate V and var(Y).
V    <- Ntrials * P * (1 - P)
VarY <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2  * mu^2 


#' And here are the Pearson residuals.
PRes2 <- (s4$count - ExpY) / sqrt(VarY)


#' Assess overdispersion.
Npar <- length(fixef(M4c)$cond) + length(fixef(M4c)$zi) +1

sum(PRes2^2) / (nrow(s4) - Npar)
# minor overdispersion.

#' What does DHARMa tell us?
testDispersion(M4c)

# Section 60: Model Validation ZIB GLMM----     



#* Subsection 60.1: Plot residuals vs fitted values and covariates----

# Residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = ExpY,
     y = PRes2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)            
#' better than the binomial GLMM!

#' What does DHARMa tell us?
E2BINqr <- simulateResiduals(fittedModel = M4c, plot = FALSE)
plotQQunif(E2BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' better than the binomial GLMM!


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E2BINqr, quantreg = TRUE, smoothScatter = FALSE) 



# Residuals vs covariates in the model.
MyVar <- c("ffun")
s4$PRes2 <- PRes2
MyMultipanel.ggp2(Z = s4, 
                  varx = MyVar, 
                  vary = "PRes2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)


plotResiduals(E2BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
#still acceptable as other parts of DHARMS still look good!

#* Subsection 60.2: Simulation for zero-inflation and one-inflation----

testZeroInflation(M3c)
#' Pretty accurate and we can skip the simulate 1000 data onwards

#' Alternative if we want to simulate 1000 data sets
#' Simulate 1000 data sets from the ZIB GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M2, NSim, size = s2$num, seed = 12345)


#' We will extract the first column, and store it in a matrix Ysim.
Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Now we have 1,000 simulated data sets from the model. What shall we do 
#' with these simulated data sets? We can calculate the number of zeros in 
#' each of the 1,000 data sets. And we can also count how often we predict
#' a number of success equal to the total (which is a proportion of 1).

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s2$num)
}




#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0.1, 0.25),
     main = "Simulation results")
points(x = sum(s2$count == s2$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#' The red dot is the % of ones in the original data set.
#' The model predicts data sets that have just similar number of ones  



#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(300, 390))
points(x = sum(s2$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)
#' The red dot is the number of ones in the original data set.
#' The model predicts data sets that have similar number of zeros
#' as the observed data set.

#' Summarising: The ZIB GLMM can cope with the excessive number of zeros,
#' and also the excessive number of ones.

# Section 61: Visualise the results of the binomial GLMM----


#* Subsection 61.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M4c)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M4c)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#' It is the fitted value for the interaction of ffun and ftemp + 95% 
#' interval. 

# Section 62: Beta glmm interpretation----

#' Execute the beta GLMM
M5c <- glmmTMB(PERT ~ ffun,
               family = beta_family(link = "logit"), 
               data = s4)

#' Here is the numerical output
summary(M5c)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M5c)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M5c, test = "Chi")
#' The three interactions (ffun * ftemp * fpop) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

glmmTMB:::Anova.glmmTMB(M5c, test.statistic = c("Chisq"))
tempMeans <- emmeans(M5c, "ftemp", data = s2)
tempMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M5c, "ffun", data = s2)
funMeans
pairs(funMeans, comparisons = TRUE)


# Section 63: Model validation----


#' Get the residuals and the fitted values.
E1 <- resid(M5c, type = "pearson")
F1 <- fitted(M5c)


#* Subsection 63.1: Plot residuals versus fitted values----
par(mfrow = c(1, 1))
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2, col = 1)


#* Subsection 63.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok?

#* Subsection 63.3 Assess overdispersion----  

# Section 64: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M5c, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M5c)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' 0.21

#' We get the same message from DHARMa.
testDispersion(M5c)
#'0.78 underdispersion

#' What does DHARMa tell us?
E3BINqr <- simulateResiduals(fittedModel = M5c, plot = FALSE)
plotQQunif(E3BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' passed the dispersion test

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E3BINqr, quantreg = TRUE, smoothScatter = FALSE) 


#' Standard model validation steps:
plotResiduals(E3BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
#' temp has a lot of within- variations 

# Section 65: Simulation----
testZeroInflation(M5c)
#' not working with M5c. Go stimulate


# The simulation function does also work for the beta GLMM:
Ysim <- simulate(M5c, 1000)

#' Now we have 1000 simulated data set. Here is the first one:
Ysim[,1]

#' We can make a histogram of this one.
hist(Ysim[,1])

Interval <- seq(0, 1, by = 0.1)
SimData.int <- cut(Ysim[,1], breaks = Interval)
Tab <- table(SimData.int)
Tab

#' Here is a visualisation of this:
plot(Tab, type = "h", xlab = "Simulated values for coverage", ylab = "Frequency")

#' How is this for the observed data?
ObsData.int <- cut(s4$PERT, breaks = Interval)
TabObs      <- table(ObsData.int)
plot(Tab, type = "h", col = 2, xlab = "Values for coverage", ylab = "Frequency")

# What about doing this 1000 times?
Tab1000 <- matrix(NA, nrow = length(Interval)-1, ncol = 1000)
rownames(Tab1000) <- names(TabObs)
for (i in 1:1000){
  SimData.int <- cut(Ysim[,i], breaks = Interval)
  Tab1000[,i] <- table(SimData.int)
}

#' This is for the first simulated data set:  
Tab1000[,1]

#' This is for the second simulated data set:  
Tab1000[,2]

#' Etc.

MyData <- data.frame(SimulatedFreq = as.vector(Tab1000),
                     ID = rep(rownames(Tab1000), 1000))
MyData$ID <- factor(MyData$ID)
head(MyData,20)

#' Plot the residuals as a boxlot
boxplot(MyData$SimulatedFreq ~ MyData$ID)
points(1:10, TabObs, col = 2, cex = 2, pch = 16)


#' Instead of a boxplot, we can also plot this as a Cleveland dotplot.
#' First make a data frame for the frequencies for the observed coverage data.
DataOriginal <- data.frame(TabObs = TabObs,
                           ID = 1:10)
DataOriginal

#' And plot the simulated frequencies and the observed ones.
p <- ggplot()
p <- p + geom_point(data = MyData,
                    aes(x = SimulatedFreq,
                        y = ID),
                    col = "purple")
p <- p + geom_point(data = DataOriginal,
                    aes(x = TabObs,
                        y = ID),
                    col = "red",
                    size = 4)
p <- p + xlab("Frequeny") + xlab("p_portocorm development success (in bins")
p
#' The red dot is for the observed data. The purple dots are for the 
#' simulated data.
#' seems like the model is good to predict numbers


# Section 66: Visualise the results of the binomial GLMM----


#* Subsection 66.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun to the
#' floor.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun), summarize,
                stage = seq(from = min(s4$stage), 
                            to   = max(s4$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ffun, 
                   data = MyData)
Xp


#'Extract parameters and parameter covariance matrix
Betas    <- fixef(M5c)$cond
CovBetas <- vcov(M5c)$cond


#'Calculate the fitted values on the predictor scale.
MyData$eta <- Xp %*% Betas
MyData$mu  <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' Calculate the SEs on the scale of the predictor function.
MyData$se    <- sqrt(diag(Xp %*% vcov(M5c)$cond %*% t(Xp)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  - 1.96 *MyData$se))
head(MyData)

#' E. Plot everything
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun),
                      position=position_dodge(width=1), alpha = 0.3)
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm formation")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = mu, colour = "red"), 
                      position=position_dodge(width=1))
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo, colour = "red"),
                         position=position_dodge(width=1))
p1 <- p1  + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#' It is the fitted value for the interaction of ffun, fpop, and ftemp + 95% 
#' interval. 

AIC(M5, M5b, M5c)
#'M5c has the lowest AIC value and the visulaise makes more sense

################ drop() #######################################
#' drop

# Section 67: Model formulation----

#' We will execute the following binomial GLMM

#' count_i ~ Binomial(Pi_i, num_i)
#' P_i us the number of success
#' E[success_i]   = num_i * Pi_i  
#' var[success_i] = num_i * Pi_i  * (1 - Pi_i)

#' count_i is the number of germinated seed on in plate i.


#'              exp(eta_i)
#' P_i = --------------------
#'           1 + exp(eta_i)

#' where
#'   eta_i = Fun_i * temp_i 

#' no random intercept.


# Section 68: Apply the binomial GLMM----

#' Execute the binomial GLMM using glmmTMB
M3d <- glmmTMB(cbind(count, Failure) ~ ffun + ftemp + fpop + ftemp*fpop + ftemp*ffun + fpop*ffun,
               data = s4,
               family = binomial)

#' Numerical outpout
summary(M3d)
#' No odd numbers.

step(M3d)
drop1(M3d, test = "Chi")
drop1(M3d, .~., test = "Chi")

#' Execute the binomial GLMM using glmmTMB
M3d1 <- glmmTMB(cbind(count, Failure) ~ ftemp*fpop + ftemp*ffun + fpop*ffun,
                data = s4,
                family = binomial)

#' Numerical outpout
summary(M3d1)
#' No odd numbers.

step(M3d1)
drop1(M3d1, test = "Chi")

# Section 69: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M3d1, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M3d1)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' Overdispersion!

#' We get the same message from DHARMa.
testDispersion(M3d1)
#1.4

# Section 70: Model Validation binomial GLMM----


#* Subsection 70.1: Plot residuals vs fitted values----

#' Get the Pearson residuals and the fitted values.
E1 <- resid(M3d1, type = "pearson")
F1 <- fitted(M3d1)


#' Plot residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
#' Ok pattern.

#* Subsection 70.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' Ok?


#' Residuals vs temp.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' not good

#' Residuals vs pop.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok

#' What does DHARMa tell us?
#' Get scaled quantile residuals.
E1Binqr <- simulateResiduals(fittedModel = M3d1, plot = FALSE)
#' warning message that we have a deviation from uniformity


#' Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(E1Binqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' That is a 'no'.


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E1Binqr, quantreg = TRUE, smoothScatter = FALSE) 
#' Problems for the larger residuals? And also for the smaller residuals
#' Those are the two clouds of residuals in our graph?
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)

plotResiduals(E1Binqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E1Binqr, form = s4$ftemp, xlab = "Rank-transformed fun") 
plotResiduals(E1Binqr, form = s4$fpop, xlab = "Rank-transformed fun") 
#' Trouble for the first and last? There is certainly a non-linear pattern
#' in the residuals with fungi, population, and temp.
#' some within-group residual distributions that are not uniform

#* Subsection 70.2: Simulation ---
testZeroInflation(M3d1)
#' we cannot cope with the zeros

#' alternative please use the below codes
#' Simulate 1000 data sets from the binomial GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M3d, NSim, size = s4$num, seed = 12345)


#' This is the first simulated data set:
YBin[,1]
#' The first column is the number of successes, and the second column is failure.
#' We will extract the first column, and store it in a matrix Ysim.

Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s4$num)
}


#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0, 1),
     main = "Simulation results")
points(x = sum(s4$count == s4$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)

#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(0, 950))
points(x = sum(s2$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)

#' Expressed as some sort of p-value:
ZerosInData <- sum(s4$count == 0)
sum(Zeros > ZerosInData) / 1000


# Section 71: Visualise the results of the binomial GLMM----


#* Subsection 71.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus distance to the
#' floor, for average values of all other covariates.



#' We will not use the predict function, because we need to do
#' the predictions manually once we reach the zero-inflation
#' model.

#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun),
                      ftemp = levels(s4$ftemp),
                      fpop = levels(s4$fpop))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ftemp*fpop + ftemp*ffun + fpop*ffun, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M3d1)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M3d1)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' Back-standardise fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + facet_grid(ftemp~fpop, scales = "fixed") + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#this result is not accurate due to overdispersion

# Section 72: ZIB GLMM---- 


#* Subsection 72.1: Model formulation----
#' A zero-inflated binomial (ZIB)  model is given by the following expression.

#'  Success_i ~ ZIB(P_i, N_i, Pi)
#'  E(success_i)   = (1 - Pi) * mu_i
#'  var(sucess_i) = (1 - Pi) * (V_i + mu_i^2) - (1 - Pi)^2 * mu_i^2

#' Pi_i is the probability of success.
#' N_i is number of seeds in a plate 
#' Pi is the probability of a false zero.

#' And:
#'   mu_i = N_i * P_i  
#'   V_i  = N_i * P_i * (1 - P_i)



#' We will use covariates to model the P_ij
#'  logit(P_ij) = Covariate stuff


#' For the probability of a false zero, we use:
#'  logit(Pi) = Intercept

#' It is possible in glmmTMB to model Pi as a function of covariates.
#' Note: If Pi = 0, then we obtain the ordinary binomial GLMM.


#' One more time: 
#'  -Pi is for the zero-inflation stuff.
#'  -Pi is constant (at least that is what we will do)
#'  -P is for the binomial part. 
#'  -P is a function of covariates and the random effects.
#'  -The logistic link function is used for both terms.



#* Subsection 72.2: Execute the ZIB GLMM----

M4d <- glmmTMB(cbind(count,Failure) ~ ffun + ftemp + fpop + ftemp*fpop + ftemp*ffun + fpop*ffun,
               family = "binomial",
               ziformula =~ 1,
               data = s4)

#' The numerical output of the NB GLM is as follows.
summary(M4d)

glmmTMB:::Anova.glmmTMB(M4d, test.statistic = c("Chisq"))
#' both results are consistent

tempMeans <- emmeans(M4c, "ftemp", data = s2)
tempMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M4c, "ffun", data = s2)
funMeans
pairs(funMeans, comparisons = TRUE)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M4d)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M4d, test = "Chi")
#' The three interactions are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.
#' 

M4d1 <- glmmTMB(cbind(count,Failure) ~ ftemp*fpop + ftemp*ffun + fpop*ffun,
                family = "binomial",
                ziformula =~ 1,
                data = s4)

#' The numerical output of the NB GLM is as follows.
summary(M4d1)

glmmTMB:::Anova.glmmTMB(M4d1, test.statistic = c("Chisq"))
#' both results are consistent

tempMeans <- emmeans(M4d1, "ftemp", data = s2)
tempMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M4c, "ffun", data = s2)
funMeans
pairs(funMeans, comparisons = TRUE)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M4d1)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M4d1, test = "Chi")
#' The three interactions (ffun * ftemp ) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

#* Subsection 72.3 Assess overdispersion----  

#' Get the fitted values
F1 <- fitted(M4d1)

#' Get the Pearson residuals
E1 <- resid(M4d1, type = "pearson")

#' Let´s calculate pearson residuals ourselves
#' Get an X matrix
X <- model.matrix(~ ftemp*fpop + ftemp*ffun + fpop*ffun,
                  data = s4)
#' Or: X <- model.matrix(M2)

#' Get the betas for the binomial part.
beta.count <- fixef(M4d1)$cond

#' Calculate eta = X * beta
eta.count <- X %*% beta.count

#' Calculate P for the binomial part.
P <- exp(eta.count) / (1 + exp(eta.count))

#' Calculate Pi (probability of a false zero).
#' Get the estimated regression part for the binary part.
gamma  <- summary(M4d1)$coefficients$zi[1,"Estimate"] 

#' And this is Pi.
Pi <- exp(gamma) / (1 + exp(gamma))
Pi
#' Pi is 0.12

#' Recall that these are the mean and variance of a ZIB GLM(M):
#'  E(Y)   <- (1 - Pi) * mu
#'  var(Y) <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2 * mu^2

#'Where:
#'    mu = N * P
#'    V  = N * P * (1 - P) 
#'    N is the  number of seeds
#'    p is the probility of success of the binomial process.


#' We first calculate E(Y).
Ntrials <- s4$num
mu      <- Ntrials * P
ExpY    <- (1 - Pi) * mu

#' Now calculate V and var(Y).
V    <- Ntrials * P * (1 - P)
VarY <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2  * mu^2 


#' And here are the Pearson residuals.
PRes2 <- (s4$count - ExpY) / sqrt(VarY)


#' Assess overdispersion.
Npar <- length(fixef(M4d1)$cond) + length(fixef(M4d1)$zi) +1

sum(PRes2^2) / (nrow(s4) - Npar)
# minor overdispersion.

#' What does DHARMa tell us?
testDispersion(M4d1)
# 1.18

# Section 73: Model Validation ZIB GLMM----     

#* Subsection 73.1: Plot residuals vs fitted values and covariates----

# Residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = ExpY,
     y = PRes2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)            
#' better than the binomial GLMM!

#' What does DHARMa tell us?
E2BINqr <- simulateResiduals(fittedModel = M4d1, plot = FALSE)
plotQQunif(E2BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' better than the binomial GLMM!


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E2BINqr, quantreg = TRUE, smoothScatter = FALSE) 



# Residuals vs covariates in the model.
MyVar <- c("ffun")
s4$PRes2 <- PRes2
MyMultipanel.ggp2(Z = s4, 
                  varx = MyVar, 
                  vary = "PRes2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)


plotResiduals(E2BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E2BINqr, form = s4$ftemp, xlab = "Rank-transformed fun") 
plotResiduals(E2BINqr, form = s4$fpop, xlab = "Rank-transformed fun") 
#still acceptable as other parts of DHARMS still look good!

#* Subsection 73.2: Simulation for zero-inflation and one-inflation----

testZeroInflation(M4d1)
#' Pretty accurate and we can skip the simulate 1000 data onwards

#' Alternative if we want to simulate 1000 data sets
#' Simulate 1000 data sets from the ZIB GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M4d1, NSim, size = s2$num, seed = 12345)


#' We will extract the first column, and store it in a matrix Ysim.
Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Now we have 1,000 simulated data sets from the model. What shall we do 
#' with these simulated data sets? We can calculate the number of zeros in 
#' each of the 1,000 data sets. And we can also count how often we predict
#' a number of success equal to the total (which is a proportion of 1).

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s2$num)
}




#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0.1, 0.25),
     main = "Simulation results")
points(x = sum(s2$count == s2$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#' The red dot is the % of ones in the original data set.
#' The model predicts data sets that have just similar number of ones  



#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(300, 390))
points(x = sum(s2$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)
#' The red dot is the number of ones in the original data set.
#' The model predicts data sets that have similar number of zeros
#' as the observed data set.

#' Summarising: The ZIB GLMM can cope with the excessive number of zeros,
#' and also the excessive number of ones.

# Section 74: Visualise the results of the binomial GLMM----


#* Subsection 74.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun),
                      fpop = levels(s4$fpop),
                      ftemp = levels(s4$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ftemp*fpop + ftemp*ffun + fpop*ffun, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M4d1)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M4d1)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + facet_grid(ftemp~fpop, scales = "fixed") + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#' It is the fitted value for the interaction of ffun and ftemp + 95% 
#' interval. 

# Section 75: Beta glmm interpretation----

#' Execute the beta GLMM
M5d <- glmmTMB(PERT ~ ffun + ftemp + fpop + ftemp*fpop + ftemp*ffun + fpop*ffun,
               family = beta_family(link = "logit"), 
               data = s4)

#' Here is the numerical output
summary(M5d)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M5d)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M5d, test = "Chi")
#' The three interactions (ffun * ftemp and ffun*ftemp) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.


#' Execute the beta GLMM
M5d1 <- glmmTMB(PERT ~  ftemp*ffun + fpop*ffun,
                family = beta_family(link = "logit"), 
                data = s4)

#' Here is the numerical output
summary(M5d1)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M5d1)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M5d1, test = "Chi")
#' The three interactions are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

glmmTMB:::Anova.glmmTMB(M5d1, test.statistic = c("Chisq"))
tempMeans <- emmeans(M5c, "ftemp", data = s2)
tempMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M5c, "ffun", data = s2)
funMeans
pairs(funMeans, comparisons = TRUE)


# Section 76: Model validation----


#' Get the residuals and the fitted values.
E1 <- resid(M5d1, type = "pearson")
F1 <- fitted(M5d1)


#* Subsection 76.1: Plot residuals versus fitted values----
par(mfrow = c(1, 1))
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2, col = 1)


#* Subsection 76.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok?

#' Residuals vs temp.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' not okay

#' Residuals vs pop.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok

#* Subsection 76.3 Assess overdispersion----  

# Section 76: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M5d1, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M5d1)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' 0.21

#' We get the same message from DHARMa.
testDispersion(M5d1)
#'0.80 a bit underdispersion

#' What does DHARMa tell us?
E3BINqr <- simulateResiduals(fittedModel = M5d1, plot = FALSE)
plotQQunif(E3BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' passed the dispersion test

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E3BINqr, quantreg = TRUE, smoothScatter = FALSE) 


#' Standard model validation steps:
plotResiduals(E3BINqr, form = s4$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E3BINqr, form = s4$ftemp, xlab = "Rank-transformed fun") 
plotResiduals(E3BINqr, form = s4$fpop, xlab = "Rank-transformed fun") 
#' temp has a lot of within- variations 

# Section 65: Simulation----
testZeroInflation(M5d1)
#' not working with M5d1. Go stimulate


# The simulation function does also work for the beta GLMM:
Ysim <- simulate(M5d1, 1000)

#' Now we have 1000 simulated data set. Here is the first one:
Ysim[,1]

#' We can make a histogram of this one.
hist(Ysim[,1])

Interval <- seq(0, 1, by = 0.1)
SimData.int <- cut(Ysim[,1], breaks = Interval)
Tab <- table(SimData.int)
Tab

#' Here is a visualisation of this:
plot(Tab, type = "h", xlab = "Simulated values for coverage", ylab = "Frequency")

#' How is this for the observed data?
ObsData.int <- cut(s4$PERT, breaks = Interval)
TabObs      <- table(ObsData.int)
plot(Tab, type = "h", col = 2, xlab = "Values for coverage", ylab = "Frequency")

# What about doing this 1000 times?
Tab1000 <- matrix(NA, nrow = length(Interval)-1, ncol = 1000)
rownames(Tab1000) <- names(TabObs)
for (i in 1:1000){
  SimData.int <- cut(Ysim[,i], breaks = Interval)
  Tab1000[,i] <- table(SimData.int)
}

#' This is for the first simulated data set:  
Tab1000[,1]

#' This is for the second simulated data set:  
Tab1000[,2]

#' Etc.

MyData <- data.frame(SimulatedFreq = as.vector(Tab1000),
                     ID = rep(rownames(Tab1000), 1000))
MyData$ID <- factor(MyData$ID)
head(MyData,20)

#' Plot the residuals as a boxlot
boxplot(MyData$SimulatedFreq ~ MyData$ID)
points(1:10, TabObs, col = 2, cex = 2, pch = 16)


#' Instead of a boxplot, we can also plot this as a Cleveland dotplot.
#' First make a data frame for the frequencies for the observed coverage data.
DataOriginal <- data.frame(TabObs = TabObs,
                           ID = 1:10)
DataOriginal

#' And plot the simulated frequencies and the observed ones.
p <- ggplot()
p <- p + geom_point(data = MyData,
                    aes(x = SimulatedFreq,
                        y = ID),
                    col = "purple")
p <- p + geom_point(data = DataOriginal,
                    aes(x = TabObs,
                        y = ID),
                    col = "red",
                    size = 4)
p <- p + xlab("Frequeny") + xlab("p_portocorm development success (in bins")
p
#' The red dot is for the observed data. The purple dots are for the 
#' simulated data.
#' seems like the model is good to predict numbers


# Section 77: Visualise the results of the binomial GLMM----


#* Subsection 77.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun to the
#' floor.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ffun = levels(s4$ffun), 
                      ftemp = levels(s4$ftemp),
                      fpop = levels(s4$fpop))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun), summarize,
                stage = seq(from = min(s4$stage), 
                            to   = max(s4$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ftemp * ffun + fpop * ffun, 
                   data = MyData)
Xp


#'Extract parameters and parameter covariance matrix
Betas    <- fixef(M5d1)$cond
CovBetas <- vcov(M5d1)$cond


#'Calculate the fitted values on the predictor scale.
MyData$eta <- Xp %*% Betas
MyData$mu  <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' Calculate the SEs on the scale of the predictor function.
MyData$se    <- sqrt(diag(Xp %*% vcov(M5d1)$cond %*% t(Xp)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  - 1.96 *MyData$se))
head(MyData)

#' E. Plot everything
#' fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun),
                      position=position_dodge(width=1), alpha = 0.3)
p1 <- p1 + xlab("Fungi") + ylab("Propotion of protocorm formation")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = mu, colour = "red"), 
                      position=position_dodge(width=1))
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ffun, 
                             ymax = SeUp, 
                             ymin = SeLo, colour = "red"),
                         position=position_dodge(width=1))
p1 <- p1  + facet_grid(ftemp~fpop, scales = "fixed") + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#' It is the fitted value for the interaction of ffun, fpop, and ftemp + 95% 
#' interval. 

AIC(M5a, M5b, M5c, M5d1)
#'M5c has the lowest AIC value and the visulaise makes more sense
#'
################ only temp #######################################

# Section 78: Model formulation----

#' We will execute the following binomial GLMM

#' count_i ~ Binomial(Pi_i, num_i)
#' P_i us the number of success
#' E[success_i]   = num_i * Pi_i  
#' var[success_i] = num_i * Pi_i  * (1 - Pi_i)

#' count_i is the number of germinated seed on in plate i.


#'              exp(eta_i)
#' P_i = --------------------
#'           1 + exp(eta_i)

#' where
#'   eta_i = Fun_i * temp_i 

#' no random intercept.


# Section 79: Apply the binomial GLMM----

#' Execute the binomial GLMM using glmmTMB
M3e <- glmmTMB(cbind(count, Failure) ~ ftemp,
               data = s4,
               family = binomial)

#' Numerical outpout
summary(M3e)
#' No odd numbers.

step(M3e)
drop1(M3e, test = "Chi")

# Section 80: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M3e, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s4)
Npar <- length(fixef(M3e)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' Overdispersion!

#' We get the same message from DHARMa.
testDispersion(M3e)


# Section 81 Model Validation binomial GLMM----


#* Subsection 81.1: Plot residuals vs fitted values----

#' Get the Pearson residuals and the fitted values.
E1 <- resid(M3e, type = "pearson")
F1 <- fitted(M3e)


#' Plot residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
#' Ok pattern.

#* Subsection 82.2: Plot Pearson residuals versus covariates----
#' Plot Pearson residuals vs each covariate in the model, and not in the model.


#' Residuals vs fun.
s4$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' Ok?

#' What does DHARMa tell us?
#' Get scaled quantile residuals.
E1Binqr <- simulateResiduals(fittedModel = M3e, plot = FALSE)
#' warning message that we have a deviation from uniformity


#' Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(E1Binqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' That is a 'no'.


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E1Binqr, quantreg = TRUE, smoothScatter = FALSE) 
#' Problems for the larger residuals? And also for the smaller residuals
#' Those are the two clouds of residuals in our graph?
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)

plotResiduals(E1Binqr, form = s4$ftemp, xlab = "Rank-transformed fun") 
#' Trouble for the first and last? There is certainly a non-linear pattern
#' in the residuals with fungi, population, and temp.
#' some within-group residual distributions that are not uniform

#* Subsection 8.2: Simulation ---
testZeroInflation(M3e)

#' alternative please use the below codes
#' Simulate 1000 data sets from the binomial GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M1, NSim, size = s4$num, seed = 12345)


#' This is the first simulated data set:
YBin[,1]
#' The first column is the number of successes, and the second column is failure.
#' We will extract the first column, and store it in a matrix Ysim.

Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s4$num)
}


#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0, 1),
     main = "Simulation results")
points(x = sum(s4$count == s4$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)

#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(0, 950))
points(x = sum(s2$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)

#' What does DHARMa tell us?
testZeroInflation(M3c)
#' we cannot cope with the zeros


#' Expressed as some sort of p-value:
ZerosInData <- sum(s4$count == 0)
sum(Zeros > ZerosInData) / 1000



# Section 83: Visualise the results of the binomial GLMM----


#* Subsection 83.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus distance to the
#' floor, for average values of all other covariates.



#' We will not use the predict function, because we need to do
#' the predictions manually once we reach the zero-inflation
#' model.

#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ftemp = levels(s4$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ftemp, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M3e)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M3e)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' Back-standardise fun:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ftemp))
p1 <- p1 + xlab("Temp") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ftemp, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ftemp, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1
#this result is not accurate due to overdispersion

# Section 84: ZIB GLMM---- 


#* Subsection 84.1: Model formulation----
#' A zero-inflated binomial (ZIB)  model is given by the following expression.

#'  Success_i ~ ZIB(P_i, N_i, Pi)
#'  E(success_i)   = (1 - Pi) * mu_i
#'  var(sucess_i) = (1 - Pi) * (V_i + mu_i^2) - (1 - Pi)^2 * mu_i^2

#' Pi_i is the probability of success.
#' N_i is number of seeds in a plate 
#' Pi is the probability of a false zero.

#' And:
#'   mu_i = N_i * P_i  
#'   V_i  = N_i * P_i * (1 - P_i)



#' We will use covariates to model the P_ij
#'  logit(P_ij) = Covariate stuff


#' For the probability of a false zero, we use:
#'  logit(Pi) = Intercept

#' It is possible in glmmTMB to model Pi as a function of covariates.
#' Note: If Pi = 0, then we obtain the ordinary binomial GLMM.


#' One more time: 
#'  -Pi is for the zero-inflation stuff.
#'  -Pi is constant (at least that is what we will do)
#'  -P is for the binomial part. 
#'  -P is a function of covariates and the random effects.
#'  -The logistic link function is used for both terms.



#* Subsection 84.2: Execute the ZIB GLMM----

M4e <- glmmTMB(cbind(count,Failure) ~ ftemp,
               family = "binomial",
               ziformula =~ 1,
               data = s4)

#' The numerical output of the NB GLM is as follows.
summary(M4e)

glmmTMB:::Anova.glmmTMB(M4e, test.statistic = c("Chisq"))
#' both results are consistent

tempMeans <- emmeans(M4e, "ftemp", data = s2)
tempMeans
pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M4c, "ffun", data = s2)
funMeans
pairs(funMeans, comparisons = TRUE)

#' Model selection
# Use classical backwards model selection using the AIC:
step(M4e)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M4e, test = "Chi")
#' The three interactions (ffun * ftemp ) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

#* Subsection 84.3 Assess overdispersion----  

#' Get the fitted values
F1 <- fitted(M4e)

#' Get the Pearson residuals
E1 <- resid(M4e, type = "pearson")

#' Let´s calculate pearson residuals ourselves
#' Get an X matrix
X <- model.matrix(~ ftemp,
                  data = s4)
#' Or: X <- model.matrix(M2)

#' Get the betas for the binomial part.
beta.count <- fixef(M4e)$cond

#' Calculate eta = X * beta
eta.count <- X %*% beta.count

#' Calculate P for the binomial part.
P <- exp(eta.count) / (1 + exp(eta.count))

#' Calculate Pi (probability of a false zero).
#' Get the estimated regression part for the binary part.
gamma  <- summary(M4e)$coefficients$zi[1,"Estimate"] 

#' And this is Pi.
Pi <- exp(gamma) / (1 + exp(gamma))
Pi
#' Pi is 0.12

#' Recall that these are the mean and variance of a ZIB GLM(M):
#'  E(Y)   <- (1 - Pi) * mu
#'  var(Y) <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2 * mu^2

#'Where:
#'    mu = N * P
#'    V  = N * P * (1 - P) 
#'    N is the  number of seeds
#'    p is the probility of success of the binomial process.


#' We first calculate E(Y).
Ntrials <- s4$num
mu      <- Ntrials * P
ExpY    <- (1 - Pi) * mu

#' Now calculate V and var(Y).
V    <- Ntrials * P * (1 - P)
VarY <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2  * mu^2 


#' And here are the Pearson residuals.
PRes2 <- (s4$count - ExpY) / sqrt(VarY)


#' Assess overdispersion.
Npar <- length(fixef(M4e)$cond) + length(fixef(M4e)$zi) +1

sum(PRes2^2) / (nrow(s4) - Npar)
# minor overdispersion.

#' What does DHARMa tell us?
testDispersion(M4e)

# Section 85: Model Validation ZIB GLMM----     



#* Subsection 85.1: Plot residuals vs fitted values and covariates----

# Residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = ExpY,
     y = PRes2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)            
#' better than the binomial GLMM!

#' What does DHARMa tell us?
E2BINqr <- simulateResiduals(fittedModel = M4e, plot = FALSE)
plotQQunif(E2BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' better than the binomial GLMM!


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E2BINqr, quantreg = TRUE, smoothScatter = FALSE) 



# Residuals vs covariates in the model.
MyVar <- c("ffun")
s4$PRes2 <- PRes2
MyMultipanel.ggp2(Z = s4, 
                  varx = MyVar, 
                  vary = "PRes2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)


plotResiduals(E2BINqr, form = s4$ftemp, xlab = "Rank-transformed fun") 
#still acceptable as other parts of DHARMS still look good!

#* Subsection 85.2: Simulation for zero-inflation and one-inflation----

testZeroInflation(M4e)
#' Pretty accurate and we can skip the simulate 1000 data onwards

#' Alternative if we want to simulate 1000 data sets
#' Simulate 1000 data sets from the ZIB GLMM.
N    <- nrow(s4) #' Sample size
NSim <- 1000     #' Number of simulated data sets
YBin <- simulate(M2, NSim, size = s2$num, seed = 12345)


#' We will extract the first column, and store it in a matrix Ysim.
Ysim <- matrix(nrow=nrow(YBin), ncol = ncol(YBin))
for (i in 1:NSim){
  Num     <- YBin[,i]
  Ysim[,i] <- Num[,1] 
}

#' Now we have 1,000 simulated data sets from the model. What shall we do 
#' with these simulated data sets? We can calculate the number of zeros in 
#' each of the 1,000 data sets. And we can also count how often we predict
#' a number of success equal to the total (which is a proportion of 1).

#' Create two vectors to store the number of success equal to 0 and Total.
Zeros <- vector(length = NSim) #' For storing the number of success equal to 0.
Ones  <- vector(length = NSim) #' For storing the number of success equal to Total. 

#' Start a loop to determine the number of 0 and 1 in each simulated data set.
for(i in 1:NSim){
  Zeros[i] <- sum(Ysim[,i] == 0)
  Ones[i]  <- sum(Ysim[,i] == s2$num)
}




#' Let's plot the number of Ones in a histogram.
hist(Ones / N, 
     xlab = "Percentage of ones",
     ylab = "Number of ones in 1000 simulated data sets",
     xlim = c(0.1, 0.25),
     main = "Simulation results")
points(x = sum(s2$count == s2$num) / N, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#' The red dot is the % of ones in the original data set.
#' The model predicts data sets that have just similar number of ones  



#' Let's plot the number of zeros in a histogram.
plot(table(Zeros),
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(300, 390))
points(x = sum(s2$count == 0),
       y = 0,
       pch = 16,
       cex = 5,
       col = 2)
#' The red dot is the number of ones in the original data set.
#' The model predicts data sets that have similar number of zeros
#' as the observed data set.

#' Summarising: The ZIB GLMM can cope with the excessive number of zeros,
#' and also the excessive number of ones.

# Section 86: Visualise the results of the binomial GLMM----


#* Subsection 86.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun.


#' A. Specify covariate values for predictions
#' Create an artificial grid of covariate values

MyData <- expand.grid(ftemp = levels(s4$ftemp))
MyData

#' Alternative, we can use ddply to summarize the data with 
#' 1000 observations
MyData <- ddply(s4, 
                .(ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                            to   = max(s2$stage), 
                            length = 50))
#' remove column MyData$stage because it is not necessary in 
#' the following section
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)


#' B. Create X matrix with expand.grid

#' Get the corresponding X matrix.
Xp <- model.matrix(~ ftemp, 
                   data = MyData)
Xp


#' C. Calculate probability of success
Betas      <- fixef(M4e)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M4e)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))


#' E. Plot everything
#' temp:
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ftemp))
p1 <- p1 + xlab("Temp") + ylab("Probability of success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ftemp, 
                          y = P), 
                      colour = "red")
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = ftemp, 
                             ymax = SeUp, 
                             ymin = SeLo),
                         colour="red")
p1 <- p1 + scale_y_continuous(labels = label_number(accuracy=0.01))
p1
#' It is the fitted value for the interaction of ftemp + 95% 
#' interval. 

AIC(M5, M4a, M5a, M5b, M5c, M5d1)
#'M5c has the lowest AIC value and the visulaise makes more sense

anova(M5, M5a, M5b)
#M5a significnatly improved fit over M5

anova(M5, M5b)
#M5b is significantly improved fit over M5

anova(M5a, M5b)
#both does not have any effect over each others