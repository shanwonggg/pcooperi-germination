
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
#'  16C and 23C
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

#* Subsection 18.2: summary table of protocorm----
#'look at the data summary in main effects
data_summary(s4, varname="Success", groupnames=c("fpop","ffun", "ftemp"))

#'conduct_pos hoc from each main effects
emmeans_M5 = emmeans(M5, ~fpop*ftemp*ffun)
pairs(tempMeans, comparisons = TRUE)

multcomp::cld(emmeans_M5,
          alpha=0.05,
          Letters=letters,
          adjust="tukey",
          decreasing = TRUE)
summary(glht(M5, mcp(ffun="Tukey")))

#' K comparison between models
AIC(M3, M4, M5)
#M5 has a lower AIC overall and provides the best fit
