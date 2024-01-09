
# Section 1: Data description----

#' Question:
#' Do seeds germinate increased proportions when cultured with Cerato
#' and Tul OTUs obtained from the corresponding sites 
#' of seed collection under their native ambient 
#' temperature (15.5 0C)? 


#' I got the pure culture of the 4 fungi isolated from 
#' large and small populations by my colleagues. 
#' These 4 fungi are the dominant OMF (Tul) and less dominant OMF (Cerato) from FF, 
#' and the dominant OMF (Cerato) and less dominant OMF (Tul) from MX. 
#' I put seeds on oat-meal agar along with each fungi separately. 
#' I have a set of control plates without any fungi. 
#' Then I made up 10 replicates from each treatment. 
#' I put 50 – 100 seeds in each plate. 
#' For the incubate temperature, one complete set (of above design) incubated 
#' at Room Temperature; another complete set (of above design)
#' at 16C + 50% humidity. 


#' The data were collected in 124 days. 


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


#' Import the data.
s1 <- read.table(file = "d124_count_stage1.txt", 
                 header = TRUE,
                 na.strings = "NA",
                 stringsAsFactors = TRUE,
                 dec = ".")

#' Check the imported data.
dim(s1)
names(s1)
str(s1)


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
s1$success <- s1$count / s1$num 

#' remove samples that are contaminated.
s2<-subset(s1, num!=0)
dim(s2)



# Section 4: Data exploration----


#* Subsection 4.1: Missing values----

colSums(is.na(s2))
#' No missing value


#* Subsection 4.2: Zeros----
#' How many observations are equal to 0?
100 * sum(s2$per == 0) / nrow(s2)
#' 11% of zeros
plot(table(s1$success))
#' plenty of zeros


#* Subsection 4.3: Outliers----

#' Make a Cleveland dotplot for the response variable (per) and the
#' covariate fun, temp, and pop. There are three categorical covariates 
#' in the data set.
#' We rename them.

s2$ffun <- factor(s2$fun) 
table(s2$ffun)

s2$ftemp <- factor(s2$temp) 
table(s2$ftemp)

s2$fpop <- factor(s2$pop) 
table(s2$fpop)

table(s2$ffun, s2$ftemp)
table(s2$ftemp, s2$fpop)
table(s2$ffun, s2$fpop)
#' balanced for an interaction term


#* Subsection 4.4: Collinearity----
#' Make a pairplot of the covariates.
ToPlot <- c("ffun", "ftemp", "fpop")
ggpairs(s2[,ToPlot])
#' balanced. non collinearity

#* Subsection 4.5: Relationships----
#' The figure below shows a graph of the response variable
#' The success of seed germination versus each of the covariates (made with MyMultipanel.ggp2
#'  which is in our support file).
#'  Patterns in fungal species in different temp
p1 <- ggplot(s2, aes(x = ffun, y = per)) + 
  geom_boxplot() + xlab("Fungi")

p2 <- ggplot(s2, aes(x = ffun, y = per)) +
  geom_boxplot() + facet_wrap(~ftemp, ncol = 2)  + xlab("Fungi")

p3 <- ggplot(s2, aes(x = ffun, y = per)) +
  geom_boxplot() + facet_wrap(~fpop, ncol = 2)  + xlab("Fungi")

plot_grid(p1, p2, p3, ncol = 1, labels = c("A", "B", "C"))
#' No clearn patterns in GC, but a bit variable in RT and FF

#' We also show a scatterplot of Depth versus the number of dolphin 
#' sightings for each year. There seems to be a weak non-linear relationship. 

ggplot(s2, aes(x = ffun, y = per)) +
  geom_boxplot() +
  facet_wrap(~ftemp, ncol = 2, scales = "free_y")
#' in GC, higher seed germination is obtained in control
#' in RT, higher when there is pff_tul

ggplot(s2, aes(x = ffun, y = per)) +
  geom_boxplot() +
  facet_wrap(~fpop, ncol = 2, scales = "free_y")
#' in FF, higher seed germination is obtained with pff_tul
#' in MX, higher seed germination is obtained in control and pff_tul


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
s2$Failure <- s2$num - s2$count

#' Execute the binomial GLMM using glmmTMB
M1 <- glmmTMB(cbind(count, Failure) ~ ftemp*fpop*ffun,
              data = s2,
              family = binomial)

#' Numerical outpout
summary(M1)
#' No odd numbers.

# Section 7: Check overdispersion of the Binomial GLMM----

#' Get the Pearson residuals
E1 <- resid(M1, type = "pearson")

#' Get the sample size and the number of parameters.
N    <- nrow(s2)
Npar <- length(fixef(M1)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' Overdispersion!

#' We get the same message from DHARMa.
testDispersion(M1)


# Section 8 Model Validation binomial GLMM----


#* Subsection 8.1: Plot residuals vs fitted values----

#' Get the Pearson residuals and the fitted values.
E1 <- resid(M1, type = "pearson")
F1 <- fitted(M1)


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
s2$E1 <- E1
p <- ggplot()
p <- p + geom_boxplot(data = s2, 
                      aes(y = E1, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' Ok


#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s2, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok?

#' Residuals vs pop.
p <- ggplot()
p <- p + geom_boxplot(data = s2, 
                      aes(y = E1, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' OK

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

plotResiduals(E1Binqr, form = s2$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E1Binqr, form = s2$fpop, xlab = "Rank-transformed temp") 
plotResiduals(E1Binqr, form = s2$ftemp, xlab = "Rank-transformed pop") 
#' Trouble for the first and last? There is certainly a non-linear pattern
#' in the residuals with fungi, population, and temp.
#' some within-group residual distributions that are not uniform

#* Subsection 8.2: Simulation ---
testZeroInflation(M1)

#' alternative please use the below codes
#' Simulate 1000 data sets from the binomial GLMM.
N    <- nrow(s2) #' Sample size
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
s2$Success = s2$count/s2$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s2, 
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

M2 <- glmmTMB(cbind(count,Failure) ~ fpop*ftemp*ffun,
              family = "binomial",
              ziformula =~ 1,
              data = s2)

#' The numerical output of the NB GLM is as follows.
summary(M2)

#' The fungi, temp, 2-way (pop/fungi or temp/fungi), and 3-way 
#' interactions are significant at the 5% level, although two-way 
#' interaction is not significant.
#' The effect of fungi and temp alone is negative and significant at 
#' the 5% level. Also, the two-way (pop/fungi) has significantly 
#' negative effects

glmmTMB:::Anova.glmmTMB(M2, test.statistic = c("Chisq"))

allMeans <- emmeans(M2, specs = c("ffun", "ftemp", "fpop"))
allMeans
allmeans_pairs <- pairs(allMeans, comparisons = TRUE)
allmeans_pairs_tb <- as.data.frame(allmeans_pairs)
View(allmeans_pairs_tb)
allmeans_pairs_tb <- write.csv(allmeans_pairs_tb, file = "/Volumes/pcooperi_seed/allmeanpairs_stage1.csv")

emmeans(M2,pairwise ~ fpop*ffun*ftemp, type =" response", adjust="fdr")

tempMeans <- emmeans(M2, "ftemp", data = s2)
tempMeans
emmeans(M2, "ftemp", type =" response", adjust="fdr")

pairs(tempMeans, comparisons = TRUE)
funMeans <- emmeans(M2, "ffun", data = s2)
funMeans
summary(glht(M2, mcp(ffun="Tukey")))
summary(glht(M4, linfct = mcp(fpop = "Tukey")))



multcomp::cld(funMeans$emmeans,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
pairs(funMeans, comparisons = TRUE)
popMeans <- emmeans(M2, "fpop", data = s2)
popMeans
pairs(popMeans, comparisons = TRUE)


#' Model selection
# Use classical backwards model selection using the AIC:
step(M2)
# Backward selection indicates that a model with the interaction is the best
# We therefore execute the optimal model

drop1(M2, test = "Chi")
#' The three interactions (ffun * ftemp * fpop) are significantly
#' important as covariates. We do not want to drop this interaction. So, keep it.
#' we leave the model as it is, which is what we will do.

#* Subsection 10.3 Assess overdispersion----  

#' Get the fitted values
F1 <- fitted(M2)

#' Get the Pearson residuals
E1 <- resid(M2, type = "pearson")

#' Let´s calculate pearson residuals ourselves
#' Get an X matrix
X <- model.matrix(~ fpop*ftemp*ffun,
                  data = s2)
#' Or: X <- model.matrix(M2)

#' Get the betas for the binomial part.
beta.count <- fixef(M2)$cond

#' Calculate eta = X * beta
eta.count <- X %*% beta.count

#' Calculate P for the binomial part.
P <- exp(eta.count) / (1 + exp(eta.count))

#' Calculate Pi (probability of a false zero).
#' Get the estimated regression part for the binary part.
gamma  <- summary(M2)$coefficients$zi[1,"Estimate"] 

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
Ntrials <- s2$num
mu      <- Ntrials * P
ExpY    <- (1 - Pi) * mu

#' Now calculate V and var(Y).
V    <- Ntrials * P * (1 - P)
VarY <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2  * mu^2 


#' And here are the Pearson residuals.
PRes2 <- (s2$count - ExpY) / sqrt(VarY)


#' Assess overdispersion.
Npar <- length(fixef(M2)$cond) + length(fixef(M2)$zi) +1

sum(PRes2^2) / (nrow(s2) - Npar)
# No overdispersion; only minor underdispersion.



#' What does DHARMa tell us?
testDispersion(M2)




#' Get the sample size and the number of parameters.
N    <- nrow(s2)
Npar <- length(fixef(M2)$cond) + 1 #' The '+1' is for the sigma

#' Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - Npar)
DispersionStatistic
#' Overdispersion!

#' We get the different message from DHARMa.
testDispersion(M2)
#' It's 0.95. Very close to 1

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
E2BINqr <- simulateResiduals(fittedModel = M2, plot = FALSE)
plotQQunif(E2BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' better than the binomial GLMM!


#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E2BINqr, quantreg = TRUE, smoothScatter = FALSE) 



# Residuals vs covariates in the model.
MyVar <- c("ffun", "fpop",
           "ftemp")
s2$PRes2 <- PRes2
MyMultipanel.ggp2(Z = s2, 
                  varx = MyVar, 
                  vary = "PRes2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)


plotResiduals(E2BINqr, form = s2$ffun, xlab = "Rank-transformed fun") 
plotResiduals(E2BINqr, form = s2$fpop, xlab = "Rank-transformed temp") 
plotResiduals(E2BINqr, form = s2$ftemp, xlab = "Rank-transformed pop") 
#still acceptable as other parts of DHARMS still look good!

#* Subsection 11.2: Simulation for zero-inflation and one-inflation----

testZeroInflation(M2)
#' Pretty accurate and we can skip the simulate 1000 data onwards

#' Alternative if we want to simulate 1000 data sets
#' Simulate 1000 data sets from the ZIB GLMM.
N    <- nrow(s2) #' Sample size
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

# Section 12: Visualise the results of the binomial GLMM----


#* Subsection 12.1: Visualise the model fit----

#' Maybe we can learn something from the model fit. Let us
#' plot the probability of success versus ffun.


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
Betas      <- fixef(M2)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))

#' If you really want to use 'predict', then use this:
#' MyData$Nest <- "NA"
#' MyData$eta  <- predict(M1, newdata = MyData, re.form =~ 0, type = "link")
#' MyData$P    <- predict(M1, newdata = MyData, re.form =~ 0, type = "response")


#' D. Calculate standard errors (SE) for predicted values
#' SE of fitted values are given by the square root of
#' the diagonal elements of: X * cov(betas) * t(X)
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M2)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE)) 
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))
#write.csv(MyData, file = "/Volumes/pcooperi_seed/predict_stage1.csv")


#' E. Plot everything
#' fun:
s2$Success = s2$count/s2$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s2, 
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
#' interval. 
#' The right bottom (RT*MX*FUN) has the highest value and largest variations
#' In comparison of both temperatures, we found higher probability 
#' of seed germinated (success) in RT. 
#' With the interaction along with temperature and pop and fun, 
#' we found higher probability of seed germination in RT * MX with 
#' pff_tul and pev_tul

#' F. Plot everything in 1 graph_future
#' fun:
s2$Success = s2$count/s2$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s2, 
                      aes(y = Success, x = ffun, colour = ftemp, 
                          shape = fpop),
                      position=position_dodge(width=0.8), alpha = 0.3)
p1 <- p1 + xlab("Fungi") + ylab("Propotion of germination success")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = ffun, 
                          y = P, colour = ftemp,
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

#' G. Plot everything in 1 graph_future
#' fun:
s2$Success = s2$count/s2$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s2, 
                      aes(y = Success, x = interaction(ffun, fpop), colour = ftemp),
                      position=position_dodge(width=0.65), alpha = 0.3)
p1 <- p1 + xlab("") + ylab("Germination(%)")
p1 <- p1 + theme(text = element_text(size = 15)) 
p1 <- p1 + geom_point(data = MyData, 
                      aes(x = interaction(ffun, fpop), 
                          y = P, colour = ftemp), 
                      position=position_dodge(width=0.65))
p1 <- p1 + geom_errorbar(data = MyData, 
                         aes(x = interaction(ffun, fpop), 
                             ymax = SeUp, 
                             ymin = SeLo, colour = ftemp),
                         position=position_dodge(width=0.65))
p1 <- p1 + scale_y_continuous(labels = label_number(accuracy=0.01), expand = c(0, 0), limits = c(0, 0.90)) + 
  scale_color_manual(values=c("black", "black")) +
  theme(legend.position="none")
p1

#' H. Plot everything in 1 graph_bar
#' 
s2$Success = s2$count/s2$num
p1 <- ggplot(data = s2, aes(x = interaction(fpop, ffun), y = Success, fill = ftemp))
p1 <- p1 + geom_bar(position=position_dodge(width=0.89), stat = "identity")
p1 <- p1 + xlab("Fungi") + ylab("Propotion of germination success")
p1 <- p1 + theme(text = element_text(size = 15)) + 
  scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#4EBAD4", "#E64A35")) +
  expand_limits(y=c(0, 0.85))
p1

#' I. Plot everything in 1 box_plot 
#' when FF and MX together
fun_mean <- function(x){
  return(round(data.frame(y=mean(x),label=mean(x,na.rm=T)),digit=2))}

s2$Success = s2$count/s2$num
p1 <- ggplot(data = s2, aes(x = interaction(ftemp, fpop), y = Success, fill = ffun))
p1 <- p1 + geom_boxplot() 
p1 <- p1 + xlab("Fungi") + ylab("Propotion of germination success")
p1 <- p1 + theme(text = element_text(size = 15)) + 
  scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#ffffff", "#ffffff","#ffffff","#ffffff","#ffffff")) +
  expand_limits(y=c(0, 0.85)) + theme(legend.position="none") + stat_summary(fun.data = fun_mean, geom="text", vjust=-0.5, position = position_dodge(width=0.7))
p1
#control, SC pev_cerato, ST pev_tul, LC pff_cerato, LT pff_tul

#' J. statistical summary from boxplot without ggplot2
p1 + stat_summary(aes(label = round(..y..,2)), fun.y = "mean", geom = "text")
aggregate(Success ~ ffun + fpop + ftemp, data = s2, mean)
boxplot(s2$Success ~ s2$ffun + s2$pop + s2$temp, las = 2, par(mar = c(12,5,4,2)+0.1))

# filling with red and blue color in boxplot
# + scale_fill_manual(values=c("#4EBAD4", "#E64A35")) 

#* Subsection 12.2: summary table of seed germination----
#'look at the data summary in main effects
data_summary(s2, varname="Success", groupnames=c("fpop","ffun", "ftemp"))

#'conduct_pos hoc from each main effects
emmeans_M2 = emmeans(M2, ~fpop*ftemp*ffun)
pairs(tempMeans, comparisons = TRUE)

multcomp::cld(emmeans_M2,
          alpha=0.05,
          Letters=letters,
          adjust="tukey",
          decreasing = TRUE)
summary(glht(M2, mcp(ffun="Tukey")))

#' K comparison between models
AIC(M1, M2)
#M2 has a lower AIC overall and provides the best fit

