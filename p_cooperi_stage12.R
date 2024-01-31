#  Section 1: Import the data and load the packages----
s2 <- read.csv("s2r.csv", header = TRUE)
s4 <- read.csv("s4r.csv", header = TRUE)
Pkgs2Load <- c("lattice", "ggplot2", "mgcv", "plyr", "MASS", 
               "pscl", "glmmTMB", "DHARMa", "cowplot", "rgl", 
               "GGally", "stringi", "performance", "emmeans", 
               "scales", "multcomp")
invisible(lapply(Pkgs2Load, library, character.only = TRUE))
# Follow outline Zuur et al. (2010) before running model
# Section 2: ZIB GLMM----
#' See what is binomial notes
Failure <- s2$num - s2$success1
M2 <- glmmTMB(cbind(success1,Failure) ~ fpop*ftemp*ffun,
              family = "binomial",
              ziformula =~ 1,
              data = s2)
glmmTMB:::Anova.glmmTMB(M2, test.statistic = c("Chisq"))
#* Subsection 2.1 Model validation--------  
F1 <- fitted(M2)
E1 <- resid(M2, type = "pearson")
X <- model.matrix(~ fpop*ftemp*ffun,
                  data = s2)
beta.count <- fixef(M2)$cond
eta.count <- X %*% beta.count
P <- exp(eta.count) / (1 + exp(eta.count))
gamma  <- summary(M2)$coefficients$zi[1,"Estimate"] 
Pi <- exp(gamma) / (1 + exp(gamma))
Pi
Ntrials <- s2$num
mu      <- Ntrials * P
ExpY    <- (1 - Pi) * mu
V    <- Ntrials * P * (1 - P)
VarY <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2  * mu^2 
PRes2 <- (s2$count - ExpY) / sqrt(VarY)
testDispersion(M2)
#' It's 0.95. Very close to 1. accept

# Section 3: Model Validation ZIB GLMM----     
#* Subsection 3.1: residuals vs fitted values and covariates----
# Residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = ExpY,
     y = PRes2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)            
#' better than the binomial GLMM!

#' Residuals vs fun.
s2$PRes2 <- PRes2
p <- ggplot()
p <- p + geom_boxplot(data = s2, 
                      aes(y = PRes2, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' ok

#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s2, 
                      aes(y = PRes2, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' okay

#' Residuals vs pop.
p <- ggplot()
p <- p + geom_boxplot(data = s2, 
                      aes(y = PRes2, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' OK
E2BINqr <- simulateResiduals(fittedModel = M2, plot = FALSE)
plotQQunif(E2BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' better than the binomial GLMM!
#' see class videos, notes, and zuur 2009. not to worry.

#' The scaled quantile residuals versus fitted values.
plotResiduals(E2BINqr, quantreg = TRUE, smoothScatter = FALSE) 
#still acceptable as other parts of DHARMS still look good! 
#see class videos and notes - the within variation is normally occurred. not to worry.
#* Subsection 3.2: Simulation for zero-inflation ----
testZeroInflation(M2)
#' Pretty accurate 
# Section 4: Visualise the results of the binomial GLMM----
MyData <- ddply(s2, 
                .(fpop, ffun, ftemp), summarize,
                stage = seq(from = min(s2$stage), 
                                  to   = max(s2$stage), 
                                  length = 50))
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)
Xp <- model.matrix(~ ffun * fpop * ftemp, 
                   data = MyData)
Xp
Betas      <- fixef(M2)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M2)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE))
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))

s2$Success= s2$count/s2$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s2, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Germination")
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

# Section 5: Beta glmm----
s4$PER <- s4$per / 100
#' See what is beta notes for the transformation.
N <- nrow(s4)
s4$PERT <- (s4$PER * (N - 1) + 0.5 ) / N
M5 <- glmmTMB(PERT ~ fpop*ffun*ftemp,
              family = beta_family(link = "logit"), 
              data = s4)
# Section 6: Model validation----
E1 <- resid(M5, type = "pearson")
F1 <- fitted(M5)
#* Subsection 6.1: residuals versus fitted values----
par(mfrow = c(1, 1))
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2, col = 1)
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
#' ok
#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' okay
#' Residuals vs pop.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' OK
# Section 6.2: overdispersion----
E1 <- resid(M5, type = "pearson")
N    <- nrow(s4)
Npar <- length(fixef(M5)$cond) + 1 
testDispersion(M5)
#'0.82. closed to 1. accept
#' see videos and notes on overdispersion and underdispersion
#'  DHARMa tell us?
E3BINqr <- simulateResiduals(fittedModel = M5, plot = FALSE)
plotQQunif(E3BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' okay. don't expect too much from DHARMa.
#' did residuals above
#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E3BINqr, quantreg = TRUE, smoothScatter = FALSE) 
#see class videos and notes - the within variation is normally occurred. not to worry.
# Section 7: Simulation Zero inflation----
Ysim <- simulate(M5, 1000)
Ysim[,1]
hist(Ysim[,1])
Interval <- seq(0, 1, by = 0.1)
SimData.int <- cut(Ysim[,1], breaks = Interval)
Tab <- table(SimData.int)
Tab
plot(Tab, type = "h", xlab = "Simulated values for coverage", ylab = "Frequency")
ObsData.int <- cut(s4$PERT, breaks = Interval)
TabObs      <- table(ObsData.int)
plot(Tab, type = "h", col = 2, xlab = "Values for coverage", ylab = "Frequency")
Tab1000 <- matrix(NA, nrow = length(Interval)-1, ncol = 1000)
rownames(Tab1000) <- names(TabObs)
for (i in 1:1000){
  SimData.int <- cut(Ysim[,i], breaks = Interval)
  Tab1000[,i] <- table(SimData.int)
}
Tab1000[,1]
MyData <- data.frame(SimulatedFreq = as.vector(Tab1000),
                     ID = rep(rownames(Tab1000), 1000))
MyData$ID <- factor(MyData$ID)
head(MyData,20)
DataOriginal <- data.frame(TabObs = TabObs,
                           ID = 1:10)
DataOriginal
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
#'  the model is good to predict numbers in zero inflation
# Section 8: Visualise the results of the beta GLMM (check the fit)----
MyData <- ddply(s4, 
                .(fpop, ffun, ftemp), summarize,
                stage = seq(from = min(s4$stage), 
                            to   = max(s4$stage), 
                            length = 50))
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)
Xp <- model.matrix(~ ffun * fpop * ftemp, 
                   data = MyData)
Xp
Betas    <- fixef(M5)$cond
CovBetas <- vcov(M5)$cond
MyData$eta <- Xp %*% Betas
MyData$mu  <- exp(MyData$eta) / (1 + exp(MyData$eta))
MyData$se    <- sqrt(diag(Xp %*% vcov(M5)$cond %*% t(Xp)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  - 1.96 *MyData$se))
View(MyData)
s4$Success = s4$count/s4$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s4, 
                      aes(y = Success, x = ffun))
p1 <- p1 + xlab("Fungi") + ylab("Protocorm")
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
#' interval. better than zero=inflated binomal and beta models. 