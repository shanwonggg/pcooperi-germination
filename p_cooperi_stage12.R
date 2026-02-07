#  Section 1: Import the data and load the packages----
s5 <- read.csv(file = "s5r.csv", 
               header = TRUE)
s4 <- read.csv("s4r.csv", header = TRUE)
Pkgs2Load <- c("lattice", "ggplot2", "mgcv", "plyr", "MASS", 
               "pscl", "glmmTMB", "DHARMa", "cowplot", "rgl", 
               "GGally", "stringi", "performance", "emmeans", 
               "scales", "multcomp", "viridis", "tidyverse", "hrbrthemes", "car")
invisible(lapply(Pkgs2Load, library, character.only = TRUE))

# Follow outline Zuur et al. (2010) before running model

# Section 2: ZIB GLMM---- (Germination)
s5$Failure <- s5$num - s5$total_count
M8 <- glmmTMB(cbind(total_count,Failure) ~ fpop*ftemp*ffun,
              family = "binomial",
              ziformula =~ 1,
              data = s5)
glmmTMB:::Anova.glmmTMB(M8, test.statistic = c("Chisq"))
#* Subsection 2.1 Model validation--------  
F1 <- fitted(M8)
E1 <- resid(M8, type = "pearson")
X <- model.matrix(~ fpop*ftemp*ffun,
                  data = s5)
beta.count <- fixef(M8)$cond
eta.count <- X %*% beta.count
P <- exp(eta.count) / (1 + exp(eta.count))
gamma  <- summary(M8)$coefficients$zi[1,"Estimate"] 
Pi <- exp(gamma) / (1 + exp(gamma))
Pi
Ntrials <- s5$num
mu      <- Ntrials * P
ExpY    <- (1 - Pi) * mu
V    <- Ntrials * P * (1 - P)
VarY <- (1 - Pi) * (V + mu^2) - (1 - Pi)^2  * mu^2 
PRes2 <- (s5$total_count - ExpY) / sqrt(VarY)
testDispersion(M8)

# Section 3: Model Validation ZIB GLMM----     
#* Subsection 3.1: residuals vs fitted values and covariates----
# Residuals vs fitted values.
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = ExpY,
     y = PRes2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)            
#' Residuals vs fun.
s5$PRes2 <- PRes2
p <- ggplot()
p <- p + geom_boxplot(data = s5, 
                      aes(y = PRes2, 
                          x = ffun))
p <- p + xlab("fun") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s5, 
                      aes(y = PRes2, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' Residuals vs pop.
p <- ggplot()
p <- p + geom_boxplot(data = s5, 
                      aes(y = PRes2, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
E2BINqr <- simulateResiduals(fittedModel = M8, plot = FALSE)
plotQQunif(E2BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' The scaled quantile residuals versus fitted values.
plotResiduals(E2BINqr, quantreg = TRUE, smoothScatter = FALSE) 

#* Subsection 3.2: Simulation for zero-inflation ----
testZeroInflation(M8)

# Section 4: Visualise the results of the binomial GLMM----
MyData <- ddply(s5, 
                .(fpop, ffun, ftemp), summarize,
                stage = seq(from = min(s5$stage), 
                            to   = max(s5$stage), 
                            length = 50))
MyData <- MyData[ , -which(names(MyData) %in% c("stage"))]
head(MyData)
Xp <- model.matrix(~ ffun * fpop * ftemp, 
                   data = MyData)
Xp
Betas      <- fixef(M8)$cond
MyData$eta <- Xp %*% Betas 
MyData$P   <- exp(MyData$eta) / (1 + exp(MyData$eta))
MyData$SE   <- sqrt(  diag(Xp %*% vcov(M8)$cond %*% t(Xp))  )
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE) / (1 + exp(MyData$eta + 1.96 * MyData$SE))
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE) / (1 + exp(MyData$eta - 1.96 * MyData$SE))

s5$success = s5$total_count/s5$num
p1 <- ggplot()
p1 <- p1 + geom_point(data = s5, 
                      aes(y = success, x = ffun))
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

# Section 5: Beta glmm----(Protocorm formation)
s4$PER <- s4$per / 100
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
#' Residuals vs temp.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = ftemp))
p <- p + xlab("temp") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 
#' Residuals vs pop.
p <- ggplot()
p <- p + geom_boxplot(data = s4, 
                      aes(y = E1, 
                          x = fpop))
p <- p + xlab("pop") + ylab("Pearson residuals")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_hline(yintercept = 0, lty = 2)
p 

# Section 6.2: overdispersion----
E1 <- resid(M5, type = "pearson")
N    <- nrow(s4)
Npar <- length(fixef(M5)$cond) + 1 
testDispersion(M5)
E3BINqr <- simulateResiduals(fittedModel = M5, plot = FALSE)
plotQQunif(E3BINqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(E3BINqr, quantreg = TRUE, smoothScatter = FALSE) 

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
