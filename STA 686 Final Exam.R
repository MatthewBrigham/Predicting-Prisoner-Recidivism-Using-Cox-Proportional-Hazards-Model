#Matthew Brigham
#STA 686
#Spring 2021
#Final Exam

library(survival)
library(Hmisc) #for linearity assumption with splines
library(rms)

#Import Data
setwd("/Users/Matt/Library/Mobile Documents/com~apple~CloudDocs/Cleveland State/Courses/STA 686 - Reliability/Exams/Final Exam")
dat = read.csv("RECID.csv")
head(dat)
str(dat)
################################################################################
################################################################################
################################################################################
#####   Problem 1
################################################################################
################################################################################
################################################################################

#Fit Univariate PH Model
S.dat<-Surv(dat$WEEK, dat$ARREST)
plot(S.dat, xlab = "Weeks to arrest", ylab = "Probability of Survival", main = "Survival Curve")

#FIN

      #Part a - curves converge - assumption met
      km.fin = survfit(S.dat ~ FIN, data = dat)
      plot(km.fin,fun="cloglog",conf.int="none",ylab="Complementary Log-Log",xlab="Time to Arrest (Weeks)", main = "Complementary log-log Plot for FIN")
      
      #Part b
      cox.fin<-coxph(S.dat~FIN,data=dat,robust=T)
      cox.fin
      cox.fin$var  #variance-covariance estimate of logHR
      sqrt(cox.fin$var) #std error of logHR
      
      #Part c
      hr.fin = exp(cox.fin$coefficients); hr.fin  #estimate HR
      exp(confint.default(cox.fin))  #95% CI for HR
      100*(hr.fin-1)  # % change in risk
      

#RACE

      #Part a - curves converge? - assumption met?
      km.race = survfit(S.dat ~ RACE, data = dat)
      plot(km.race,fun="cloglog",conf.int="none",ylab="Complementary Log-Log",xlab="Time to Arrest (Weeks)", main = "Complementary log-log Plot for RACE")
      
      #Part b
      cox.race<-coxph(S.dat~RACE,data=dat,robust=T)
      cox.race
      cox.race$var  #variance-covariance estimate of logHR
      sqrt(cox.race$var) #std error of logHR
      
      #Part c
      hr.race = exp(cox.race$coefficients); hr.race  #estimate HR
      exp(confint.default(cox.race))  #95% CI for HR
      100*(hr.race-1)  # % change in risk

      
################################################################################
################################################################################
################################################################################
#####   Problem 2
################################################################################
################################################################################
################################################################################
      
#AGE
      
      #Part A - KM Curve
      cutpoints.age <- quantile(dat$AGE, probs=c(0,0.5,1)) # Dichotomize age
      dat$AGE.cut <- cut(dat$AGE, cutpoints.age, include.lowest=T) # Assign each inmate to a specific category
      
      km.age <- survfit(S.dat ~ AGE.cut, data = dat)
      plot(km.age, conf.int = F, ylab = "Probability of Survival", xlab="Time to Arrest in Weeks", main = "KM Estimate AGE")      

      #Part B - PH Assumption - Converge - assumption met
      plot(km.age,fun="cloglog",conf.int="none",ylab="Complementary Log-Log",xlab="Time to Arrest (Weeks)", main="Complementary log-log Plot for AGE")
      
      #Part C - Univariate Model
      cox.age <- coxph(S.dat ~ AGE, data = dat, robust = T)
      cox.age
      cox.age$var      
      sqrt(cox.age$var)      

      #Part D
      hr.age = exp(cox.age$coefficients); hr.age  #estimate HR
      exp(confint.default(cox.age))  #95% CI for HR
      100*(hr.age-1)  # % change in risk
      
      #Part E - Nonlinearity Assumption
      rcspline.plot(x = as.numeric(dat$AGE), y = as.numeric(dat$ARREST), event = as.numeric(dat$WEEK), model="cox", nk=3, ylim=c(-5,5))
      cox.age.ass = cph(S.dat ~ rcs(AGE,3), data=dat, x = TRUE, y = TRUE)
      cox.age.ass <- robcov(cox.age.ass)
      anova(cox.age.ass)
      
#PRIO
      
      #Part A - KM Curve
      cutpoints.prio <- quantile(dat$PRIO, probs=c(0,0.5,1)) # Dichotomize age
      dat$PRIO.cut <- cut(dat$PRIO, cutpoints.prio, include.lowest=T) # Assign each inmate to a specific category
      
      km.prio <- survfit(S.dat ~ PRIO.cut, data = dat)
      plot(km.prio, conf.int = F, ylab = "Probability of Survival", xlab="Time to Arrest in Weeks", main = "KM Estimate PRIO")      
      
      #Part B - PH Assumption - Converge - assumption met
      plot(km.prio,fun="cloglog",conf.int="none",ylab="Complementary Log-Log",xlab="Time to Arrest (Weeks)", main="Complementary log-log Plot for PRIO")
      
      #Part C - Univariate Model
      cox.prio <- coxph(S.dat ~ PRIO, data = dat, robust = T)
      cox.prio
      cox.prio$var      
      sqrt(cox.prio$var)      
      
      #Part D
      hr.prio = exp(cox.prio$coefficients); hr.prio  #estimate HR
      exp(confint.default(cox.prio))  #95% CI for HR
      100*(hr.prio-1)  # % change in risk
      
      #Part E - Nonlinearity Assumption
      rcspline.plot(x = as.numeric(dat$PRIO), y = as.numeric(dat$ARREST), event = as.numeric(dat$WEEK), model="cox", nk=3, ylim=c(-5,5))
      cox.prio.ass = cph(S.dat ~ rcs(PRIO,3), data=dat, x = TRUE, y = TRUE)
      cox.prio.ass <- robcov(cox.prio.ass)
      anova(cox.prio.ass)
      
      
################################################################################
################################################################################
################################################################################
#####   Problem 3
################################################################################
################################################################################
################################################################################
      
      #Fit full effect model - using 7 covariates
      cox.full <- cph(S.dat ~ FIN+AGE+RACE+WEXP+MAR+PARO+PRIO, data=dat, x=TRUE, y=TRUE)
      cox.full = robcov(cox.full)
      
      #Part A
      vif(cox.full) 
      
      #Part B
      exp(cox.full$coefficients) #HR for each covariate
      exp(confint.default(cox.full)) #95% CI
      cox.full
      
################################################################################
################################################################################
################################################################################
#####   Problem 4
################################################################################
################################################################################
################################################################################
      
      #Fit reduced model
      cox.reduced <- cph(S.dat ~ FIN+AGE+PRIO, data=dat, x=TRUE, y=TRUE, surv = T)
      cox.reduced
      #Likelihood Ratio Test - p-value = 0.3772 > 0.05 use reduced (same)
      lrtest(cox.full, cox.reduced)
      
################################################################################
################################################################################
################################################################################
#####   Problem 5
################################################################################
################################################################################
################################################################################
      
      #Use cox.reduced
      
      #Parts A and B
      summary(dat$FIN)
      summary(dat$AGE)
      summary(dat$PRIO)
      
        #Categorical:
        finplot.4a = survplot(cox.reduced, FIN = c(0,1), AGE = 23, PRIO = 2); finplot.4a #FIN
        title(main="Estimated Adjusted Survival Curves")
        legend(40, 0.35, c("FIN=0", "FIN=1"), lty = c(1,2))
        
        # #Continuous
        # ageplot.4a = survplot(cox.reduced, FIN = 1, AGE = 23, PRIO = 2); ageplot.4a #AGE
        # abline(v=c(13, 26, 39))
        # abline(h=c(0.958, 0.885, 0.825))
        # 
        # prioplot.4a = survplot(cox.reduced, FIN = 1, AGE = 23, PRIO = 2); prioplot.4a #PRIO
        # abline(v=c(13, 26, 39))
        # abline(h=c(0.968, 0.915, 0.867))
        
      #Part B
      fin0plot.4b = survplot(cox.reduced, FIN = 0, AGE = 23, PRIO = 2)
      fin0plot.4b  
      abline(v=c(13, 26, 39))
      abline(h=c(0.953, 0.875, 0.802))
      title(main="Predicted Adjusted Survival Curves for FIN = 0")
      
      fin1plot.4b = survplot(cox.reduced, FIN = 1, AGE = 23, PRIO = 2)
      fin1plot.4b  
      abline(v=c(13, 26, 39))
      abline(h=c(0.965, 0.908, 0.855))
      title(main="Predicted Adjusted Survival Curves for FIN = 1")
      
      # ageplot.4b = survplot(cox.reduced, FIN = 1, AGE = 23, PRIO = 2); 
      # ageplot.4b
      # abline(v=c(13, 26, 39))
      # abline(h=c(0.965, 0.905, 0.855))
      # 
      # prioplot.4b = survplot(cox.reduced, FIN = 1, AGE = 23, PRIO = 2) 
      # prioplot.4b 
      # abline(v=c(13, 26, 39))
      # abline(h=c(0.968, 0.905, 0.855))
      
################################################################################
################################################################################
################################################################################
#####   Problem 6
################################################################################
################################################################################
################################################################################
      
      #Part A - Is AGE a confounder?  No - coeffiecient is not significant
      lm.age.confound <- lm(AGE ~ FIN, data=dat)
      summary(lm.age.confound)
      
      #Part B - Is PRIO a confounder?  No - coeffiecient is not significant
      lm.prio.confound <- lm(PRIO ~ FIN, data=dat)
      summary(lm.prio.confound)
      
################################################################################
################################################################################
################################################################################
#####   Problem 7
################################################################################
################################################################################
################################################################################
      
      #Part A - FIN-AGE Interaction Plot suggests potential interaction, but anova says no interaction. Go with anova
      cox.finage.int <- cph(S.dat ~ FIN*AGE, data = dat, x=TRUE, y=TRUE)
      cox.finage.int <- robcov(cox.finage.int)
      cox.finage.int
      anova(cox.finage.int)
      
      exp(cox.finage.int$coefficients[2]) #estimate HR for AGE given no FIN
      
      exp(cox.finage.int$coefficients[2] + cox.finage.int$coefficients[3]) #estimate HR for AGE given FIN
      
      # Plot the predicted log(HR) for treatment (FIN = 0 and 1) over 
      # AGE. This is the best way to illustrate an interaction.
      
      pred.agefin.int <- predict(cox.finage.int, newdata = dat, se.fit=T, conf.int=0.95)
      index0 <- which(dat$FIN==0)
      age0 <- dat$AGE[index0]
      age1 <- dat$AGE[-index0]
      plot(age0, pred.agefin.int$linear.predictors[index0], type='l', xlab="AGE", ylab="Log HR")
      lines(age1,pred.agefin.int$linear.predictors[-index0], col="red", lty=2)
      legend(35, 0.5, c("FIN=0", "FIN=1"), lty = c(1,1), col=1:2, text.col = 1:2)
      title("Interaction plot between FIN and AGE")
      
      #Part B - FIN-PRIO - graph suggest possible interaction except it isnt statistically significant
      cox.finprio.int <- cph(S.dat ~ FIN*PRIO, data = dat, x=TRUE, y=TRUE)
      cox.finprio.int<-robcov(cox.finprio.int)
      cox.finprio.int
      anova(cox.finprio.int)
      
      exp(cox.finprio.int$coefficients[2]) #estimate HR for PRIO given no FIN
      
      exp(cox.finprio.int$coefficients[2] + cox.finprio.int$coefficients[3]) #estimate HR for PRIO given FIN
      
      # Plot the predicted log(HR) for treatment (FIN = 0 and 1) over 
      # PRIO. This is the best way to illustrate an interaction.
      
      pred.finprio.int <- predict(cox.finprio.int, newdata=dat, se.fit=T, conf.int=0.95)
      index0 <- which(dat$FIN==0)
      prio0 <- dat$PRIO[index0]
      prio1 <- dat$PRIO[-index0]
      plot(prio0, pred.finprio.int$linear.predictors[index0], type='l', xlab="PRIO", ylab="Log HR")
      lines(prio1,pred.finprio.int$linear.predictors[-index0], col="red", lty=2)
      legend(0, 2, c("FIN=0", "FIN=1"), lty = c(1,1), col=1:2, text.col = 1:2)
      title("Interaction plot between FIN and PRIO")
      
################################################################################
################################################################################
################################################################################
#####   Problem 8
################################################################################
################################################################################
################################################################################
      
      #Part A - not much difference in Dxy index.orig and incex.corrected - good thing, low prediction error
      validate(cox.reduced, B=100, dxy=T)
      
      #Part B 
      cal3 = calibrate(cox.reduced, B=100, u=13, m=100) #past 3 months
      plot(cal3)
      title(main="Calibration Plot for 3 Months")
      
      #Part C
      cal9 = calibrate(cox.reduced, B=100, u=39, m=100) #past 9 months
      plot(cal9)
      title(main="Calibration Plot for 9 Months")
      
      #Part D
      #No evidence of overfitting
        
      