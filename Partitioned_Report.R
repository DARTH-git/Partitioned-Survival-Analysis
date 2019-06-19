rm(list=ls())

if (!require(survival)) install.packages('survival', repos = "http://cran.us.r-project.org"); library(survival)
if (!require(flexsurv)) install.packages('flexsurv', repos = "http://cran.us.r-project.org"); library(flexsurv)
if (!require(MASS)) install.packages('MASS',         repos = "http://cran.us.r-project.org"); library(MASS)
if (!require(DAAG)) install.packages('DAAG',         repos = "http://cran.us.r-project.org"); library(DAAG)
if (!require(knitr)) install.packages('knitr',       repos = "http://cran.us.r-project.org"); library(knitr)


source("SurvFunctions_final.R")


data      <- read.csv('data.cens.csv', header = T)
head_data <- head(data)
kable(head_data)


KM.pfs <- survfit(Surv(time = pfs, event = event.pfs) ~ 1, data = data) # Kaplan-Meier fit for PFS
plot(KM.pfs, main = paste("Progression-Free Survival"), ylab = "Survival probability", xlab = "Years") # plot KM for PFS


KM.os  <- survfit(Surv(time = os,  event = event.os)  ~ 1, data = data) # Kaplan-Meier fit for OS
plot( KM.os,  main = paste("Overall Survival"),   ylab = "Survival probability", xlab = "Years") # plot KM for OS


v.n            <- c("Stable","Progressed","Dead") # state names
n.s            <- length(v.n)    # number of states
n.t            <- 60             # number of cycles to run
c.l            <- 1 / 12         # cycle length (a month)
n.psa          <- 10000          # number of PSA simulations
# Costs and utilities  
c.S    <- 400                     # cost of remaining one cycle stable
c.P    <- 100                     # cost of remaining one cycle progressed
c.D    <- 0                       # cost of remaining one cycle dead
u.S    <- 0.8                     # utility when stable 
u.P    <- 0.5                     # utility when progressed
u.D    <- 0                       # utility when dead
d.r    <- 0.03                    # discount rate per cycle, same for costs and effectiveness
times  <- seq(from = 0, to = n.t, by = c.l)  # sequence of times to be considered in the model

v.dwc  <- 1 / (1 + d.r) ^ (times) # calculate discount weights for costs for each cycle based on discount rate d.r
v.dwe  <- 1 / (1 + d.r) ^ (times) # calculate discount weights for effectiveness for each cycle based on discount rate d.r


fit.pfs  <- fit.fun(time = "pfs", event = "event.pfs", data = data) 


fit.os  <- fit.fun( time = "os", event = "event.os", data = data)  


best.pfs <- fit.pfs$Weibull   
plot(KM.pfs, ylab = "Survival Probability", xlab = "Time",  main = paste ("True vs Fitted PFS")) 
lines(best.pfs,  col = 2, t = times, lty = 2)     


best.os <- fit.os$Weibull 
plot(KM.os, ylab = "Survival Probability", xlab = "Time",  main = c("True vs Fitted OS"))
lines(best.os,  col = 2, t = times, lty = 2)    


surv <-partsurv(best.pfs, best.os, title = "all", time= times)

M.tr <- as.matrix(surv$trace)                       
matplot(M.tr, type = 'l', lty=1)                      
legend("right", v.n, col=1:n.s, lty=rep(1,n.s), bty='n')


# Mean Costs and QALYs per cycle
v.tc <- M.tr %*% c(c.S, c.P, c.D)  # calculate expected costs by multiplying m.M with the cost vector for the different health states   
v.tu <- M.tr %*% c(u.S, u.P, u.D)  # calculate expected QALYs by multiplying m.M with the utilities for the different health states  

# Discounted Mean Costs and QALYs
v.tc.d <-  t(v.tc) %*% v.dwc   # Discount costs by multiplying the cost vector with discount weights (v.dwc) 
v.te.d <-  t(v.tu) %*% v.dwe *c.l   # Discount QALYS by multiplying the QALYs vector with discount weights (v.dwe)

results <- data.frame( "Total Discounted Cost" = v.tc.d, 
                       "Total Discounted QALYs" = v.te.d, 
                       check.names = F)
kable(results)


# generate correlated values from a mutlivariate normal distribution for PFS and OS.
psa.pfs <- normboot.flexsurvreg(x = best.pfs,n.psa)
psa.os  <- normboot.flexsurvreg(x = best.os,n.psa)

# vectors to store the total costs and total utilities from each simulation
TC <- TE <- c()

for (i in 1:n.psa) {
  # at each simulation, draw a sample of the parameters (shape and scale) from the best-fitting survivor distribution for PSA and OS
  pfs.surv <- best.pfs$dfns$p(times, psa.pfs[i,1], psa.pfs[i,2], lower.tail = F)
  os.surv  <- best.pfs$dfns$p(times, psa.os[i,1],   psa.os[i,2], lower.tail = F)
  
  prog            <- os.surv - pfs.surv          # estimate the probability of remaining in the progressed state
  prog[prog < 0]  <- 0                           # in cases where the probability is negative replace with zero
  s               <- pfs.surv                    # probability of remaining stable
  d               <- 1 - os.surv                 # probability of being dead
  
  M.tr <- as.matrix(data.frame(S = s, P = prog, D = d))
  
  # Compute Cost-Effectiveness Outcomes
  # Mean Costs and QALYs per cycle
  v.tc <- M.tr %*% c(c.S, c.P, c.D)  # calculate expected costs by multiplying m.M with the cost vector for the different health states   
  v.tu <- M.tr %*% c(u.S, u.P, u.D)  # calculate expected QALYs by multiplying m.M with the utilities for the different health states  
  
  # Discounted Mean Costs and QALYs
  v.tc.d <-  t(v.tc) %*% v.dwc   # Discount costs  by multiplying the cost vector with discount weights (v.dwc) 
  v.te.d <-  t(v.tu) %*% v.dwe * c.l   # Discount QALYS  by multiplying the QALYs vector with discount weights (v.dwe)
  
  TC <- c(TC, v.tc.d)
  TE <- c(TE, v.te.d)
}


results <- data.frame( "Mean Total Discounted Cost" = mean(TC), 
                       "Mean Total Discounted QALYs" = mean(TE), 
                       check.names = F)

kable(results)

