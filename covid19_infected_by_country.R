#install.packages("ggplot2")
#library(ggplot2)
#install.packages("tidyverse")
install.packages("purrr")

library(tidyr)

#install.packages("deSolve")

library(purrr)

library(dplyr)

#theme_minimal()

# Germany
# Infected <- c(16, 18, 21, 26, 53, 66, 117, 150, 188, 240, 349, 534, 684, 847, 1112, 1460, 1884, 2369, 3062, 3795, 4838, 6012, 7156, 8198, 10999)

#China
#Infected <- c(45, 62, 121, 198, 291, 440, 571, 830, 1287, 1975, 2744, 4515)

# UK
Infected <- c(2,2,3,3,4,8,8,9,13,13,19,23,35,40,51,85,114,160,206,271,321,373,456,590,797,1061,1391,1543,1950,2626,3269)

# SW
# Infected <- c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,7,7,12,14,15,21,35,94,101,161,203,248,355,500,599,814,961,1022,1103,1190,1279,1439,1639,1763,1934,2046,2286,2526,2840,3069,3447,3700,4028,4435,4947,5568,6131,6443,6830,7206,7693,8419,8419)

# FR
# Infected <- c(0,0,2,3,3,3,4,5,5,5,6,6,6,6,6,6,6,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,14,18,38,57,100,130,191,204,285,377,653,949,1126,1209,1784,2281,2281,3661,4469,4499,6633,7652,9043,10871,12612,14282,16018,19856,22304,25233,29155,32964,37575,40174,44550,52128,56989,59105,64338,89953,92839,98010,109069,112950,112950)

# RW
#Infected <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,5,7,8,8,17,17,19,36,40,41,50,54,60,70,70,75,82,84,89,102,104,105,105,110,110)
Day <- 1:(length(Infected)) # length of the infection period
N <- 66000000 # pupulation of the UK
#N <- 10230000 # pupulation of the SW
#N <- 66990000 # pupulation of the FR
#N <- 12210000 # pupulation of the RW

old <- par(mfrow = c(1, 2))
#plot(Day, Infected, type ="b") # bubbles
plot(Day, Infected, type ="b")
plot(Day, Infected, log = "y") # type b by default
abline(lm(log10(Infected) ~ Day))
title("Confirmed COVID-19 cases in UK, lin & log scales", outer = TRUE, line = -2)

## Let's extrapolate to maybe 90 days, via epidemiologic modeling
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

#SIR <- function(t, y, parms) { # SIR <- function(time, state, parameters), with time=t, state=y
#  beta <- parms[1]
#  gamma <- parms[2]
#  S <- y[1]
#  I <- y[2]
#  R <- y[3]
#  return(list(c(S = -beta/N * S * I, I = beta/N * S * I - gamma * I, gamma * I)))
#}

library(deSolve)
init <- c(S = N-Infected[1], I = Infected[1], R = 0)
RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  sum((Infected - fit)^2)
}

# general optimazion nethod based on Nelder-Mead
# https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
Opt <- optim(c(0.5, 0.5), RSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1)) # optimize with some sensible conditions
Opt$message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
# OR
## [1] "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
print(Opt_par)
#      beta     gamma 
# 0.6235803 0.3764198 
# with beta, β (infection rate) = rate of contacts per time unit (say day) an infected individual meets a specific other person in the population
# and 
# and base reproduction rate R_0 = beta*N/gamma
# and incubation (or infection) period gamma = constant recovery rate γ
# where γ = 1/D, with D the incubation period (days)

## Note. these results for beta and gamma differ highly from what is expected theoritically
# i.e beta <- 4.5e-07 and gamma <- 1/5 from previous studies
# sources: https://staff.math.su.se/hoehle/blog/2020/03/16/flatteningthecurve.html

# R0 for the beta and gamma values
# Altogether, this leads to an R0 of 2.25, which roughly corresponds to the R0 of SARS-CoV-2
#R0 <- Opt_par["beta"]*N/Opt_par["gamma"]
R0 <- Opt_par["beta"]/Opt_par["gamma"] # R0/N
print(R0)
# 1.656609


# Grid where to evaluate
max_time <- 150 # here predictions up to 150 days, i.e + (length(Day)-max_time) beyond length(Day)
times <- seq(0, max_time, by=0.01)

# time in days, prediction to up to 90 days (highly unlikely fact will confirm the hypothesis though)
#fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))
fit <- data.frame(ode(y = init, times = times, func = SIR, parms = Opt_par))
col <- 1:3 # colour
print(fit)


#sink('fit0_09042020.txt'); fit; sink()
#write.csv(fit, "fit0_09042020.csv")

#matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 2, lty = 1, col = col, log = "y")
matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 2, lty = 1, col = col)
matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 2, lty = 1, col = col, log = "y")
#title("Infected Cases 2019-nCoV in UK (worst case), SIR model", outer = TRUE, line = -2)

points(Day, Infected)
legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 1, col = col, inset = 0.05)
title("Predicted 2019-nCoV in UK (worst case), SIR", outer = TRUE, line = -2)


# Function to compute the final size.
s_inf <- function(R0) {
  f_target <- function(x) { x - exp(-R0*(1-x)) }
  result <- uniroot(f_target, lower=1e-12, upper=1-1e-12)$root
  return(result)
}

# Final proportion of infected.
1 - s_inf(R0)
# [1] 0.6709214

## We can use the above equation to verify that the larger $R_0$, the larger is the final size of the outbreak:
R0_grid <- c(1.25, 1.5, 1.75,2, 2.25, 2.5, 3)
map_dbl( R0_grid, ~ 1-s_inf(.x)) %>% setNames(R0_grid) %>% scales::percent(accuracy=1)
#1.25   1.5  1.75     2  2.25   2.5     3 
#"37%" "58%" "71%" "80%" "85%" "89%" "94%"

