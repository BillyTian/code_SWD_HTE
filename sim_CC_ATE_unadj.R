library(lme4)

source("functions_calc_ss.R")
source("functions_gendata.R")
source("functions_empirical.R")

delta_OTE <- 0.3
mu_X <- 1

b3j <- c(0.2, 0.4, 0.6, 0.8, 1.6, 3.2)
beta4 <- 0.1

N <- c(rep(20, 16), rep(50, 16))
J <- rep(c(rep(3, 4), rep(4, 4), rep(5, 4), rep(6, 4)), 2)
alpha0 <- rep(c(rep(0.015, 2), rep(0.1, 2)), 8)
alpha1 <- rep(c(rep(0.01, 2), rep(0.05, 2)), 8)
alpha2 <- rep(c(rep(0.2, 2), rep(0.5, 2)), 8)
rho0 <- rep(c(0.2, 0.5), 16)
table <- cbind(N, J, alpha0, alpha1, alpha2, rho0)

I <- numeric(32)
nclusters_per_arm <- numeric(32)
pred.power <- numeric(32)
analytical_var <- numeric(32)

for (i in 1:nrow(table)){
  N.input <- as.numeric(table[i,][1])
  J.input <- as.numeric(table[i,][2])
  alpha0.input <- as.numeric(table[i,][3])
  alpha1.input <- as.numeric(table[i,][4])
  alpha2.input <- as.numeric(table[i,][5])
  rho0.input <- as.numeric(table[i,][6])
  
  nclusters_per_arm[i] <- calc_CC_ATE(eff=delta_OTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, alpha2=alpha2.input, rho0=rho0.input)[1]
  I[i] <- calc_CC_ATE(eff=delta_OTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, alpha2=alpha2.input, rho0=rho0.input)[2]
  pred.power[i] <- calc_CC_ATE(eff=delta_OTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, alpha2=alpha2.input, rho0=rho0.input)[3]
  analytical_var[i] <- calc_CC_ATE(eff=delta_OTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, alpha2=alpha2.input, rho0=rho0.input)[4]
}

table <- cbind(table, nclusters_per_arm, I, analytical_var)

empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_CC_ATE_unadj(nullcase=T, parameter=table[i,]))
}

empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_CC_ATE_unadj(parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))