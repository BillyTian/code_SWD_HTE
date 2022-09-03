
library(nlme)

source("functions_calc_ss.R")
source("functions_gendata.R")
source("functions_empirical.R")

delta_OTE <- 0.3
mu_X <- 1

N <- c(rep(20, 24), rep(50, 24))
J <- rep(c(rep(3, 6), rep(4, 6), rep(5, 6), rep(6, 6)), 2)
alpha0 <- rep(c(rep(0.015, 3), rep(0.1, 3)), 8)
alpha1 <- rep(c(rep(0.01, 3), rep(0.05, 3)), 8)
rho0 <- rep(c(0.15, 0.3, 0.5), 16)
rho1 <- rep(c(0.1, 0.15, 0.3), 16)
table <- cbind(N, J, alpha0, alpha1, rho0, rho1)

I <- numeric(48)
nclusters_per_arm <- numeric(48)
pred.power <- numeric(48)
analytical_var <- numeric(48)

for (i in 1:nrow(table)){
  N.input <- as.numeric(table[i,][1])
  J.input <- as.numeric(table[i,][2])
  alpha0.input <- as.numeric(table[i,][3])
  alpha1.input <- as.numeric(table[i,][4])
  rho0.input <- as.numeric(table[i,][5])
  rho1.input <- as.numeric(table[i,][6])
  
  nclusters_per_arm[i] <- calc_CS_ATE(eff=delta_OTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, rho0=rho0.input, rho1=rho1.input)[1]
  I[i] <- calc_CS_ATE(eff=delta_OTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, rho0=rho0.input, rho1=rho1.input)[2]
  pred.power[i] <- calc_CS_ATE(eff=delta_OTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, rho0=rho0.input, rho1=rho1.input)[3]
  analytical_var[i] <- calc_CS_ATE(eff=delta_OTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, rho0=rho0.input, rho1=rho1.input)[4]
}

table <- cbind(table, nclusters_per_arm, I, analytical_var)

empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_CS_ATE(nullcase=T, parameter=table[i,]))
}

empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_CS_ATE(parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
