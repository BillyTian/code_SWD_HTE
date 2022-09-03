

library(nlme)

source("functions_calc_ss.R")
source("functions_gendata.R")
source("functions_empirical.R")

delta_HTE <- 0.2

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
analytical_var4 <- numeric(48)

for (i in 1:nrow(table)){
  N.input <- as.numeric(table[i,][1])
  J.input <- as.numeric(table[i,][2])
  alpha0.input <- as.numeric(table[i,][3])
  alpha1.input <- as.numeric(table[i,][4])
  rho0.input <- as.numeric(table[i,][5])
  rho1.input <- as.numeric(table[i,][6])
  
  nclusters_per_arm[i] <- calc_CS_HTE(eff=delta_HTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, rho0=rho0.input, rho1=rho1.input)[1]
  I[i] <- calc_CS_HTE(eff=delta_HTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, rho0=rho0.input, rho1=rho1.input)[2]
  pred.power[i] <- calc_CS_HTE(eff=delta_HTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, rho0=rho0.input, rho1=rho1.input)[3]
  analytical_var4[i] <- calc_CS_HTE(eff=delta_HTE, N=N.input, J=J.input, alpha0=alpha0.input, alpha1=alpha1.input, rho0=rho0.input, rho1=rho1.input)[4]
}

table <- cbind(table, nclusters_per_arm, I, analytical_var4)

empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_CS_HTE(nullcase=T, parameter=table[i,]))
}

empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_CS_HTE(parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))

















parameter <- table[7,]

I <- as.numeric(parameter[8])
J <- as.numeric(parameter[2])
N <- as.numeric(parameter[1])
alpha0 <- as.numeric(parameter[3])
alpha1 <- as.numeric(parameter[4])
rho0 <- as.numeric(parameter[5])
rho1 <- as.numeric(parameter[6])

beta2 <- 0.6
beta4 <- delta_HTE

set.seed(520)
simdata <- gendata(beta2=beta2, beta4=beta4, 
                   rho0=rho0, rho1=rho1, alpha0=alpha0, alpha1=alpha1, 
                   I=I, J=J, N=N)

mu1 <- mean(simdata[simdata$time==1,]$X)
mu2 <- mean(simdata[simdata$time==2,]$X)
mu3 <- mean(simdata[simdata$time==3,]$X)
mu4 <- mean(simdata[simdata$time==J,]$X)
D_mu <- diag(c(mu1, mu2, mu3, mu4))

eta1 <- mean((simdata[simdata$time==1,]$X)^2)
eta2 <- mean((simdata[simdata$time==2,]$X)^2)
eta3 <- mean((simdata[simdata$time==3,]$X)^2)
eta4 <- mean((simdata[simdata$time==J,]$X)^2)
D_eta <- diag(c(eta1, eta2, eta3, eta4))

otx1 <- c(0,0,0,0)

trtSeq <- matrix(0, ncol=J, nrow=J-1)
trtSeq[upper.tri(trtSeq)] <- 1
W_noN <- NULL
for (i in 1:(J-1)){
  W_noN <- c(W_noN, rep(trtSeq[i,], I/(J-1)))
}
W_noN_matrix <- matrix(W_noN, byrow=T, nrow=I)


Omega_diag <- NULL
for (i in 1:J){
  Omega_diag <- c(Omega_diag, var(W_noN_matrix[,i]))
}
Omega <- diag(Omega_diag, nrow=J)

B <- matrix(NA, nrow=I, ncol=J)
for (i in 1:I){
  for (j in 1:J){
    B[i,j] <- W_noN_matrix[i,j]-mean(W_noN_matrix[,j])
  }
}

Sigma <- matrix(0, nrow=J, ncol=J)
for (k in 1:J){
  for (l in 1:J){
    BB <- NULL
    for (i in 1:I){
      BB <- c(BB, (B[i,]%*%t(B[i,]))[k,l])
    }
    Sigma[k,l] <- mean(BB)
  }
}

Omega%*%D_mu
Sigma%*%D_mu

Omega%*%D_eta
Sigma%*%D_eta

Omega%*%D_eta
Sigma%*%D_eta


S_r1 <- as.matrix(cbind(diag(1, nrow=4), otx1, D_mu, otx1))
S_r2 <- as.matrix(c(otx1, sum(diag(Omega)), otx1, sum(diag(Omega%*%D_mu))))
S_r3 <- as.matrix(cbind(D_mu, otx1, D_eta, otx1))
S_r4 <- as.matrix(c(otx1, sum(diag(Omega%*%D_mu)), otx1, sum(diag(Omega%*%D_eta))))
S <- N*rbind(S_r1, t(S_r2), S_r3, t(S_r4))

cluster_period_mean <- NULL
for (i in 1:(I*J)){
  cluster_period_mean <- c(cluster_period_mean, 
                         mean(simdata[simdata$T_id==i,]$X))
}
time <- rep(seq(1,J,1), I)
cpdata <- data.frame(cbind(cluster_period_mean, time))

tau1 <- mean((cpdata[cpdata$time==1,]$cluster_period_mean)^2)
tau2 <- mean((cpdata[cpdata$time==2,]$cluster_period_mean)^2)
tau3 <- mean((cpdata[cpdata$time==3,]$cluster_period_mean)^2)
tau4 <- mean((cpdata[cpdata$time==4,]$cluster_period_mean)^2)
D_tau <- diag(c(tau1, tau2, tau3, tau4))

M_r1 <- as.matrix(cbind(diag(1, nrow=4), otx1, D_mu, otx1))
M_r2 <- as.matrix(c(otx1, sum(diag(Omega)), otx1, sum(diag(Omega%*%D_mu))))
M_r3 <- as.matrix(cbind(D_mu, otx1, D_tau, otx1))
M_r4 <- as.matrix(c(otx1, sum(diag(Omega%*%D_mu)), otx1, sum(diag(Omega%*%D_tau))))
M <- N^2*rbind(M_r1, t(M_r2), M_r3, t(M_r4))


cpdata$I_id <- rep(1:I, each=J)

Lambda <- matrix(0, nrow=J, ncol=J)
for (j in 1:J){
  for (k in 1:J){
    CC <- NULL
    for (i in 1:I){
      C <- as.matrix(cpdata[cpdata$I_id==i,]$cluster_period_mean)
      CC <- cbind(CC, (C%*%t(C))[j,k])
    }
    Lambda[j,k] <- mean(CC)
  }
}





ones <- matrix(1, nrow=J, ncol=J)
Q_r1 <- as.matrix(cbind(ones, otx1, D_mu, otx1))
Q_r2 <- as.matrix(c(otx1, sum(diag(Sigma%*%ones)), otx1, t(rep(1,J))%*%Sigma%*%D_mu%*%rep(1,J)))
Q_r3 <- as.matrix(cbind(D_mu%*%ones, otx1, Lambda, otx1))
Q_r4 <- as.matrix(c(otx1, t(rep(1,J))%*%Sigma%*%D_mu%*%rep(1,J), otx1, sum(diag(Sigma%*%Lambda))))
Q <- N^2*rbind(Q_r1, t(Q_r2), Q_r3, t(Q_r4))

alpha0 <- 0.015
alpha1 <- 0.01
lambda1 <- 1-alpha0
lambda2 <- 1+(N-1)*alpha0-N*alpha1
lambda3 <- 1+(N-1)*alpha0+(J-1)*alpha1

c <- 1/lambda1
d <- -(lambda2-lambda1)/(N*lambda1*lambda2)
h <- -(lambda3-lambda2)/(J*N*lambda2*lambda3)
U <- c*S+d*M+h*Q
U_inv <- solve(U)
U_inv_lr <- U_inv[6:10, 6:10]
U_inv_lr

1/N* t(rep(1,J))%*%((c+d*N)*Omega+h*N*Sigma)%*%rep(1,J)/ ( sum(diag(Sigma%*%(c*D_eta+d*N*D_tau+h*N*Lambda)))*t(rep(1,J))%*%((c+d*N)*Omega+h*N*Sigma)%*%rep(1,J) - (t(rep(1,J))%*%((c+d*N)*Omega+h*N*Sigma)%*%D_mu%*%rep(1,J))^2 )



########################################################################
D_mu <- diag(1, nrow=4)
D_tau <- diag(2, nrow=4)
D_eta <- diag(3, nrow=4)
Omega <- matrix(1.5, nrow=4, ncol=4)
Sigma <- matrix(2.5, nrow=4, ncol=4)
Lambda <- matrix(3.5, nrow=4, ncol=4)

1/N* t(rep(1,J))%*%((c+d*N)*Omega+h*N*Sigma)%*%rep(1,J)/ ( sum(diag(Sigma%*%(c*D_eta+d*N*D_tau+h*N*Lambda)))*t(rep(1,J))%*%((c+d*N)*Omega+h*N*Sigma)%*%rep(1,J) - (t(rep(1,J))%*%((c+d*N)*Omega+h*N*Sigma)%*%D_mu%*%rep(1,J))^2 )

c <- 1/lambda1
d <- -(lambda2-lambda1)/(N*lambda1*lambda2)
h <- -(lambda3-lambda2)/(J*N*lambda2*lambda3)

F1 <- (c+d*N)*diag(1, nrow=4)+h*N*matrix(1, ncol=4, nrow=4)
F2 <- as.matrix(rep(0,J))
F3 <- t(F2)
F4 <- t(as.matrix(rep(1,J))) %*% ((c+d*N)*Omega+h*N*Sigma) %*% as.matrix(rep(1,J))
F_matrix <- rbind(cbind(F1, F2), cbind(F3, F4))

G1 <- F1 %*% D_mu
G2 <- F2
G3 <- t(F2)
G4 <- t(as.matrix(rep(1,J))) %*% ((c+d*N)*Omega+h*N*Sigma) %*% D_mu %*% as.matrix(rep(1,J))
G_matrix <- rbind(cbind(G1, G2), cbind(G3, G4))

H1 <- c*D_eta + d*N*D_tau + h*N*Lambda
H2 <- F2
H3 <- t(F2)
H4 <- sum(diag( Sigma %*% (c*D_eta + d*N*D_tau + h*N*Lambda) )) 
H_matrix <- rbind(cbind(H1, H2), cbind(H3, H4))

U <- N*rbind(cbind(F_matrix, G_matrix), cbind(t(G_matrix), H_matrix))

solve(U)
########################################################################