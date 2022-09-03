
#Function to generate the data with the cross-sectional design
gendata_CS <- function(beta2, beta4, rho0, rho1, alpha0, alpha1, sigma2x=1, sigma2y=1, I, J, N){
  # Argument:
  # beta2: true parameter of the main effect of treatment
  # beta4: true parameter of the interaction effect of treatment and the univariate covariate
  # rho0: within-period covariate ICC
  # rho1: between-period covariate ICC
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # sigma2x: total variance of the single covariate X
  # sigma2y: total variance of outcome Y
  # I: number of clusters 
  # J: number of periods in each cluster
  # N: number of individuals in each period
  # 
  # Output:
  # data: simulated data frame
  
  nca <- I/(J-1)
  N_id <- seq(1, (I*J*N), 1)
  T_id <- rep(1:(I*J), each=N)
  I_id <- rep(1:I, each=N*J)
  time <- rep(rep(1:J, each=N), I)
  data <- data.frame(cbind(N_id, T_id, I_id, time))
  
  #Generate beta1j and beta3j
  beta1j <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75)[1:J]
  beta3j <- c(0.2, 0.4, 0.6, 0.8, 1.6, 3.2)[1:J]
  data$beta1 <- rep(rep(beta1j, each=N), I)
  data$beta3 <- rep(rep(beta3j, each=N), I)
  
  #Generate the treatment assignment indicators based on design parameters
  trtSeq <- matrix(0, ncol=J, nrow=J-1)
  trtSeq[upper.tri(trtSeq)] <- 1
  W_assignment <- NULL
  for (i in 1:(J-1)){
    W_assignment <- c(W_assignment, rep(trtSeq[i,], nca))
  }
  data$W <- rep(W_assignment, each=N)
  
  #Generate X (one continuous individual-level covariate)
  Xmean <- 1
  c <- rnorm(I*J*N, 0, sqrt(sigma2x*(1-rho0)))
  b <- rep(rnorm(I*J, 0, sqrt(sigma2x*(rho0-rho1))), each=N)
  a <- rep(rnorm(I, 0, sqrt(sigma2x*rho1)), each=N*J)
  data$X <- Xmean+a+b+c
  
  #Generate Y under null and alternative
  epsilon <- rnorm(I*J*N, 0, sqrt(sigma2y*(1-alpha0)))
  u <- rep(rnorm(I*J, 0, sqrt(sigma2y*(alpha0-alpha1))), each=N)
  gamma <- rep(rnorm(I, 0, sqrt(sigma2y*alpha1)), each=N*J)
  
  data$Y <- data$beta1 + beta2*data$W + data$beta3*data$X + beta4*data$W*data$X + gamma+u+epsilon

  return(data)
}

#Function to generate the data with the closed-cohort design
gendata_CC <- function(beta2, beta4, rho0, alpha0, alpha1, alpha2, sigma2x=1, sigma2y=1, I, J, N){
  # Argument:
  # beta2: true parameter of the main effect of treatment
  # beta4: true parameter of the interaction effect of treatment and the univariate covariate
  # rho0: within-period covariate ICC
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # sigma2x: total variance of the single covariate X
  # sigma2y: total variance of outcome Y
  # I: number of clusters 
  # J: number of periods in each cluster
  # N: number of individuals in each period
  # 
  # Output:
  # data: simulated data frame
  
  nca <- I/(J-1)
  N_id <- seq(1, (I*J*N), 1)
  T_id <- rep(1:(I*J), each=N)
  I_id <- rep(1:I, each=N*J)
  time <- rep(rep(1:J, each=N), I)
  
  #To fit random effect v, create v_id
  ik_matrix <- matrix(1:(N*I), nrow=I, ncol=N, byrow=T)
  v_id_matrix <- NULL
  for (i in 1:J){
    v_id_matrix <- cbind(v_id_matrix, ik_matrix)
  }
  v_id <- NULL
  for (i in 1:I){
    v_id <- c(v_id, v_id_matrix[i,])
  }
  
  data <- data.frame(cbind(N_id, T_id, I_id, time, v_id))
  
  #Generate beta1j and beta3j
  beta1j <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75)[1:J]
  beta3j <- c(0.2, 0.4, 0.6, 0.8, 1.6, 3.2)[1:J]
  data$beta1 <- rep(rep(beta1j, each=N), I)
  data$beta3 <- rep(rep(beta3j, each=N), I)
  
  #Generate the treatment assignment indicators based on design parameters
  trtSeq <- matrix(0, ncol=J, nrow=J-1)
  trtSeq[upper.tri(trtSeq)] <- 1
  W_assignment <- NULL
  for (i in 1:(J-1)){
    W_assignment <- c(W_assignment, rep(trtSeq[i,], nca))
  }
  data$W <- rep(W_assignment, each=N)
  
  #Generate X (one continuous individual-level covariate)
  Xmean <- 1
  
  #Create a_i, repeat for periods
  a <- rep(rnorm(I, 0, sqrt(sigma2x*rho0)), each=N)
  c <- rnorm(I*N, 0, sqrt(sigma2x*(1-rho0)))
  X_oneperiod <- Xmean+a+c
  X_matrix <- matrix(X_oneperiod, ncol=N, nrow=I, byrow=T)
  
  X_vector <- NULL
  for (i in 1:I){
    X_vector <- c(X_vector, rep(X_matrix[i,], J))
  }
  data$X <- X_vector
  
  #Generate Y under null and alternative
  #Create v_ik
  v_matrix <- matrix(rnorm(I*N, 0, sqrt(sigma2y*(alpha2-alpha1))),
                     nrow=I, ncol=N, byrow=T)
  v_big_matrix <- NULL
  for (i in 1:J){
    v_big_matrix <- cbind(v_big_matrix, v_matrix)
  }
  v <- NULL
  for (i in 1:I){
    v <- c(v, v_big_matrix[i,])
  }
  
  epsilon <- rnorm(I*J*N, 0, sqrt(sigma2y*(1-alpha0+alpha1-alpha2)))
  u <- rep(rnorm(I*J, 0, sqrt(sigma2y*(alpha0-alpha1))), each=N)
  gamma <- rep(rnorm(I, 0, sqrt(sigma2y*alpha1)), each=N*J)
  
  data$Y <- data$beta1 + beta2*data$W + data$beta3*data$X + beta4*data$W*data$X + gamma+u+epsilon+v
  
  return(data)
}