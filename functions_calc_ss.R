#Function to calculate predicted number of clusters to get at least (1-beta) power, the actual power that can be reached, and the analytical variance 
#under the cross-sectional design, HTE test
calc_CS_HTE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, rho0, rho1, alpha=0.05, beta=0.2){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # rho0: within-period covariate ICC
  # rho1: between-period covariate ICC
  # alpha: type I error rate
  # beta: type II error rate
  # 
  # Output:
  # nca: number of clusters per arm (distinct sequence)
  # I: predicted number of clusters
  # Var4: analytical variance for HTE
  # power: actual predicted power
  
  nca <- 0
  power <- 0
  while (power < 1-beta){
    nca <- nca + 1
    I <- (J-1)*nca
    
    #make the staggered table
    trtSeq <- matrix(0, ncol=J, nrow=J-1)
    trtSeq[upper.tri(trtSeq)] <- 1
    
    W_noN <- NULL
    for (i in 1:(J-1)){
      W_noN <- c(W_noN, rep(trtSeq[i,], nca))
    }
    W_noN_matrix <- matrix(W_noN, byrow=T, nrow=I)
    
    #Follow the definition in Li et al. 2018
    #W_ij, cluster-period as the unit
    U <- sum(W_noN)
    W <- sum(colSums(W_noN_matrix)^2)
    V <- sum(rowSums(W_noN_matrix)^2)
    
    lambda1 <- 1-alpha0
    lambda2 <- 1+(N-1)*alpha0-N*alpha1
    lambda3 <- 1+(N-1)*alpha0+(J-1)*N*alpha1
    zeta1 <- 1-rho0
    zeta2 <- 1+(N-1)*rho0-N*rho1
    zeta3 <- 1+(N-1)*rho0+(J-1)*N*rho1
    
    Var4 <- sigma2y/sigma2x * I*J^2/
      ( (I*U-W)*J*(J*(N-1)*zeta1/lambda1 + (J-1)*zeta2/lambda2 + zeta3/lambda3) + (U^2+I*J*U-J*W-I*V)*(1/lambda2-1/lambda3)*(zeta3-zeta2) )
    
    power <- pnorm( sqrt(eff^2/Var4)-qnorm(1-alpha/2) )
  }
  return(c(nca, I, power, Var4))
}


#Function to calculate predicted number of clusters to get at least (1-beta) power, the actual power that can be reached, and the analytical variance
#under the closed-cohort design, HTE test
calc_CC_HTE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, alpha2, rho0, alpha=0.05, beta=0.2){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # rho0: within-period covariate ICC
  # alpha: type I error rate
  # beta: type II error rate
  # 
  # Output:
  # nca: number of clusters per arm (distinct sequence)
  # I: predicted number of clusters
  # Var4: analytical variance for HTE
  # power: actual predicted power
  
  nca <- 0
  power <- 0
  while (power < 1-beta){
    nca <- nca + 1
    I <- (J-1)*nca

    trtSeq <- matrix(0, ncol=J, nrow=J-1)
    trtSeq[upper.tri(trtSeq)] <- 1
    
    W_noN <- NULL
    for (i in 1:(J-1)){
      W_noN <- c(W_noN, rep(trtSeq[i,], nca))
    }
    W_noN_matrix <- matrix(W_noN, byrow=T, nrow=I)

    U <- sum(W_noN)
    W <- sum(colSums(W_noN_matrix)^2)
    V <- sum(rowSums(W_noN_matrix)^2)
    
    lambda1 <- 1-alpha0+alpha1-alpha2
    lambda2 <- 1-alpha0-(J-1)*(alpha1-alpha2)
    lambda3 <- 1+(N-1)*(alpha0-alpha1)-alpha2
    lambda4 <- 1+(N-1)*alpha0+(J-1)*(N-1)*alpha1+(J-1)*alpha2
    
    #zeta1 <- 0
    zeta2 <- J*(1-rho0)
    #zeta3 <- 0
    zeta4 <- J*(1+(N-1)*rho0)
    
    kappa1 <- (N-1)*zeta2/lambda2 + zeta4/lambda4
    kappa3 <- (alpha2+(N-1)*alpha1)*(1+(N-1)*rho0)/(lambda3*lambda4) + (N-1)*(alpha2-alpha1)*(1-rho0)/(lambda1*lambda2) 
    
    Var4 <- sigma2y/sigma2x * (I*J)/
      ( (I*U-W)*kappa1 + (U^2+I*J*U-J*W-I*V)*J*kappa3 )
    
    power <- pnorm( sqrt(eff^2/Var4)-qnorm(1-alpha/2) )
  }
  return(c(nca, I, power, Var4))
}


#Function to calculate predicted number of clusters to get at least (1-beta) power, the actual power that can be reached, and the analytical variance
#under the cross-sectional design, ATE test
calc_CS_ATE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, rho0, rho1, alpha=0.05, beta=0.2){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # rho0: within-period covariate ICC
  # rho1: between-period covariate ICC
  # alpha: type I error rate
  # beta: type II error rate
  # 
  # Output:
  # nca: number of clusters per arm (distinct sequence)
  # I: predicted number of clusters
  # Var: analytical variance for ATE
  # power: actual predicted power
  
  nca <- 1
  power <- 0
  while (power < 1-beta){
    nca <- nca + 1
    I <- (J-1)*nca
    
    #make the staggered table
    trtSeq <- matrix(0, ncol=J, nrow=J-1)
    trtSeq[upper.tri(trtSeq)] <- 1
    
    W_noN <- NULL
    for (i in 1:(J-1)){
      W_noN <- c(W_noN, rep(trtSeq[i,], nca))
    }
    W_noN_matrix <- matrix(W_noN, byrow=T, nrow=I)
    
    U <- sum(W_noN)
    W <- sum(colSums(W_noN_matrix)^2)
    V <- sum(rowSums(W_noN_matrix)^2)
    
    lambda1 <- 1-alpha0
    lambda2 <- 1+(N-1)*alpha0-N*alpha1
    lambda3 <- 1+(N-1)*alpha0+(J-1)*N*alpha1
    
    Var <- sigma2y/N * (I*J*lambda2*lambda3)/
      ( (U^2+I*J*U-J*W-I*V)*lambda3-(U^2-I*V)*lambda2 )
    
    power <- pt(qt(1-alpha/2, I-2), I-2, ncp=eff/sqrt(Var), lower.tail = F)
  }
  return(c(nca, I, power, Var))
}

#Function to calculate predicted number of clusters to get at least (1-beta) power, the actual power that can be reached, and the analytical variance
#under the closed-cohort design, ATE test
calc_CC_ATE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, alpha2, rho0, alpha=0.05, beta=0.2){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # rho0: within-period covariate ICC
  # alpha: type I error rate
  # beta: type II error rate
  # 
  # Output:
  # nca: number of clusters per arm (distinct sequence)
  # I: predicted number of clusters
  # Var: analytical variance for ATE
  # power: actual predicted power  
  
  nca <- 1
  power <- 0
  while (power < 1-beta){
    nca <- nca + 1
    I <- (J-1)*nca

    trtSeq <- matrix(0, ncol=J, nrow=J-1)
    trtSeq[upper.tri(trtSeq)] <- 1
    
    W_noN <- NULL
    for (i in 1:(J-1)){
      W_noN <- c(W_noN, rep(trtSeq[i,], nca))
    }
    W_noN_matrix <- matrix(W_noN, byrow=T, nrow=I)
    
    U <- sum(W_noN)
    W <- sum(colSums(W_noN_matrix)^2)
    V <- sum(rowSums(W_noN_matrix)^2)
    
    lambda1 <- 1-alpha0+alpha1-alpha2
    lambda2 <- 1-alpha0-(J-1)*(alpha1-alpha2)
    lambda3 <- 1+(N-1)*(alpha0-alpha1)-alpha2
    lambda4 <- 1+(N-1)*alpha0+(J-1)*(N-1)*alpha1+(J-1)*alpha2
    
    Var <- sigma2y/N * (I*J*lambda3*lambda4)/( (U^2+I*J*U-J*W-I*V)*lambda4-(U^2-I*V)*lambda3 )
    
    power <- pt(qt(1-alpha/2, I-2), I-2, ncp=eff/sqrt(Var), lower.tail = F)

  }
  return(c(nca, I, power, Var))
}
