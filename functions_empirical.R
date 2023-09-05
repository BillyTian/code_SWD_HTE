
#Function to compute empirical power or empirical type I error under cross-sectional design, HTE test 
empirical_CS_HTE <- function(nullcase=F, parameter, nsims=5000){
  # Argument:
  # nullcase: calculate empirical type I error or power. The former if nullcase=TRUE
  # parameter: the combination of design parameters of interest 
  # nsims: number of data replications
  
  I <- as.numeric(parameter[8])
  J <- as.numeric(parameter[2])
  N <- as.numeric(parameter[1])
  alpha0 <- as.numeric(parameter[3])
  alpha1 <- as.numeric(parameter[4])
  rho0 <- as.numeric(parameter[5])
  rho1 <- as.numeric(parameter[6])
  
  beta2 <- 0.6
  beta4 <- delta_HTE
  
  if (nullcase==T){
    beta4 <- 0
  }
  
  pvalue <- NULL
  count <- NULL
  
  for (i in 1:nsims){
    set.seed(0402+2022*i)
    simdata <- gendata_CS(beta2=beta2, beta4=beta4, rho0=rho0, rho1=rho1, alpha0=alpha0, alpha1=alpha1, I=I, J=J, N=N)
    simdata$time <- factor(simdata$time)
    
    fit <- try(lme(Y ~ -1 + time + W + time:X + W:X, data=simdata, random=list(I_id = ~ 1, T_id = ~ 1)), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    
    R <- c(rep(0,2*J+1),1)
    beta <- fit$coef$fixed
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 1-pchisq(test.stat, 1)
  }
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  
  return(c(empirical, error.rate))
}


#Function to compute empirical power or empirical type I error under closed-cohort design, HTE test 
empirical_CC_HTE <- function(nullcase=F, parameter, nsims=5000){
  # Argument:
  # nullcase: calculate empirical type I error or power. The former if nullcase=TRUE
  # parameter: the combination of design parameters of interest 
  # nsims: number of data replications
  
  I <- as.numeric(parameter[8])
  J <- as.numeric(parameter[2])
  N <- as.numeric(parameter[1])
  alpha0 <- as.numeric(parameter[3])
  alpha1 <- as.numeric(parameter[4])
  alpha2 <- as.numeric(parameter[5])
  rho0 <- as.numeric(parameter[6])
  
  beta2 <- 0.6
  beta4 <- delta_HTE
  
  if (nullcase==T){
    beta4 <- 0
  }
  
  pvalue <- NULL
  count <- NULL
  for (i in 1:nsims){
    set.seed(0402+2022*i)
    simdata <- gendata_CC(beta2=beta2, beta4=beta4, rho0=rho0, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, I=I, J=J, N=N)
    simdata$time <- factor(simdata$time)
    
    fit <- try(lmer(Y ~ -1 + time + W + time:X + W:X + (1|I_id) + (1|I_id:T_id) + (1|I_id:v_id), data=simdata), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    
    R <- c(rep(0,2*J+1),1)
    beta <- fixef(fit)
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 1-pchisq(test.stat, 1)
  }
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  
  return(c(empirical, error.rate))
}


#Function to compute empirical power or empirical type I error under cross-sectional design, ATE test 
empirical_CS_ATE <- function(nullcase=F, parameter, nsims=5000){
  # Argument:
  # nullcase: calculate empirical type I error or power. The former if nullcase=TRUE
  # parameter: the combination of design parameters of interest 
  # nsims: number of data replications
  
  I <- as.numeric(parameter[8])
  J <- as.numeric(parameter[2])
  N <- as.numeric(parameter[1])
  alpha0 <- as.numeric(parameter[3])
  alpha1 <- as.numeric(parameter[4])
  rho0 <- as.numeric(parameter[5])
  rho1 <- as.numeric(parameter[6])
  
  beta4 <- 0.2
  beta2 <- delta_OTE-mu_X*beta4
  
  if (nullcase==T){
    beta2 <- -mu_X*beta4
  }
  
  pvalue <- NULL
  count <- NULL
  for (i in 1:nsims){
    set.seed(0704+2022*i)
    simdata <- gendata_CS(beta2=beta2, beta4=beta4, rho0=rho0, rho1=rho1, alpha0=alpha0, alpha1=alpha1, I=I, J=J, N=N)
    simdata$time <- factor(simdata$time)
    
    fit <- try(lme(Y ~ -1 + time + W + time:X_centered + W:X_centered, data=simdata, random=list(I_id = ~ 1, T_id = ~ 1)), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    
    R <- c(rep(0,J), 1, rep(0,J+1))
    beta <- fit$coef$fixed
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 1-pf(test.stat, 1, I-2)
  }
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  return(c(empirical, error.rate))
}


#Function to compute empirical power or empirical type I error under closed-cohort design, ATE test 
empirical_CC_ATE <- function(nullcase=F, parameter, nsims=5000){
  # Argument:
  # nullcase: calculate empirical type I error or power. The former if nullcase=TRUE
  # parameter: the combination of design parameters of interest 
  # nsims: number of data replications
  
  I <- as.numeric(parameter[8])
  J <- as.numeric(parameter[2])
  N <- as.numeric(parameter[1])
  alpha0 <- as.numeric(parameter[3])
  alpha1 <- as.numeric(parameter[4])
  alpha2 <- as.numeric(parameter[5])
  rho0 <- as.numeric(parameter[6])
  
  beta4 <- 0.1
  beta2 <- delta_OTE-mu_X*beta4
  
  if (nullcase==T){
    beta2 <- -mu_X*beta4
  }
  
  pvalue <- NULL
  count <- NULL
  for (i in 1:nsims){
    set.seed(0704+2022*i)
    simdata <- gendata_CC(beta2=beta2, beta4=beta4, rho0=rho0, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, I=I, J=J, N=N)
    simdata$time <- factor(simdata$time)
    
    fit <- try(lmer(Y ~ -1 + time + W + time:X_centered + W:X_centered + (1|I_id) + (1|I_id:T_id) + (1|I_id:v_id), data=simdata), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    
    R <- c(rep(0,J),1,rep(0,J+1))
    beta <- fixef(fit)
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 1-pf(test.stat, 1, I-2)
  }
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  
  return(c(empirical, error.rate))
}


#Function to compute empirical power or empirical type I error under cross-sectional design, ATE test, unadjusted
empirical_CS_ATE_unadj <- function(nullcase=F, parameter, nsims=5000){
  # Argument:
  # nullcase: calculate empirical type I error or power. The former if nullcase=TRUE
  # parameter: the combination of design parameters of interest 
  # nsims: number of data replications
  
  I <- as.numeric(parameter[8])
  J <- as.numeric(parameter[2])
  N <- as.numeric(parameter[1])
  alpha0 <- as.numeric(parameter[3])
  alpha1 <- as.numeric(parameter[4])
  rho0 <- as.numeric(parameter[5])
  rho1 <- as.numeric(parameter[6])
  
  beta2 <- delta_OTE-mu_X*beta4
  
  if (nullcase==T){
    beta2 <- -mu_X*beta4
  }
  
  pvalue <- NULL
  count <- NULL
  for (i in 1:nsims){
    set.seed(0704+2022*i)
    simdata <- gendata_CS(beta2=beta2, beta4=beta4, rho0=rho0, rho1=rho1, alpha0=alpha0, alpha1=alpha1, I=I, J=J, N=N)
    simdata$time <- factor(simdata$time)
    
    fit <- try(lme(Y ~ -1 + time + W, data=simdata, random=list(I_id = ~ 1, T_id = ~ 1)), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    
    R <- c(rep(0,J), 1)
    beta <- fit$coef$fixed
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 1-pf(test.stat, 1, I-2)
  }
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  
  return(c(empirical, error.rate))
}


#Function to compute empirical power or empirical type I error under closed-cohort design, ATE test, unadjusted
empirical_CC_ATE_unadj <- function(nullcase=F, parameter, nsims=5000){
  # Argument:
  # nullcase: calculate empirical type I error or power. The former if nullcase=TRUE
  # parameter: the combination of design parameters of interest 
  # nsims: number of data replications
  
  I <- as.numeric(parameter[8])
  J <- as.numeric(parameter[2])
  N <- as.numeric(parameter[1])
  alpha0 <- as.numeric(parameter[3])
  alpha1 <- as.numeric(parameter[4])
  alpha2 <- as.numeric(parameter[5])
  rho0 <- as.numeric(parameter[6])
  
  beta2 <- delta_OTE-mu_X*beta4
  
  if (nullcase==T){
    beta2 <- -mu_X*beta4
  }
  
  pvalue <- NULL
  count <- NULL
  for (i in 1:nsims){
    set.seed(0704+2022*i)
    simdata <- gendata_CC(beta2=beta2, beta4=beta4, rho0=rho0, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, I=I, J=J, N=N)
    simdata$time <- factor(simdata$time)
    
    fit <- try(lmer(Y ~ -1 + time + W + (1|I_id) + (1|I_id:T_id) + (1|I_id:v_id), data=simdata), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1

    R <- c(rep(0,J),1)
    beta <- fixef(fit)
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 1-pf(test.stat, 1, I-2)
  }
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  
  return(c(empirical, error.rate))
}



