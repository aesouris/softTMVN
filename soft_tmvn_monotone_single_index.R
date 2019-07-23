## Written by Anirban Bhattacharya and Allyson Larsen
## Used for "The Soft Multivariate Truncated Distribution" paper by Allyson Larsen, Debdeep Pati, and Anirban Bhattacharya: https://arxiv.org/abs/1807.09155

## Samples from the posterior of the monotone single index model
## See section 5 in "The Soft Multivariate Truncated Distribution" 

## monotone single index model:
## y_i = f(x_i^T %*% alpha) + epsilon_i, i = 1,...,n, x_i, alpha in R^p, y_i in R
## epsilon_i ~ N(0, sigma^2)
## x_i are centered, standardized and x_i = x_i/C where c = max_i{||x_i||}
## f(x_i^T %*% alpha) = B_M^i %*% theta
## B_M^i = [B_{M,0}(x_i^T %*% alpha), ..., B_{M,M}(x_i^T %*% alpha)]
## B_{M,j}(t) = M choose j * 1/2* ((t+1)/2)^j *(1-(t+1)/2)^{M-j} - the bernstein polynomials with support [-1,1]
## B_{M,j}(t) = 1/(M+1)*1/2*p_j((t+1)/2), p_j(u) is the beta(j+1, M-j+1) density function
## psi_0 = theta_0, psi_1 = theta_1-theta_0, ..., phi_M = theta_M - theta_{M-1}
## psi_0 ~ N(0, 25)
## psi_k ~ N(0, 25)1{psi_k \geq 0}, k = 1, ..., M
## alpha = beta/||beta||
## beta ~ N(0,0.1^2)
## sigma^2 ~ IB(2.1,1.1) E[] = 1, Var() = 10

## msim_tmvn_gibbs samples assuming psi has tMVN distribution, msim_soft_gibbs samples assuming psi has soft tMVN distribution

##Updated 7/1/19

library(BayesLogit) #sample from polya gamma distribution: rpg
library(TruncatedNormal) #botev's tMVN sampler
library(tmvtnorm)

# calculates B_{M,j}(t)
Btilde <- function(M, j, t) {
  #M - degree of bernstein polynomials, number of basis functions: M+1
  #j - indicates which of the M+1 bernstein polynomials are being evaluated 
  #t - the value the polynomial is being evaluated at, in [-1,1]
  fun <- 1/(2*(M+1))*dbeta((t+1)/2, j+1, M-j+1)
  return(fun)
}

# caluculates f(t) for a vector of t values
BigB <- function(M, t, theta) {
  #M - degree of bernstein polynomials, number of basis functions: M+1
  #t - the value the polynomial is being evaluated at, in [-1,1]
  #theta - the vector of theta
  
  values <- sapply(0:M, Btilde, M = M, t = t) #rows are over t, columns are over 0:M, values[i,j] is Btilde(M,j,t[i])
  fun <- values %*% theta
  return(as.vector(fun))
}

#calculates D_alpha matrix
D_mat_fun <- function(M, t) {
  #t is a vector of length n, t in [-1,1]
  #M is the degree of the bernstein polynomials
  
  values <- sapply(M:0, Btilde, M = M, t = t) #rows are over t, columns are over M:0, values[i,j] is Btilde(M,(M-j),t[i])
  D_mat <- apply(values, 1, cumsum)
  
  return(t(D_mat[(M+1):1,]))
}

#gibbs sampler assuming psi has soft tMVN prior
msim_soft_gibbs <- function(X, Y, M, acc_sig, eta = 100, BURNIN=1000, MCMC=5000, THIN=1) {
  #input: data(x_i,y_i) i = 1,..,n, x_i in R^p, y_i in R
  #input: M, the degree of the Bernstein polynomials
  #input: acc_sig, the sigma value that controls the proposed value of beta in the Metropolis step
  #input: eta = 100, the hyperparameter of the soft tMVN distribution
  #input: BURNIN, MCMC, THIN, the amount of BURNIN, after burnin MCMC samples, and thinning parameter for the Gibbs sampler

  #global constants:
  n <- nrow(X)
  p <- ncol(X)
  EFF <- MCMC/THIN
  NITER <- MCMC + BURNIN
  Im1 <- diag(M+1) #The (M+1)x(M+1) identity matrix
  
  #set hyperparameter values
  a <- 2.1
  b <- 1.1
  kappa <- c(0,rep(1/2,M))
  
  #initialize values:
  beta <- rnorm(p,0,1)
  alpha <- beta/sqrt(sum(beta^2))
  sig_sq <- 1 #note: using sigma^2 (var) not sigma (sd)
  psi <- rep(1, M+1)
  
  #Get D_alpha
  t <- X%*%alpha
  D_alpha <- D_mat_fun(M, t)
  
  #caluculate ||y - D_alpha %*% psi||_2^2
  ydpsi <- Y - D_alpha %*% psi
  norm_val <- sum(ydpsi^2)
  
  #output vectors
  alpha_mat <- matrix(NA, nrow = p, ncol = EFF)
  sig_sq_vec <- rep(NA, EFF)
  psi_mat <- matrix(NA, nrow = M+1, ncol = EFF)
  accept_prob_vec <- rep(NA,EFF)
  change_vec <- rep(NA, NITER)
  
  #gibbs sampler:
  for(gloop in 1:NITER) {
    
    ## sample from sigma^2|-
    sig_sq <- 1/rgamma(1, shape = a + n/2, rate = b + 1/2*norm_val)

    ## sample from psi|-
    #caclulate D_alpha^T %*% D_alpha
    DTD <- t(D_alpha)%*%D_alpha
    
    # sample w|psi,-
    omega <- rpg(M, 1, eta*psi[2:(M+1)])

    #sample psi|omega,-
    #use Rue (2001) sampling scheme: N(Q^{-1}b, Q^{-1})
    Q_mat <- DTD/sig_sq + 1/25 * Im1 + eta^2*diag(c(0, omega))
    b_vec <- 1/sig_sq*t(D_alpha)%*%Y + eta*kappa
    
    R <- chol(Q_mat) #Q = R^T %*% R = L %*% L^T
    z <- rnorm(M+1, 0, 1) #z ~ N(0, I)
    v <- solve(t(R), b_vec) #v = (R^T)^{-1}b
    mu <- solve(R, v) #mu = R^{-1} v = R^{-1} (R^T)^{-1}b = Q^{-1}b
    nm <- solve(R, z) #nm = R^{-1}z ~ n(0, q^{-1})
    psi <- mu+nm #~N(Q^{-1}b, Q^{-1})
    
    ## sample from alpha|-
    #get proposed values
    beta_star <- rnorm(p, beta, acc_sig)
    alpha_star <- beta_star/sqrt(sum(beta_star^2))
    
    #calculate D_alpha_star
    t_star <- X%*%alpha_star
    D_alpha_star <- D_mat_fun(M, t_star)
    
    #caluculate ||y - D_alpha^* %*% psi||_2^2
    ydpsi_star <- Y - D_alpha_star %*% psi
    norm_val_star <- sum(ydpsi_star^2)
    
    #calculate the acceptance prob
    accept_prob <- min(exp(-1/(2*sig_sq) * norm_val_star - 1/2*sum(beta_star^2) + 1/(2*sig_sq) * norm_val + 1/2*sum(beta^2)),1)
    
    #Metropolis step
    u <- runif(1,0,1)
    if(u <= accept_prob) {
      beta <- beta_star
      alpha <- alpha_star
      ydpsi <- ydpsi_star
      norm_val <- norm_val_star
      t <- t_star
      D_alpha <- D_alpha_star
      change_vec[gloop] <- 1
    }
    
    #store output
    if(gloop > BURNIN && gloop%%THIN==0) {
      alpha_mat[,(gloop-BURNIN)/THIN] <- alpha
      sig_sq_vec[(gloop-BURNIN)/THIN] <- sig_sq
      psi_mat[,(gloop-BURNIN)/THIN] <- psi
      accept_prob_vec[(gloop-BURNIN)/THIN] <- accept_prob
      
    }
  }
  
  return(list('alpha' = alpha_mat, 'sigma_sq' = sig_sq_vec, 'psi' = psi_mat, 'acceptance' = accept_prob_vec, 'change' = change_vec))
  
}

#Gibbs sampler assuming psi has tMVN prior from Botev
msim_tmvn1_gibbs <- function(X, Y, M, acc_sig, BURNIN=1000, MCMC=5000, THIN=1) {
  #input: data(x_i,y_i) i = 1,..,n, x_i in R^p, y_i in R
  #input: M, the degree of the Bernstein polynomials
  #input: acc_sig, the sigma value that controls the proposed value of beta in the Metropolis step
  #input: BURNIN, MCMC, THIN, the amount of BURNIN, after burnin MCMC samples, and thinning parameter for the Gibbs sampler
  #output:
  
  #global constants:
  n <- nrow(X)
  p <- ncol(X)
  EFF <- MCMC/THIN
  NITER <- MCMC + BURNIN
  Im1 <- diag(M+1) #The (M+1)x(M+1) identity matrix
  
  #set hyperparameter values
  a <- 2.1
  b <- 1.1
  kappa <- c(0,rep(1/2,M))
  lower <- c(-Inf, rep(0, M)) #lower bound for tMVN
  upper <- rep(Inf, M+1) #upper bound for tMVN
  
  #initialize values:
  beta <- rnorm(p,0,1)
  alpha <- beta/sqrt(sum(beta^2))
  sig_sq <- 1 #note: using sigma^2 (var) not sigma (sd)
  psi <- rep(1, M+1)
  
  #Get D_alpha
  t <- X%*%alpha
  D_alpha <- D_mat_fun(M, t)
  
  #caluculate ||y - D_alpha^* %*% psi||_2^2
  ydpsi <- Y - D_alpha %*% psi
  norm_val <- sum(ydpsi^2)
  
  #output vectors
  alpha_mat <- matrix(NA, nrow = p, ncol = EFF)
  sig_sq_vec <- rep(NA, EFF)
  psi_mat <- matrix(NA, nrow = M+1, ncol = EFF)
  accept_prob_vec <- rep(NA,EFF)
  change_vec <- rep(NA, NITER)
  
  #gibbs sampler:
  for(gloop in 1:NITER) {
    
    ## sample from sigma^2|-
    sig_sq <- 1/rgamma(1, shape = a + n/2, rate = b + 1/2*norm_val)
    
    ## sample from psi|-
    DTD <- t(D_alpha)%*%D_alpha
    Sigma <- solve(DTD/sig_sq + 1/25 * Im1)
    mu <- 1/sig_sq*Sigma %*%t(D_alpha)%*%Y
    
    nm <- mvrandn(lower-mu,upper-mu,Sigma,1) #tMVN(0,Sigma)
    psi <- mu+nm #tMVN(mu, Sigma)
    
    ## sample from alpha|-
    #get proposed values
    beta_star <- rnorm(p, beta, acc_sig)
    alpha_star <- beta_star/sqrt(sum(beta_star^2))
    
    #calculate D_alpha_star
    t_star <- X%*%alpha_star
    D_alpha_star <- D_mat_fun(M, t_star)
    
    #caluculate ||y - D_alpha^* %*% psi||_2^2
    ydpsi_star <- Y - D_alpha_star %*% psi
    norm_val_star <- sum(ydpsi_star^2)

    #calculate the acceptance prob
    accept_prob <- min(exp(-1/(2*sig_sq) * norm_val_star - 1/2*sum(beta_star^2) + 1/(2*sig_sq) * norm_val + 1/2*sum(beta^2)),1)
    
    #Metropolis step
    u <- runif(1,0,1)
    if(u <= accept_prob) {
      beta <- beta_star
      alpha <- alpha_star
      ydpsi <- ydpsi_star
      norm_val <- norm_val_star
      t <- t_star
      D_alpha <- D_alpha_star
      change_vec[gloop] <- 1
    }
    
    #store output
    if(gloop > BURNIN && gloop%%THIN==0) {
      alpha_mat[,(gloop-BURNIN)/THIN] <- alpha
      sig_sq_vec[(gloop-BURNIN)/THIN] <- sig_sq
      psi_mat[,(gloop-BURNIN)/THIN] <- psi
      accept_prob_vec[(gloop-BURNIN)/THIN] <- accept_prob
      
    }
  }
  
  return(list('alpha' = alpha_mat, 'sigma_sq' = sig_sq_vec, 'psi' = psi_mat, 'acceptance' = accept_prob_vec, 'change' = change_vec))
  
}
