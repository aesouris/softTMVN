# Written by Anirban Bhattacharya and Allyson Souris
# Used for "The Soft Multivariate Truncated Distribution" paper by Allyson Souris, Debdeep Pati, and Anirban Bhattacharya: https://arxiv.org/abs/1807.09155

# -- function to fit sample truncated multivariate normal distribution from probit-Gaussian motivation -- #

# See section 4.2 in "The Soft Multivariate Truncated Distribution" 

# samples theta ~ N_{n+p}(0, Sigma) 1_{C}(theta)
# where Sigma has the form [I_n + X Lambda X^T & X Lambda \\ Lambda X^T & Lambda]
# Lambda = diag{ xi_j^{-1} }_{j=1}^{p}
# C = C_1 x ... x C_n x R^p
# C_i = (0,infty) if y_i = 1 and (-infty, 0) if y_i = 0
# eta_s is the hyperparameter for the logit approximation: called eta in paper

# motivated by a probit regression model:
# y_i = 1{z_i > 0}, i = 1, ..., n
# z_i ~ N(X_i^T * beta, 1)
# beta_j ~ N(0, xi_j^{-1}), j = 1, ..., p
# then theta = [z, beta]^T

library(BayesLogit) #sample from polya gamma distribution: rpg

softprobit_ridge = function(y,X,eta_s,xi,BURNIN=1000,MCMC=10000,THIN=5){
  # y is nx1 binary vector, X is nxp design matrix
  # eta_s > 0 is the constant in the sigmoidal function approximating the indicator
  # xi is the prior precision for beta
  
  # global constants
  ptm = proc.time()
  n = nrow(X)
  p = ncol(X)
  NITER = BURNIN + MCMC
  EFF = MCMC/THIN
  I_n=diag(n)
  
  kappa = y-0.5
  G = X%*%((1/xi)*t(X))                       # G = X %*% Lambda %*% X'
  
  # initialize theta = [z,beta]
  theta = rep(0,n+p)
  z = theta[1:n]
  beta = theta[(n+1):(n+p)]
  
  # output files  
  betaout = matrix(0,nrow=p,ncol=EFF)  # only save the beta part 
  zout = matrix(0,nrow=n,ncol=EFF)  # only save the z part 
  thetaout = matrix(0,nrow=n+p,ncol=EFF)  # save the all of theta 
  
  # begin MCMC
  for(g in 1:NITER) {
    
    # -- sample omega_i for i = 1,...,n -- #
    omega = rpg(n, 1, eta_s*z)
    omegahalf = sqrt(omega)
    
    # -- sample theta -- #
    
    # calculate u
    delta0 = rnorm(n)
    betasamp = (1/sqrt(xi))*rnorm(p)             # betasamp ~ N(0,xi^{-1} I_p)
    zsamp = X%*%betasamp + delta0                # (zsamp,betasamp) ~ N(0,Sigma)
    #u = c(zsamp,betasamp)
    
    # calculate v
    delta = rnorm(n)
    v = eta_s*(omegahalf*zsamp) + delta
    
    # calculate w
    alpha = (1/omegahalf)*kappa
    tempmat = (I_n + G)*omegahalf        # tempmat = Omega^(1/2)*(I + (1/xi) XX')
    mat = (eta_s^2)*t(t(tempmat)*omegahalf)         # mat = Phi Sigma Phi'
    w = solve((I_n+mat),(alpha-v))
    
    # set theta by calculating z and beta
    z = zsamp + eta_s*t(tempmat)%*%w
    beta = betasamp + (eta_s/xi)*t(X*omegahalf)%*%w    # u[n+1:p] + eta_s*(1/xi)*X' Omega^(1/2) w
    
    
    
     
    if(g > BURNIN && g%%THIN==0) {
      betaout[,(g-BURNIN)/THIN] = beta
      zout[,(g-BURNIN)/THIN] = z
      thetaout[,(g-BURNIN)/THIN] = c(z, beta)
    }
    
  }
  tt = proc.time() - ptm
  postmean = apply(betaout,1,mean)
  return(list('elp_time' = tt,'beta_matrix' = betaout,'postmean_beta' = postmean, 'z_matrix' = zout, 'theta_matrix' = thetaout))    
}

