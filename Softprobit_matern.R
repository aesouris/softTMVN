# Written by Anirban Bhattacharya and Allyson Souris
# Used for "The Soft Multivariate Truncated Distribution" paper by Allyson Souris, Debdeep Pati, and Anirban Bhattacharya: https://arxiv.org/abs/1807.09155

# -- function to fit Gaussian Process Probit model with Matern Karnel -- #


# See section 4.1 in "The Soft Multivariate Truncated Distribution" 

# samples theta ~ N_{n}(0, K) 1_{C}(theta)
# where K is the Matern Kernel matrix
# C = C_1 x ... x C_n
# C_i = (0,infty) if y_i = 1 and (-infty, 0) if y_i = 0
# eta is the hyperparameter for the logit approximation

# motivated by a probit Gaussian Process model:
# Y_i = 1{Z(s_i) >0} 
# Z(s_i) ~ GP(0,K) where K is the matern kernel
# s_1, ..., s_n are observed locations

library(BayesLogit) #sample from polya gamma distribution: rpg)

Matern = function(x, y ,l, nu){
  #function provided by Dr. Bhattacharya
  #Calculates the matern kernel matrix
  ifelse(abs(x-y)>0, (sqrt(2*nu)*abs(x-y)/l)^nu/(2^(nu-1)*gamma(nu))*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu), 1.0)
} 

softprobit_matern = function(y,s,eta,nu = 3/5,BURNIN=1000,MCMC=10000,THIN=5){
  # y is nx1 binary vector, s is nx1 vector of observed locations
  # eta > 0 is the constant in the sigmoidal function approximating the indicator
  # nu is the smoothness parameter for the matern kernel
  
  # global constants
  ptm = proc.time()
  n = length(s)
  NITER = BURNIN + MCMC
  EFF = MCMC/THIN

  kappa = y-0.5
  
  #create the covariance from the matern kernel
  Sigma <- matrix(NA, nrow = n, ncol = n)
  for(i in 1:length(s)) {
    for(j in i:length(s)) {
      m = Matern(s[i], s[j], 1, nu)
      Sigma[i,j] = Sigma[j,i] = m
    }
  }
  
  #get Sigma^{-1}
  R_Sig = chol(Sigma) #Sigma = t(R_Sig) %*% R_Sig where R_Sig is upper triangular, Sigma^-1 = R_Sig^-1 %*% t(R_Sig)^-1
  Rinv = solve(R_Sig)
  SigInv = Rinv %*% t(Rinv)

  # initialize theta = z
  z = rep(0,n)

  # output files  
  zout = matrix(0,nrow=n,ncol=EFF) 

  # begin MCMC
  for(g in 1:NITER) {
    
    # -- sample omega_i for i = 1,...,n -- #
    omega = rpg(n, 1, eta*z)

    # -- sample z -- #
    #Use Rue 2001 sampling scheme

    x = rnorm(n, 0, 1)                # x~N(0,1)
    Q = SigInv + eta^2*diag(omega)    # z|- ~ N(Q^{-1}kappa*eta, Q^{-1})
    R_Q = chol(Q)                     # Q = t(R) %*% R
    v = solve(t(R_Q), eta*kappa)          # v = t(R)^{-1}kappa
    w = solve(R_Q, v)                 # w = R^{-1}v
    u = solve(R_Q, x)                 # u = R^{-1}x
    z = w + u                         # theta ~ N(Q^{-1}kappa, Q^{-1})

    
    if(g > BURNIN && g%%THIN==0) {
      zout[,(g-BURNIN)/THIN] = z
    }
    
  }
  tt = proc.time() - ptm
  postmean = apply(zout,1,mean)
  return(list('elp_time' = tt,'z_matrix' = zout,'postmean_z' = postmean))    
}

