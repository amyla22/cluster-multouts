#############################################################################
###########  multOutData() is a function to simulate multiple outcomes data
###########    within a full simulation for our multiple outcomes model with
###########    clustering of outcomes based on DP prior. The data simulated
###########    here closely follows the simulated data method in Xiao et al.
#############################################################################
library(MASS)


multOutData <- function(betaf, bD, betaDf, s2bO, s2rD, s2r, n, letter, domAssn, seed, seedfile, Dstart){
  ### Function Argurments:
  ## Parameter Value specifications: betaf, bD, betaDf, s2rD, s2bO, s2r
  # J = no. of outcomes
  # n = no. of subjects
  # D = no. of prespecified domains
  # pF = no. of overall fixed effects covariates
  # pDF = no. of covariates w/ domain specific fixed effects
  # betaf, betaDf, s2rD, s2r = parameter value specifications
  # domAssn = true assignment of outcomes into domains
  # seed = seed for simulating data set
  # nmis = number of total missing observed outcomes --> these will be randomly sampled
  # Dstart = number of domains to use as a starting point (if DPP=FALSE this must be equal to K)
  ## Could optionally add the following if considering multiple exposures:
  # pF = no. of overall fixed effects covariates
  # pD = no. of covariates w/ domain specific random effects
  # pO =  no. of covariates w/ outcome specific random effects
  simseed <- seedfile[seed,]
  cat('random seed: ',simseed , '\n')
  set.seed(simseed)
  
  # ********************************************************
  # *********     DOMAIN ASSIGNMENTS          *********
  # ********************************************************
  J = length(domAssn) # number of outcomes
  D = length(unique(domAssn)) # number of domains
  ## Membership matrix (JxD) with j,d^th entry = 1 if outcome j is in d 
  #   Each row will sum to 1 since each outcome can be in only 1 domain
  if(D>1){
    Lambda = cbind(1*(domAssn==1),model.matrix( ~ as.factor(domAssn)  ,contrasts='contr.sum')[,-1]) 
  }else{
    Lambda = cbind(rep(1,J))
  }
  nsubjects = n
  
  # ********************************************************
  # *********       OVERALL FIXED EFFECTS          *********
  # ********************************************************
  # Make these values be somewhat realistic, but stronger
  # Note: this part should correspond to the structure of F, S_D and S_O.
  # slope for exposure (\beta_x)
  beta = betaf
  pF = length(betaf)
  
  # ********************************************************
  # *********    DOMAIN-SPECIFIC FIXED EFFECTS     *********
  # ********************************************************
  ncovs = nrow(betaDf)
  #cat("Domain-specific fixed effects for  covariates:\n")
  betaDF = betaDf
  #print(betaDF)
  
  # ********************************************************
  # *********   DOMAIN-SPECIFIC RANDOM EFFECTS     *********
  # ********************************************************
  #cat("Domain-specific random effects for hg and covariates:\n")
  b.D = bD
  pD = length(bD)/D
  #print(b.D)
  
  # ********************************************************
  # *********  OUTCOME-SPECIFIC RANDOM EFFECTS      ********
  # ********************************************************
  b.O = rnorm(J, mean = 0, sd = sqrt(s2bO))
  pO = length(s2bO)
  #cat("Outcome-specific random effects for hg:\n")
  #print(b.O)
  var(b.O)
  
  Beta <- beta + Lambda%*%b.D + b.O   
  
  BetaDF <- Lambda%*%t(betaDF)
  
  ## Variance components for r AND r.D
  sig.r = sqrt(s2r)   # sqrt(var) of overall subject-specific rand effect
  #cat("sig.r1=",sig.r,"\n")
  
  sig.rDvec = sqrt(s2rD)
  #cat("sig.r2vec=",sig.rDvec,"\n")
  
  #### Construct the data covariance structure
  ## ---> QUESTION: Where is the domain-specific Hg RE???
  if(D>1){
    Sigma <- matrix(0, J+1+ncovs+1+D, J+1+ncovs+1+D)
    diag(Sigma)[(J+1+ncovs+1)+(1:D)] <- sig.rDvec^2 
  }else{
    Sigma <- matrix(0, J+1+ncovs+1, J+1+ncovs+1)
  }
  diag(Sigma)[1:(J+1+ncovs)] <- 1  # outcome-specific RE, overall Hg effect, and domain-specific RE (scaled data)
  diag(Sigma)[J+1+ncovs+1] <- sig.r^2
  
  ### Calculate pairwise correlations between outcomes:
  for(j in 1:(J-1))
    for(k in (j+1):J)
      Sigma[j,k] <- Beta[j]*Beta[k] + t(BetaDF[j,])%*%BetaDF[k,] + sig.r^2 + sig.rDvec[Lambda[j,]==1]^2*t(Lambda[j,])%*%Lambda[k,]
  
  ### The pairwise correlations between each outcome and the exposure is the sum of the coefficients (betaf + bD):
  for(j in 1:J)
    Sigma[j,J+1] <- Beta[j]
  
  
  ### The pairwise correlations between each outcome and each covariate, k, is the coefficient (betaDf_k):
  for(j in 1:J)
    for(k in 1:ncovs)
      Sigma[j,J+1+k] <-  BetaDF[j,k]
  
  
  ### The pairwise correlations between each outcome and the subject-specific random effect is sigr2:
  for(j in 1:J) 
    Sigma[j,J+1+ncovs+1] <- sig.r^2
  
  ### The pairwise correlations between each outcome and the domain-specific subject random effect is sigr2D:
  if(D>1){
    for(j in 1:J)
      for(d in 1:D)
        Sigma[j,J+1+ncovs+1+d] <- Lambda[j,d]*sig.rDvec[d]^2
  }
  
  ### fill in lower triangular part of matrix because symmetric:
  if(D>1){
    for(j in 1:(J+1+ncovs+D))
      for(k in (j+1):(J+1+ncovs+D+1))
        Sigma[k,j] <- Sigma[j,k]
  }else{
    for(j in 1:(J+1+ncovs))
      for(k in (j+1):(J+1+ncovs+1))
        Sigma[k,j] <- Sigma[j,k]
  }
  rawdat = mvrnorm(n=nsubjects, mu = rep(0,nrow(Sigma)), Sigma=Sigma)
  
  Y = matrix(rawdat[,1:J],ncol=1)
  Sf = SD = SO = matrix(rawdat[,J+1],ncol=1)
  SDf = matrix(rawdat[,(J+1+1):(J+1+ncovs)],nrow = n, ncol = ncovs)
  y.jlabs = rep(1:J,each=nsubjects)
  y.ilabs = rep(1:nsubjects,J)
  
  ### Priors set corresponding to Thurston et al., 
  prior <- list(beta0f=rep(0,pF), SigInv.betaf=diag(.00001,pF), 
                      beta0Df=rep(0,pDF), SigInv.betaDf=diag(.00001,pDF), 
                      sig2.bO=c(0.5,0.00005), sig2.bD=c(0.5,0.00005), sig2.eps=c(0.5,0.0005),
                      sig2.r=c(0.5,0.00005), sig2.rD=c(0.5,0.00005), alpha=1)
  
  inits <- list(D.init=Dstart, dvec=rep(1:Dstart,length.out=J), betaF = rep(0,pF), bO = rep(0,pO*J), 
                r = rep(0,n), sig2bO=rep(.05,pO), sig2bD=rep(.05,pD), sig2eps=rep(.05,J),  sig2r=.05,
                betaDf = lapply(1:Dstart, function(i) as.vector(matrix(0,nrow=Dstart,ncol=pDF)[i,,drop=FALSE])), 
                bD = lapply(1:Dstart, function(i) as.vector(matrix(0,nrow=Dstart,ncol=pD)[i,,drop=FALSE])), 
                rD=lapply(1:Dstart, function(i) as.vector(matrix(0,nrow=Dstart,ncol=n)[i,,drop=FALSE])),
                sig2rD=lapply(1:Dstart, function(i) as.vector(matrix(0.05,nrow=Dstart,ncol=1)[i,,drop=FALSE])))
  
  
  fileName <- paste('J',J,'n',n,'D',D,letter,'seed',seed,'txt',sep='.')
  return(list(params=list(betaf=betaf, bD=bD, betaDf=betaDf, s2rD=s2rD, s2r=s2r, n=n, J=J, D=D, domAssn=domAssn, seed=seed),
              covs=list(Sf=Sf, SO=SO, SD=SD, SDf=SDf), 
              Y=Y, y.jlabs=y.jlabs, y.ilabs=y.ilabs, priors=prior, inits=inits, fileName=fileName))
  
}

