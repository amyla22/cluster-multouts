###########################################################################################################
###########################################################################################################
##### Document of FINAL functions for clustering multiple outcomes into domains via DPP
###########################################################################################################

library(MASS)
library(MCMCpack)
library(mvtnorm)
library(cluster)
library(mcmcse)

###########################################################################################################
##### clust.multouts() is a function which takes the data, priors and starting values for an MCMC designed 
#####         to simultaneously fit the multiple outcomes model and group the outcomes into domains. The
#####         grouping of outcomes into domains is executed using the Dirichlet process prior; hence, the
#####         number of groups is not prespecified.
###########################################################################################################
#####         The function takes the (1) the data vector, Y, stacked by outcomes, 
#####         (2) covariate matrices for each of the ones specified in Luo et al. 2014: Sf, SDf, SD, SO  
#####         (3/4) the outcome labels y.jlabs, and the observation labels, y.ilabs, for the Y vector 
#####         (5/6)  priors as a list, priors, and the initial values as a list, inits
#####         (7/8) min/max number of iterations to perform the mcmc sampler (9) prints iteration num if T
#####         (10) the FIXED domain assignments (11) the min ESS for bO.draws (all J of these)
#####         (12) a T/F to indicate whether to introduce a hierarchical prior on sig2rD
###########################################################################################################
###########################################################################################################

cluster.multouts <- function(Y, Sf, SDf, SD, SO, y.jlabs, y.ilabs, priors, inits, m, tune, nburn, min.niter, min.niter2, max.niter, 
                             domAssn, printIter, printSumms, sim.seed ){
  #### y.jlabs = distinguish outcomes
  #### y.ilabs = distinguish subjects
  #### printIter = TRUE to show the iteration for the overall MCMC
  #### printSumms = TRUE to show the clustering process for each outcome at each iteration
  #### sim.seed = random seed used in simulating everything but the covariate matrices in the function multOutData() 
  #### m = # of proposed new states in draw.domA8()
  #### tune = # of Gibbs updates for launch split and launch merge in draw.dom.JN() 
  ## accepted: 1-accepted a split, 2-rejected a split, 3-accepted a merge, 4-rejected a merge
  
  # Define matrices/lists/vectors to store parameters which are the size of the maximum number of iterations
  niter <- max.niter
  ##########################################
  ########## Define the following:
  nobs.vec <- table(y.jlabs)  # number of observations for each outcome (same if all have complete obs on all outcomes)
  pF <- ncol(Sf) 
  pDF<- ncol(SDf) 
  pD <- ncol(SD) 
  pO <- ncol(SO)
  J <- nlevels(factor(y.jlabs)) # number of outcomes
  n <- nlevels(factor(y.ilabs))# number of subjects/outcome
  y.mislabs <- apply(matrix(Y,ncol=1), 1, is.na) 
  nmis <- sum(y.mislabs)
  
  
  ##########################################
  ########## Matrices to Store Parameter Draws: 
  ## Non-domain specific parameters:
  betaf.draws <- matrix(NA, nrow=niter, ncol=pF)
  bO.draws <- matrix(NA, nrow=niter, ncol=pO*J)
  sig2bO.draws <- matrix(NA, nrow=niter, ncol=pO)
  r.draws <- matrix(NA, nrow=niter, ncol=n)
  sig2r.draws <- rep(NA, niter)
  sig2eps.draws <- matrix(NA, nrow=niter, ncol=J)
  sig2bD.draws <- matrix(NA, nrow=niter, ncol=pD)
  
  ## Domain-specific parameters: 
  betaDf.draws <- list()
  bD.draws <- list()
  rD.draws <- list()
  sig2rD.draws <- list()
  
  
  ## Domain assignments:
  d.draws <- matrix(NA,nrow=niter,ncol=J)
  accepted <- rep(NA,niter)
  
  ##########################################
  ########## Initial Values (entered from inits=list()) :
  betaf.draws[1,] <- inits$betaF
  bO.draws[1,] <- inits$bO
  colnames(bO.draws) <- rep(1:J,each=pO)
  sig2bO.draws[1,] <- inits$sig2bO
  r.draws[1,] <-  inits$r
  sig2r.draws[1] <- inits$sig2r
  sig2eps.draws[1,] <- inits$sig2eps
  sig2bD.draws[1,] <-inits$sig2bD
  
  ## need lists here.
  betaDf.draws[[1]] <- inits$betaDf
  bD.draws[[1]] <- inits$bD
  rD.draws[[1]] <- inits$rD
  sig2rD.draws[[1]] <- inits$sig2rD
  d.draws[1,] <- inits$dvec
  K <- nlevels(factor(d.draws[1,]))
  
  ##########################################
  ########## Prior Values (entered from inits=list()) :
  pr.beta0f <- priors$beta0f
  pr.Sigma0Invf <- priors$SigInv.betaf
  pr.beta0Df <- priors$beta0Df
  pr.Sigma0InvDf <- priors$SigInv.betaDf 
  pr.A0.bO <- priors$sig2.bO[1]
  pr.B0.bO <- priors$sig2.bO[2]
  pr.A0.bD <- priors$sig2.bD[1]
  pr.B0.bD <- priors$sig2.bD[2] 
  pr.A0.eps <- priors$sig2.eps[1]
  pr.B0.eps <- priors$sig2.eps[2]
  pr.A0.r <- priors$sig2.r[1]
  pr.B0.r <- priors$sig2.r[2]
  pr.A0.rD <- priors$sig2.rD[1]
  pr.B0.rD <- priors$sig2.rD[2] 
  pr.alpha <- priors$alpha
  
  MeHg.coef <- matrix(NA, nrow=niter, J)
  MeHg.coef[1,] <- betaf.draws[1,] + bO.draws[1,] + unlist(bD.draws[[1]][d.draws[1,]])
  MeHg.ess <- 2000 #minESS(J)
  
  ##########################################
  ##########  Membership Matrix - j rows represent the outcomes, d columns represent the membership of the
  #  jth outcome in the dth domain
  #   ---> this matrix will be updated at each iteration and may change in column length
  if(nlevels(factor(d.draws[1,]))>1){
    Lambda <- cbind(1*(d.draws[1,]==1),model.matrix( ~ as.factor(d.draws[1,])  ,contrasts='contr.sum')[,-1])
  }else{
    ## If all outcomes grouped into the same domain, need to have just a vector of 1s
    Lambda <- d.draws[1,]
  }
  RDvec <- kronecker(Lambda,diag(1,n))%*%unlist(rD.draws[[1]])
  
  ##########################################
  ##########  Covariate matrices which will NOT change with changing domains:   
  Fmat.F <- kronecker(rep(1, J),Sf)  
  Zmat.O <- kronecker(diag(1,J), SO) 
  
  t = 2
  while(t%in%2:niter){
    # cat('\n t=',t,' K=',K,'; d.draws=',d.draws[t-1,],'\n')
    if(printIter==TRUE){
      if(t%%1000==0) cat('\n t=',t,' K=',K,'; d.draws=',d.draws[t-1,],'\n')
      if(t>3){
        if(is.same(d.draws[t-2,],d.draws[t-1,])!=J) cat('\n t=',t,' K=',K,'; d.draws=',d.draws[t-1,],'\n')
      }
    } 
    ## Determine whether we have singletons in our data (also will not need several domain-specific parameters):
    Ktemp <- K - sum(table(d.draws[t-1,])==1)
    if(Ktemp<K){
      ## identify numeric value of the singleton domain(s)
      singleton <- as.numeric(which(table(d.draws[t-1,])==1))
    }else{
      singleton <- NULL
    } 
    #cat('\n t=',t,' K=',K,'; d.draws=',d.draws[t-1,],'\n')
    # ---> Redefine covariate matrices at each step based on outcome membership:
    Fmat.DF <- kronecker(Lambda,SDf) 
    Zmat.D <- kronecker(Lambda, SD)
    Rvec <- t(t(kronecker(c(rep(1,J)),r.draws[t-1,])))
    
    ## label outcomes according to domain assignment (may not need to keep track of this):
    d.labs <- rep(d.draws[t-1,],nobs.vec)  #rep(d.draws[t,],nobs.vec)
    ## number of outcomes in each domain:
    njD <- table(d.draws[t-1,])  #table(d.draws[t,]) 
    
    ########################
    ### ---> Draw missing values
    means.Y <- Fmat.F%*%(betaf.draws[t-1,]) + Fmat.DF%*%unlist(betaDf.draws[[t-1]]) + Zmat.D%*%unlist(bD.draws[[t-1]]) + Zmat.O%*%(bO.draws[t-1,]) + RDvec + Rvec
    sds.Y <- rep(sqrt(sig2eps.draws[t-1,]),each=n)
    Ymis <- rnorm(nmis, mean=means.Y[y.mislabs], sd=sds.Y[y.mislabs])
    Y[y.mislabs] <- Ymis    
    
    ###########################################################################################################
    #### ---> Gibbs Steps for Domain-Specific Parameters:
    ## EDIT: added in D_epsilon,d for means.Df..
    means.Df <- lapply(seq_along(1:K), FUN=function(k) pr.Sigma0InvDf%*%pr.beta0Df + 
                         t(kronecker(rep(1,njD[k]),SDf))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))
                       %*%(as.matrix(Y[d.labs==k],ncol=1) - kronecker(rep(1,njD[k]),(Sf%*%betaf.draws[(t-1),])) - kronecker(rep(1,njD[k]),(SD%*%(bD.draws[[t-1]][[k]])))
                           - kronecker(diag(1,njD[k]),SO)%*%bO.draws[(t-1),(colnames(bO.draws) %in% y.jlabs[d.labs==k])]
                           - as.matrix(kronecker(rep(1,njD[k]),r.draws[(t-1),]) - kronecker(rep(1,njD[k]),(rD.draws[[t-1]][[k]])),ncol=1)) )
    var.betaDf <- lapply(seq_along(1:K), FUN=function(k) solve(pr.Sigma0InvDf + sum(1/1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])])*t(SDf)%*%SDf) )#t(kronecker(rep(1,njD[k]),SDf))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%kronecker(rep(1,njD[k]),SDf)) )
    betaDf.draws[[t]] <- lapply(X=seq_along(1:K), FUN = function(k) mvrnorm(1, mu=unlist(var.betaDf[[k]])%*%unlist(means.Df[[k]]) , Sigma =unlist(var.betaDf[[k]])))  
    
    if(K>1){
      means.bD <- lapply(seq_along(1:K), FUN=function(k) t(kronecker(rep(1,njD[k]),SD))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%
                           (as.matrix(Y[d.labs==k],ncol=1) - kronecker(rep(1,njD[k]),(Sf%*%betaf.draws[(t-1),])) - kronecker(rep(1,njD[k]),(SDf%*%betaDf.draws[[t]][[k]])) 
                            - kronecker(diag(1,njD[k]),SO)%*%bO.draws[(t-1),(colnames(bO.draws) %in% y.jlabs[d.labs==k])] - as.matrix(kronecker(rep(1,njD[k]),r.draws[(t-1),]) 
                                                                                                                                      - kronecker(rep(1,njD[k]),(rD.draws[[t-1]][[k]])),ncol=1)) )
      var.bD <- lapply(seq_along(1:K), FUN=function(k) solve(diag(1/sig2bD.draws[(t-1),],pD) + t(kronecker(rep(1,njD[k]),SD))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%kronecker(rep(1,njD[k]),SD)) )
      bD.draws[[t]] <- lapply(X=seq_along(1:K), FUN = function(k) mvrnorm(1, mu=unlist(var.bD[[k]])%*%unlist(means.bD[[k]]) , Sigma =unlist(var.bD[[k]])))  
      if(!is.null(singleton)) bD.draws[[t]][singleton] <- 0
      
      bDd.sums <- apply(matrix(unlist(bD.draws[[t]]), nrow=pD),1,function(x) sum(x^2))
      if(!is.null(singleton)){
        sig2bD.draws[t,] <- rinvgamma(pD, shape=pr.A0.bD + (K-length(singleton))/2, scale=pr.B0.bD + .5*bDd.sums)
      }else{
        sig2bD.draws[t,] <- rinvgamma(pD, shape=pr.A0.bD + K/2, scale=pr.B0.bD + .5*bDd.sums)
      }
      
      var.rD <- lapply(seq_along(1:K), FUN=function(k) solve((1/sig2rD.draws[[t-1]][[k]]) + sum(1/sig2eps.draws[t-1,d.draws[t-1,]==k])) )
      separ.rD <- by((Y - Fmat.F%*%(betaf.draws[t-1,]) - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]])  - Zmat.O%*%(bO.draws[t-1,]) - Rvec),
                     INDICES=factor(y.jlabs),FUN=function(x) 1*x)
      means.rDbyout <- lapply(seq(1:J), FUN=function(j) (1/sig2eps.draws[t-1,j])*unlist(separ.rD[[j]]) )
      means.rDmat <- matrix(unlist(means.rDbyout),ncol=J,byrow=F)
      means.rD <- lapply(seq_along(1:K), FUN=function(k) 
        if( sum(d.draws[t-1,]==k)==1){
          means.rDmat[,k] 
        }else{
          rowSums(means.rDmat[,d.draws[t-1,]==k])
        }        )
      rD.draws[[t]] <- lapply(seq_along(1:K), FUN=function(k) rnorm(n, mean=unlist(var.rD[[k]])%*%unlist(means.rD[[k]]) , sd=sqrt(unlist(var.rD[[k]])))) 
      if(!is.null(singleton)) for(i in 1:length(singleton)) rD.draws[[t]][[singleton[i]]] <- rep(0,n)
      
      sig2rD.draws[[t]] <- lapply(seq(1:K), FUN=function(k)  rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(unlist(rD.draws[[t]][[k]]))%*%unlist(rD.draws[[t]][[k]])))
    }else{
      bD.draws[[t]] <- 0
      sig2bD.draws[t,] <- rinvgamma(pD, shape=pr.A0.bD + K/2, scale=pr.B0.bD )
      rD.draws[[t]] <- list(rep(0,n))
      sig2rD.draws[[t]] <- rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD )
    }
    RDvec <- kronecker(Lambda,diag(1,n))%*%unlist(rD.draws[[t]])
    ###########################################################################################################
    
    ###########################################################################################################
    #### ---> Gibbs Steps for non-domain specific parameters:
    ## Changed to sum over the n-subjects
    VbetafInv <- solve(pr.Sigma0Invf + sum(1/sig2eps.draws[t-1,])*t(Sf)%*%Sf)
    means.betaf1 <- VbetafInv%*%(pr.Sigma0Invf%*%pr.beta0f +(rep(1/sig2eps.draws[t-1,],each=n)*rep(Sf,J))%*%(Y - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]]) - Zmat.O%*%(bO.draws[(t-1),]) - Rvec - RDvec))
    betaf.draws[t,] <- rnorm(1, mean=means.betaf1,  sd=sqrt(VbetafInv)) 
    
    VbOInv <- by(sig2eps.draws[t-1,], INDICES=factor(1:J), FUN=function(x) solve(diag(1/sig2bO.draws[t-1,],pO)+(1/x)*t(SO)%*%SO)) 
    means.bO <- lapply(seq_along(1:J), FUN=function(j) unlist(VbOInv[[j]])%*%(1/sig2eps.draws[t-1,j])*t(SO)%*%
                         (Y[y.jlabs==j]- Sf%*%(betaf.draws[t,]) - SDf%*%betaDf.draws[[t]][[d.draws[t-1,j]]] - SD%*%bD.draws[[t]][[d.draws[t-1,j]]] - r.draws[t-1,] - rD.draws[[t]][[d.draws[t-1,j]]]))                     
    bO.draws[t,] <- unlist(lapply(seq_along(1:J), FUN=function(j) mvrnorm(1, mu=unlist(means.bO[[j]]), Sigma=unlist(VbOInv[[j]])) ))
    
    VrInv <- solve(1/sig2r.draws[t-1] + sum(1/sig2eps.draws[t-1,]))
    mean.diffr <- (Y - Fmat.F%*%(betaf.draws[t,]) - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]])  - Zmat.O%*%(bO.draws[t,]) - RDvec)
    mean.diffrlist <- by(mean.diffr, INDICES=factor(y.jlabs), FUN=function(x) 1*x)
    means.rbyout <-  lapply(seq(1:J), FUN=function(j) (1/sig2eps.draws[t-1,j])*unlist(mean.diffrlist[[j]]))
    means.r <- apply(matrix(unlist(means.rbyout),nrow=J,byrow=T),2,function(x) VrInv*sum(x))
    r.draws[t,] <- rnorm(n, mean=means.r, sd=sqrt(VrInv))
    
    Rvec <- t(t(kronecker(c(rep(1,J)),r.draws[t,]))) ## will NOT change with changing domains
    
    
    ### ---> Remaining variance parameters:
    bOj.sums <- apply(matrix(bO.draws[t,], nrow=pO),1,function(x) sum(x^2))
    sig2bO.draws[t,] <- rinvgamma(pO, shape=pr.A0.bO+J/2, scale=pr.B0.bO+.5*bOj.sums)
    
    sig2r.draws[t] <- rinvgamma(1, shape=pr.A0.r + n/2, scale=pr.B0.r + .5*t(r.draws[t,])%*%r.draws[t,])
    
    mean.diff2 <- (Y - Fmat.F%*%(betaf.draws[t,]) - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]])  - Zmat.O%*%(bO.draws[t,]) - RDvec - Rvec)
    sig2eps.draws[t,] <- rinvgamma(J, shape=pr.A0.eps + n/2, scale=pr.B0.eps + .5*aggregate(mean.diff2, by=list(outcomes=y.jlabs), FUN=function(x) sum(x^2))[,'V1'])
    ###########################################################################################################
    
    
    ###########################################################################################################
    ##### Clustering the Outcomes into domains:
    fixed.mean <- Fmat.F%*%(betaf.draws[t,]) + Zmat.O%*%(bO.draws[t,]) + Rvec 
    
    d.olddraws <- d.draws[(t-1),]
    
    if(t%%5==0){
      draw.d <- draw.JN.d( Y=Y, tune=tune, d.olddraws=d.olddraws, y.jlabs=y.jlabs, fixed.mean=fixed.mean, alpha=pr.alpha,
                           SD=SD, SDf=SDf, sig2eps=sig2eps.draws[t,],sig2rD=sig2rD.draws[[t]], sig2bD=sig2bD.draws[t,], 
                           r.D=rD.draws[[t]], b.D=bD.draws[[t]], beta.Df = betaDf.draws[[t]], 
                           pr.Sigma0InvDf=pr.Sigma0InvDf, pr.beta0Df=pr.beta0Df,pr.A0.rD=pr.A0.rD,pr.B0.rD=pr.B0.rD, m=m, printSumms=printSumms )    
      #print(draw.d$accepted)
      accepted[t] <- draw.d$accepted
    }else{
      draw.d <- draw.domA8( Y=Y, n=n, J=J, d.olddraws=d.olddraws, y.jlabs=y.jlabs, fixed.mean=fixed.mean, alpha=pr.alpha,
                            SD=SD, SDf=SDf, sig2eps=sig2eps.draws[t,],sig2rD=sig2rD.draws[[t]], sig2bD=sig2bD.draws[t,], 
                            r.D=rD.draws[[t]], b.D=bD.draws[[t]], beta.Df = betaDf.draws[[t]], 
                            pr.Sigma0InvDf=pr.Sigma0InvDf, pr.beta0Df=pr.beta0Df,pr.A0.rD=pr.A0.rD,pr.B0.rD=pr.B0.rD, 
			    m=m, printSumms=printSumms, t=t )    
      accepted[t] <- NA
     
    }
    
    d.draws[t,] <- draw.d$d.new
    K <- max(d.draws[t,])
    betaDf.draws[[t]] <- draw.d$betaDf.new
    bD.draws[[t]] <- draw.d$bD.new
    rD.draws[[t]] <- draw.d$rD.new 
    sig2rD.draws[[t]] <- draw.d$sig2rD.new 
    
    if(nlevels(factor(d.draws[t,]))>1){
      Lambda <- cbind(1*(d.draws[t,]==1),model.matrix( ~ as.factor(d.draws[t,])  ,contrasts='contr.sum')[,-1])  
    }else{
      ## If all outcomes grouped into the same domain, need to have just a vector of 1s
      Lambda <- d.draws[t,]
    }
    
    MeHg.coef[t,] <- betaf.draws[t,] + bO.draws[t,] + unlist(bD.draws[[t]][d.draws[t,]])
    
    niterT <- t
    
    if(t>min.niter){
      if(t%%200==0){
        ess.test <- multiESS(MeHg.coef[1:t,])
        # cat('\nbO.ess=',ess.test)
        # cat('\nESS>ess.test:',sum(ess.test>=bOESS))
        if(ess.test>=MeHg.ess){
           cat('\nConvergence based on',MeHg.ess,' is reached at t=',t)  
	         t <- niter           
        }
      }
    }
    t <- t + 1 
  }
  niterT.1 <- niterT
  ess.test1 <- ess.test
  ### Adaptive Step 1:
  ## Use Least-squares clustering to determine best fixed domain assigments for stage 2.
  
  burnin <- nburn
  d.drawsMinusBurnin <- d.draws[c((burnin+1):niterT),]
  K.draws <- apply(d.drawsMinusBurnin, 1, max)
  if(length(unique(K.draws))==1){
    if(max(unique(K.draws))!=1){
	d.drawsNewAllK <- renumber(d.drawsMinusBurnin, unique(K.draws))
    } else {
    	d.drawsNewAllK <-list(newds=d.drawsMinusBurnin, gppatterns=nrow(d.drawsMinusBurnin))
    }
  } else {
    d.drawsNewAllK <- renumber(d.drawsMinusBurnin)
  }
  ## Compute distance matrix, compute LS errors, pick clustering that minimizes LS:
  S.hat <- 1 - dissMatrix(d.drawsNewAllK$newds)
  LSclustresults <- dahlLSclust(d.drawsNewAllK$newds, S.hat)
  dom.assnsLS <- d.drawsNewAllK$newds[which(LSclustresults==min(LSclustresults))[1],]
  ## Set Kpost = final number of domains:
  Kpost <- max(dom.assnsLS)
  cat('The final posterior domain assignment (after least squares clustering) is')
  print(by(domAssn, dom.assnsLS, FUN=function(x) as.character(x)))
 
  ## Identify last iteration which had same domain assignment as dom.assnLS
  niterLast <- which(LSclustresults==min(LSclustresults))[length(which(LSclustresults==min(LSclustresults)))] + burnin
  
  ## Summaries for clustering to save to output (NOTE: This is after the burnin!!)
  #  (number of iterations matching the least squares clustering)
  #  (number of different patterns of size Kpost)
  #  (number of unique patterns visited after relabeling)
  niterMatch <- sum(apply(d.drawsNewAllK$newds, 1,function(x) sum(x==dom.assnsLS)==J))
  nKPatts <- sum(unlist(lapply(1:length(d.drawsNewAllK$gppatterns), 
                               function(x) max(d.drawsNewAllK$gppatterns[[x]])==Kpost)))
  nPatts <- length(d.drawsNewAllK$gppatterns)
  
  ## Use the last iteration's values which had the same domain assignment as starting poitns. This avoids a burnin.
  inits1 <-  list(betaF = c(betaf.draws[niterLast,]), bO = c(bO.draws[niterLast,]),r = c(r.draws[niterLast,]), 
                  sig2bO=sig2bO.draws[niterLast], sig2bD=sig2bD.draws[niterLast], sig2eps=sig2eps.draws[niterLast,], 
                  sig2rD=unlist(sig2rD.draws[[niterLast]]), sig2r=sig2r.draws[niterLast],
                  betaDf = betaDf.draws[[niterLast]], 
                  bD = bD.draws[[niterLast]], 
                  rD=rD.draws[[niterLast]])
  
  ## Clear out parameters with same names.
  rm(inits, d.draws, d.drawsNewAllK, betaDf.draws, betaf.draws, bO.draws, bD.draws, rD.draws, r.draws, sig2r.draws,
     sig2eps.draws, sig2rD.draws, sig2bD.draws, sig2bO.draws)
  cat('\nStarting Stage 2 with K=',Kpost,'domains.\n')
  ### Adaptive Step 1:
  clust1 <- fixDom.multouts( Y=Y, Sf=Sf, SDf=SDf, SD=SD, SO=SO, y.jlabs=y.jlabs, y.ilabs=y.ilabs, 
                            priors=priors, inits=inits1, min.niter=min.niter2, max.niter = max.niter, 
                            printIter=printIter, dom.assns=dom.assnsLS ) 
  
  ess.test2 <- clust1$ess.test
  niterT.2 <- clust1$niterT.2
  
  ## -- Do not need to remove burnin from all coefficient parameters since using initial values from last chains ending point.
  ## -- Place the parameter draws into matrices, selecting only the MCMC draws with "Kpost" domains
  betaDfdrawsK <- unlist(lapply(seq_along(1:niterT.2), FUN=function(k) matrix(unlist(clust1$betaDf.draws[[k]]),ncol=pDF*Kpost)))
  betaDfdrawsKmat <- matrix(betaDfdrawsK,nrow=niterT.2,byrow=T)
  colnames(betaDfdrawsKmat) <- paste(paste('betaDf',seq(1,pDF),sep=''),rep(1:Kpost,each=pDF),sep='_')
  #s2rDdrawsKmat <- clust1$sig2rD.draws[1:niterT.2,]
  
  ################################ START HERE ######################################################
  ### -- Place all coefficient parameter draws into 1 large matrix
  ##  We want the est, se, and 95% post int for these parameters
  param.drawsFull=cbind(clust1$betaf.draws, clust1$sig2bD.draws, clust1$sig2bO.draws)[1:niterT.2,]
  colnames(param.drawsFull) <- c('betaf','sig2bD', 'sig2bO')
  ##  We want the est, se for these parameters
  param.drawsRed=cbind(clust1$sig2r.draws[1:niterT.2], clust1$sig2eps.draws[1:niterT.2,], betaDfdrawsKmat, clust1$sig2rD.draws[1:niterT.2,])
  colnames(param.drawsRed) <- c('sig2r', paste('sig2eps',1:J,sep='.'), colnames(betaDfdrawsKmat), paste('sig2rD',1:Kpost,sep='.'))
  
  ### Parameter estimates, SEs, and 95% posterior intervals
  results.end1 <- summAll(mcmc.results=param.drawsFull, dom.assns=dom.assnsLS, full=TRUE )
  ### Parameter estimates, and SEs
  results.end2 <- summAll(mcmc.results=param.drawsRed, dom.assns=dom.assnsLS, full=FALSE )
  results.doms <- matrix(dom.assnsLS,nrow=1)
  colnames(results.doms) <- paste('d',1:length(dom.assnsLS),sep='.')
  
  cat('\nEnd of round****************\n')
  full.results <- data.frame(results.doms, 'ess.test1'=ess.test1, 'niterT.1'=niterT.1, 'Kpost'=Kpost, 'niterMatch'=niterMatch,
                             'nPatterns'=nPatts, 'nKpatterns'=nKPatts, 'ess.test2'=ess.test2, 
                             'niterT.2'=niterT.2, results.end1, results.end2)
  
  return(full.results)
}


##############################################################################################################
#####     FUNCTION TO DRAW THE DOMAIN INDICATOR VARIABLE IN NONCONJUGATE MODEL USING SPLIT-MERGE PROCEDURE 
##############################################################################################################  
new.goes2 <- function(x, y){
  if(max(x) == y) {
    return(2)  
  }else{
    return(1)  
  }
}


# sig2eps=sig2eps.draws[t,];sig2rD=sig2rD.draws[[t]];sig2bD=sig2bD.draws[t]
# r.D=rD.draws[[t]]; b.D=bD.draws[[t]]; beta.Df = betaDf.draws[[t]]
# pr.Sigma0InvDf; pr.beta0Df; pr.A0.rD; pr.B0.rD
draw.JN.d <- function(Y, tune,  d.olddraws, y.jlabs, fixed.mean, alpha, mean.diffnoDs, SD, SDf,
                      sig2eps=sig2eps.draws[t,],sig2rD=sig2rD.draws[[t]],sig2bD=sig2bD.draws[t,], 
                      r.D=rD.draws, b.D=bD.draws[[t]], beta.Df = betaDf.draws[[t]],
                      pr.Sigma0InvDf, pr.beta0Df, pr.A0.rD, pr.B0.rD, printSumms ){
  
  ## tune = (tuning param) # of intermediate RGSS to modify the launch_split and launch_merge states
  ## accepted: 1-accepted a split, 2-rejected a split, 3-accepted a merge, 4-rejected a merge
  
  K = max(d.olddraws)
  qniter <- tune[1]
  rniter <- tune[2]
  
  pDF <- ncol(SDf); pD <- ncol(SD)
  n <- nrow(SDf)
  J <- length(unique(y.jlabs))
  
  outs <- sample(1:J, size=2, replace=F)
  # outs=c(1,4); d_1=2; d_13=3
  d.Lsplit <- d.Lmerge <-  d.olddraws[which(d.olddraws==d.olddraws[outs[1]]|d.olddraws==d.olddraws[outs[2]])]
  ## S = set out outcomes in same group as outcome 1 and outcome 2 sampled in 'outs'
  S <- which(d.olddraws==d.olddraws[outs[1]]|d.olddraws==d.olddraws[outs[2]])[!which(d.olddraws==d.olddraws[outs[1]]|d.olddraws==d.olddraws[outs[2]])%in%(outs)] 
  ## S.ij = S + outcome i and j from 'outs'
  S.ij <- which(d.olddraws==d.olddraws[outs[1]]|d.olddraws==d.olddraws[outs[2]])
  ## Kcur = total current number of domains
  Kcur <- length(unique(d.olddraws))
  
  
  ################ Step 3 in Paper ###################
  ##### CONSTRUCT THE LAUNCH STATES:
  if(d.olddraws[outs[1]]==d.olddraws[outs[2]]){ 
    ### If i and j are in the same domain, 
    #    - split them into two by changing outcome 1 and keeping 2 the same for launch_split
    #    - make no changes for the launch_merge
    d.Lsplit[S.ij==outs[1]] <-  di_new <-  Kcur + 1 # need for launch_split
    dj <- d.olddraws[outs[2]]
    if(printSumms==TRUE) cat('Split: outs=',outs,' dij=',c(di_new,dj))
    
  }else{
    ### If i and j are in separate domains, 
    #    - let them keep their original labels for launch_split
    #    - assign all outcomes to the j^th domain (outcome 2) for launch_merge
    di_new <- d.olddraws[outs[1]] # need for launch_split
    dj <- d.olddraws[outs[2]]
    d.Lmerge <- rep(dj, length(S.ij))
    if(printSumms==TRUE) cat('Merge: outs=',outs,' dij=',c(di_new,dj))
    
  }
  ## Allocate the remaining outcomes in S to either of the two components w.p. 1/2
  d.Lsplit[S.ij%in%S] <- sample(c(di_new,dj), size=sum(S.ij%in%S),replace=T)
  
  # set.seed(97)
  ## Draw new parameters from the PRIOR for each component of luanch_split and only 1 for the launch_merge:
  # (1) sig2rD.new ~ IG(a, b)Ind(0,10) where a and b are hyperparameters and Ind(0,10) is the truncation
#   sig2rD.Lsplit <- list(rinvgamma(1, shape=pr.A0.rD, scale=pr.B0.rD), rinvgamma(1, shape=pr.A0.rD, scale=pr.B0.rD))
#   sig2rD.Lmerge <- rinvgamma(1, shape=pr.A0.rD, scale=pr.B0.rD)
  sig2rD.Lsplit <- list(sig2rD[[dj]], sig2rD[[dj]])
  sig2rD.Lmerge <- sig2rD[[dj]]
  
  # (2) rD.Lsplit|sig2rD.Lsplit ~ N(0, sqrt(sig2rD.Lsplit)) 
#   rD.Lsplit <- list(rnorm(n, mean=0, sd=sqrt(sig2rD.Lsplit[[1]])), rnorm(n, mean=0, sd=sqrt(sig2rD.Lsplit[[2]])))
#   rD.Lmerge <- rnorm(n, mean=0, sd=sqrt(sig2rD.Lmerge))
  rD.Lsplit <- list(r.D[[dj]], r.D[[dj]])
  rD.Lmerge <- r.D[[dj]]
  
  # (3) bD.Lsplit ~ N(0, sqrt(sig2bD.curr)) 
#   bD.Lsplit <- list(rnorm(pD, mean=0, sd=sqrt(sig2bD)), rnorm(pD, mean=0, sd=sqrt(sig2bD)))
#   bD.Lmerge <- rnorm(pD, mean=0, sd=sqrt(sig2bD) )
  bD.Lsplit <- list(b.D[[dj]], b.D[[dj]])
  bD.Lmerge <- b.D[[dj]]
  
  # (4) betaDF.Lsplit ~ N(pr.beta0Df, pr.Sigma0InvDf) -- independent in the prior
#   betaDf.Lsplit <- list(matrix(rmvnorm(1, mean=pr.beta0Df, sigma = solve(pr.Sigma0InvDf)),ncol=1),
#                         matrix(rmvnorm(1, mean=pr.beta0Df, sigma = solve(pr.Sigma0InvDf)),ncol=1))
#   betaDf.Lmerge <- matrix(rmvnorm(1, mean=pr.beta0Df, sigma = solve(pr.Sigma0InvDf)),ncol=1)
  betaDf.Lsplit <- list(beta.Df[[dj]],beta.Df[[dj]])
  betaDf.Lmerge <- beta.Df[[dj]]
  
  
  #### UPDATE THE LAUNCH_SPLIT STATES WITH q AND r RESTRICTED GIBBS SAMPLING SCANS (RGSS):
  meandiffs <- Y - fixed.mean
  
  d.ij <- c(di_new, dj)
  for(i in 1:qniter){
    # i = 1
    # i = i + 1
    # cat('\ni=',i,'dsplit=',d.Lsplit)
    # This will force outcome i's domain to always be in the first component of the parameter list 
    #  (recall that in launch split if d_i=d_j originally then d_i = K+1)
    
    
    meansDfLsplit <- lapply(seq_along(1:2), FUN=function(k) pr.Sigma0InvDf%*%pr.beta0Df + t(kronecker(rep(1,sum(d.Lsplit==d.ij[k])),SDf))%*%
                              diag(rep(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]],each=n))%*%(matrix(meandiffs[y.jlabs%in%S.ij[d.Lsplit==d.ij[k]]],ncol=1)
                                                                                      - kronecker(rep(1,sum(d.Lsplit==d.ij[k])),(SD%*%(bD.Lsplit[[k]]))) - matrix(kronecker(rep(1,sum(d.Lsplit==d.ij[k])),(rD.Lsplit[[k]])),ncol=1 ))  )                                                       
    vDfLsplit <- lapply(seq_along(1:2), FUN=function(k) solve(pr.Sigma0InvDf + sum(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]])*t(SDf)%*%SDf)) ## optimize other sampler!!!
    betaDf.Lsplit <- lapply(X=seq_along(1:2), FUN = function(k) mvrnorm(1, mu=unlist(vDfLsplit[[k]])%*%unlist(meansDfLsplit[[k]]) , Sigma =unlist(vDfLsplit[[k]])))  
    
    meansbDLsplit <- lapply(seq_along(1:2), FUN=function(k) t(kronecker(rep(1,sum(d.Lsplit==d.ij[k])),SD))%*%diag(rep(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]],each=n))%*%
                              (matrix(meandiffs[y.jlabs%in%S.ij[d.Lsplit==d.ij[k]]],ncol=1)
                               - kronecker(rep(1,sum(d.Lsplit==d.ij[k])),(SDf%*%(betaDf.Lsplit[[k]]))) - matrix(kronecker(rep(1,sum(d.Lsplit==d.ij[k])),(rD.Lsplit[[k]])),ncol=1 ))  )                                                       
    vbDLsplit<- lapply(seq_along(1:2), FUN=function(k) solve(diag(1/sig2bD,pD) + sum(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]])*t(SD)%*%SD) )
    bD.Lsplit <- lapply(X=seq_along(1:2), FUN = function(k) rnorm(pD, mean=unlist(vbDLsplit[[k]])%*%unlist(meansbDLsplit[[k]]) , sd = sqrt(vbDLsplit[[k]])))  
    
    sig2rD.Lsplit <- lapply(seq(1:2), FUN=function(k)  
      rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(unlist(rD.Lsplit[[k]]))%*%unlist(rD.Lsplit[[k]])))
    
    vrDLsplit <- lapply(seq_along(1:2), FUN=function(k) solve((1/sig2rD.Lsplit[[k]]) + sum(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]])) )
    separ.rDsplit <- lapply(seq_along(1:length(S.ij)), FUN=function(j) (1/sig2eps[S.ij[j]])*(matrix(meandiffs[y.jlabs%in%S.ij[j]],ncol=1) -
                                                                                               SDf%*%(betaDf.Lsplit[[new.goes2(d.Lsplit,d.Lsplit[j])]]) - SD%*%(bD.Lsplit[[new.goes2(d.Lsplit,d.Lsplit[j])]])) )
    means.rDsplitmat <- matrix(unlist(separ.rDsplit),ncol=length(S.ij),byrow=F)
    means.rDLsplit <- lapply(seq_along(1:2), FUN=function(k)
      if( sum(d.Lsplit==d.ij[k])==1){
        means.rDsplitmat[,k] 
      }else{
        rowSums(means.rDsplitmat[,d.Lsplit==d.ij[k]]) 
      } )
    rD.Lsplit <- lapply(seq_along(1:2), FUN=function(k) rnorm(n, mean=unlist(vrDLsplit[[k]])%*%(unlist(means.rDLsplit[[k]]) ), sd=sqrt(unlist(vrDLsplit[[k]])))) 
    
    
    if(length(d.Lsplit)>2){
      lmeans <- lapply(seq_along(1:2),FUN=function(k) 
        matrix(fixed.mean[y.jlabs%in%S],ncol=1 ) + kronecker(rep(1,length(S)),SDf)%*%betaDf.Lsplit[[k]] + 
          kronecker(rep(1,length(S)),SD)%*%bD.Lsplit[[k]] + matrix(kronecker(rep(1,length(S)),rD.Lsplit[[k]]),ncol=1) )
      lsds <- rep(sqrt(sig2eps[S]),each=n)
      ldens <- lapply(seq_along(1:2), FUN=function(k) by(dnorm(Y[y.jlabs%in%S],mean=lmeans[[k]],sd=lsds,log=T),  
                                                         INDICES = factor(y.jlabs[y.jlabs%in%S]), FUN = sum))
      
      phi <- lapply(seq_along(S), FUN=function(k) -max(ldens[[1]][k],ldens[[2]][k]))
      dens <- lapply(seq_along(1:2), FUN = function(k) exp(ldens[[k]]+unlist(phi)))
      
      probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
      # cat(' launch split: probi=',probi)
      d.Lsplit[S.ij%in%S] <- unlist(lapply(seq_along(S), FUN=function(x) sample(d.ij, size=1, prob=c(probi[[x]],1-probi[[x]]))))
    }
    ########### DO I NEED AN ELSE HERE -- e.g., if there are only two outcomes do they stay separate or is
    ###########  there a possibility that they can merge even in the launch split state...?
  }
  
  for(i in 1:rniter){
    if(K==1){
      vDfLmerge <- solve(pr.Sigma0InvDf + Reduce('+',lapply(1:length(S.ij),FUN=function(k) (1/sig2eps[S.ij[k]])*t(SDf)%*%SDf)))
      meanbetaDf.Lmerge <- pr.Sigma0InvDf%*%pr.beta0Df + t(SDf)%*%
        Reduce('+',lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                             by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij], INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]]))  
      betaDf.Lmerge <- matrix(rmvnorm(1, mean=vDfLmerge%*%meanbetaDf.Lmerge, sigma=vDfLmerge),ncol=1)
      
    }else{
      sig2rD.Lmerge <- rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge)
      
      vrD.Lmerge <- solve((1/sig2rD.Lmerge) + sum(1/sig2eps[S.ij]))
      meanrD.Lmerge <-  lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                                  by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SDf)%*%betaDf.Lmerge - t(t(rep(SD,length(S.ij))))%*%bD.Lmerge ,
                                     INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]])
      rD.Lmerge <- rnorm(n, mean=vrD.Lmerge%*%rowSums(matrix(unlist(meanrD.Lmerge),nrow=n,byrow=F)), sd=sqrt(vrD.Lmerge))
      
      vbD.Lmerge <- solve(1/sig2bD + sum(apply(matrix((1/sig2eps[S.ij]),ncol=length(S.ij)),2,function(x) x*t(SD)%*%SD)))
      meanbD.Lmerge <- t(SD)%*%rowSums(matrix(unlist(lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                                                               by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SDf)%*%betaDf.Lmerge - matrix(kronecker(rep(1,length(S.ij)),rD.Lmerge),ncol=1) , INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]])),nrow=n,byrow=F))    
      bD.Lmerge <- rnorm(pD, mean=vbD.Lmerge%*%meanbD.Lmerge, sd = sqrt(vbD.Lmerge))
      
      vDfLmerge <- solve(pr.Sigma0InvDf + Reduce('+',lapply(1:length(S.ij),FUN=function(k) (1/sig2eps[S.ij[k]])*t(SDf)%*%SDf)))
      meanbetaDf.Lmerge <- pr.Sigma0InvDf%*%pr.beta0Df + t(SDf)%*%
        Reduce('+',lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                             by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SD)%*%bD.Lmerge - matrix(kronecker(rep(1,length(S.ij)),rD.Lmerge),ncol=1) , INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]]))  
      betaDf.Lmerge <- matrix(rmvnorm(1, mean=vDfLmerge%*%meanbetaDf.Lmerge, sigma=vDfLmerge),ncol=1)
    }
  }
  
  
  
  if(d.olddraws[outs[1]]==d.olddraws[outs[2]]){
    ### PROPOSE A SPLIT:
    d.split <- d.Lsplit
    ## Assign all k in S to either K+1 or based on launch_split 
    ##  We can re-use the final probability calculated in the launch_split procedure:   lmeans.split <- lapply(seq_along(1:2),FUN=function(k) matrix(fixed.mean[y.jlabs%in%S],ncol=1 ) + kronecker(rep(1,length(S)),SDf)%*%betaDf.Lsplit[[k]] + kronecker(rep(1,length(S)),SD)%*%bD.Lsplit[[k]] + matrix(kronecker(rep(1,length(S)),rD.Lsplit[[k]]),ncol=1) )
    #    probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
    
    if(length(d.Lsplit)>2){
      d.split[S.ij%in%S] <- unlist(lapply(seq_along(S), FUN=function(x) sample(d.ij, size=1, prob=c(probi[[x]],1-probi[[x]]))))
    }
    
    meansDfsplit <- lapply(seq_along(1:2), FUN=function(k) pr.Sigma0InvDf%*%pr.beta0Df + t(kronecker(rep(1,sum(d.split==d.ij[k])),SDf))%*%
                             diag(rep(1/sig2eps[S.ij[d.split==d.ij[k]]],each=n))%*%(matrix(meandiffs[y.jlabs%in%S.ij[d.split==d.ij[k]]],ncol=1)
                                                                                    - kronecker(rep(1,sum(d.split==d.ij[k])),(SD%*%(bD.Lsplit[[k]]))) 
                                                                                    - matrix(kronecker(rep(1,sum(d.split==d.ij[k])),(rD.Lsplit[[k]])),ncol=1 ))  )                                                       
    vDfsplit <- lapply(seq_along(1:2), FUN=function(k) solve(pr.Sigma0InvDf + sum(1/sig2eps[S.ij[d.split==d.ij[k]]])*t(SDf)%*%SDf)) ## optimize other sampler!!!
    betaDf.split <- lapply(X=seq_along(1:2), FUN = function(k) mvrnorm(1, mu=unlist(vDfsplit[[k]])%*%unlist(meansDfsplit[[k]]) , Sigma =unlist(vDfsplit[[k]])))  
    
    meansbDsplit <- lapply(seq_along(1:2), FUN=function(k) t(kronecker(rep(1,sum(d.split==d.ij[k])),SD))%*%diag(rep(1/sig2eps[S.ij[d.split==d.ij[k]]],each=n))%*%
                             (matrix(meandiffs[y.jlabs%in%S.ij[d.split==d.ij[k]]],ncol=1)
                              - kronecker(rep(1,sum(d.split==d.ij[k])),(SDf%*%(betaDf.split[[k]]))) - matrix(kronecker(rep(1,sum(d.split==d.ij[k])),(rD.Lsplit[[k]])),ncol=1 ))  )                                                       
    vbDsplit<- lapply(seq_along(1:2), FUN=function(k) solve(diag(1/sig2bD,pD) + sum(1/sig2eps[S.ij[d.split==d.ij[k]]])*t(SD)%*%SD) )
    bD.split <- lapply(X=seq_along(1:2), FUN = function(k) rnorm(pD, mean=unlist(vbDsplit[[k]])%*%unlist(meansbDsplit[[k]]) , sd = sqrt(vbDsplit[[k]])))  
    
    
    
    vrDsplit <- lapply(seq_along(1:2), FUN=function(k) solve((1/sig2rD.Lsplit[[k]]) + sum(1/sig2eps[S.ij[d.split==d.ij[k]]])) )
    separ.rDsplit <- lapply(seq_along(1:length(S.ij)), FUN=function(j) (1/sig2eps[S.ij[j]])*(matrix(meandiffs[y.jlabs%in%S.ij[j]],ncol=1) -
                                                                                               SDf%*%(betaDf.split[[new.goes2(d.split,d.split[j])]]) - SD%*%(bD.split[[new.goes2(d.split,d.split[j])]]))) 
    means.rDsplitmat <- matrix(unlist(separ.rDsplit),ncol=length(S.ij),byrow=F)
    means.rDsplit <- lapply(seq_along(1:2), FUN=function(k)
      if( sum(d.split==d.ij[k])==1){
        means.rDsplitmat[,k] 
      }else{
        rowSums(means.rDsplitmat[,d.split==d.ij[k]]) 
      } )
    rD.split <- lapply(seq_along(1:2), FUN=function(k) rnorm(n, mean=unlist(vrDsplit[[k]])%*%(unlist(means.rDsplit[[k]]) ), sd=sqrt(unlist(vrDsplit[[k]])))) 
    
    sig2rD.split <- lapply(seq(1:2), FUN=function(k)  
      rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(unlist(rD.split[[k]]))%*%unlist(rD.split[[k]])))
    
    
    if(printSumms==TRUE) cat('betadf=',dmvnorm(beta.Df[[dj]], mean=betaDf.Lmerge, sigma=vDfLmerge, log=T))
    if(printSumms==TRUE & K>1) cat('\nbd=', dnorm(b.D[[dj]], mean=bD.Lmerge, sd=sqrt(vbD.Lmerge), log=T),
                                   '\nrd=' ,sum(dnorm(r.D[[dj]], mean=rD.Lmerge, sd=sqrt(vrD.Lmerge),log=T)),
                                   '\ns2rd=', actuar::dinvgamma(sig2rD[[dj]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge,log=T ))
    
    if(K==1){
      lqcurr.split <- dmvnorm(beta.Df[[dj]], mean=betaDf.Lmerge, sigma=vDfLmerge, log=T) 
      lpr.curr <- (log(factorial(length(S.ij) -1) )
                   + dmvnorm(beta.Df[[dj]],  mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf), log=T) ) 
    }else{
      lqcurr.split <- (dmvnorm(beta.Df[[dj]], mean=betaDf.Lmerge, sigma=vDfLmerge, log=T) + dnorm(b.D[[dj]], mean=bD.Lmerge, sd=sqrt(vbD.Lmerge), log=T)
                       + sum(dnorm(r.D[[dj]], mean=rD.Lmerge, sd=sqrt(vrD.Lmerge),log=T)) + actuar::dinvgamma(sig2rD[[dj]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge,log=T ) )
      
      lpr.curr <- (log(factorial(length(S.ij) -1) )
                   + dmvnorm(beta.Df[[dj]],  mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf), log=T) + dnorm(b.D[[dj]], mean=0, sd=sqrt(sig2bD), log=T)
                   + sum(dnorm(r.D[[dj]],mean=rep(0,n), sd=sqrt(sig2rD[[dj]]), log=T) ) + actuar::dinvgamma(sig2rD[[dj]], shape=pr.A0.rD , scale=pr.B0.rD, log=T ) ) 
    }
    lpr.split <- (log(alpha) + log(factorial(sum(d.split==d.ij[1]) -1) + factorial(sum(d.split==d.ij[2]) -1) ) 
                  + Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(betaDf.split[[k]], mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf), log=T) 
                                       + dnorm(bD.split[[k]], mean=0, sd=sqrt(sig2bD), log=T)
                                       + sum(dnorm(rD.split[[k]], mean=rep(0,n), sd=sqrt(sig2rD.split[[k]]), log=T) )
                                       + actuar::dinvgamma(sig2rD.split[[k]], shape=pr.A0.rD, scale=pr.B0.rD , log=T) )) )
    
    if(length(d.Lsplit)>2){
      splmeans <- lapply(seq_along(1:2),FUN=function(k) matrix(fixed.mean[y.jlabs%in%S],ncol=1 ) + kronecker(rep(1,length(S)),SDf)%*%betaDf.split[[k]] + kronecker(rep(1,length(S)),SD)%*%bD.split[[k]] + matrix(kronecker(rep(1,length(S)),rD.split[[k]]),ncol=1) )
      splsds <- rep(sqrt(sig2eps[S]),each=n)
      spldens <- lapply(seq_along(1:2), FUN=function(k) by(dnorm(Y[y.jlabs%in%S],mean=splmeans[[k]],sd=splsds,log=T),  INDICES = factor(y.jlabs[y.jlabs%in%S]), FUN = sum))
      
      phi <- lapply(seq_along(S), FUN=function(k) -max(spldens[[1]][k],spldens[[2]][k]))
      dens <- lapply(seq_along(1:2), FUN = function(k) exp(spldens[[k]]+unlist(phi)))
      
      probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
      
      probs <- list(probi, (1-probi))
      
      lqsplit.curr <- (Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(betaDf.split[[k]], mean=betaDf.Lsplit[[k]], sigma=vDfLsplit[[k]], log=T) 
                                          + dnorm(bD.split[[k]], mean=bD.Lsplit[[k]], sd=sqrt(vbDLsplit[[k]]), log=T)
                                          + sum(dnorm(rD.split[[k]], mean=rD.Lsplit[[k]], sd=sqrt(vrDLsplit[[k]]), log=T) )
                                          + actuar::dinvgamma(sig2rD.split[[k]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[k]])%*%rD.Lsplit[[k]], log=T ) ))
                       + Reduce(sum, lapply(seq(1,length(S)), FUN = function(k) log(probs[[apply(matrix(d.split[S.ij%in%S],nrow=1),2 , function(x) which(d.ij==x))[k]]][[k]])
                       ) ) )
    } else {
      lqsplit.curr <- Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(betaDf.split[[k]], mean=betaDf.Lsplit[[k]], sigma=vDfLsplit[[k]], log=T) 
                                         + dnorm(bD.split[[k]], mean=bD.Lsplit[[k]], sd=sqrt(vbDLsplit[[k]]), log=T)
                                         + sum(dnorm(rD.split[[k]], mean=rD.Lsplit[[k]], sd=sqrt(vrDLsplit[[k]]), log=T) )
                                         + actuar::dinvgamma(sig2rD.split[[k]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[k]])%*%rD.Lsplit[[k]], log=T ) ))
      
    }               
    
    llik.curr <- sum(dnorm(Y[y.jlabs%in%S.ij], mean=matrix(fixed.mean[y.jlabs%in%S.ij],ncol=1 ) + kronecker(rep(1,length(S.ij)),SDf)%*%beta.Df[[dj]] 
                           + kronecker(rep(1,length(S.ij)),SD)%*%b.D[[dj]] + matrix(kronecker(rep(1,length(S.ij)),r.D[[dj]]),ncol=1), 
                           sd = rep(sqrt(sig2eps[S.ij]),each=n) ,log=T))
    
    newlabsdij <- as.numeric(factor(d.split,levels=d.ij))
    
    llik.split <-  Reduce(sum, lapply(seq(1,2), FUN = function(k)  dnorm(Y[y.jlabs[y.jlabs%in%(S.ij[newlabsdij==k])]], mean=matrix(fixed.mean[y.jlabs[y.jlabs%in%(S.ij[newlabsdij==k])]],ncol=1 ) + 
                                                                           kronecker(rep(1,sum(newlabsdij==k)),SDf)%*%betaDf.split[[k]]
                                                                         + kronecker(rep(1,sum(newlabsdij==k)),SD)%*%bD.split[[k]] + matrix(kronecker(rep(1,sum(newlabsdij==k)),rD.split[[k]]),ncol=1), 
                                                                         sd = rep(sqrt(sig2eps[unique(y.jlabs[y.jlabs%in%(S.ij[newlabsdij==k])])]),each=n) ,log=T)))
    if(printSumms==TRUE) cat(' lqcur.spl=',lqcurr.split,'\n lqspl.cur=',lqsplit.curr,'\n lpr.spl=',lpr.split,'\n lpr.cur=',lpr.curr,'\n llik.spl=',llik.split, '\n llik.cur=',llik.curr)
    if(is.infinite(exp((lqcurr.split - lqsplit.curr) + (lpr.split - lpr.curr) + (llik.split - llik.curr)))){
      cat('Infinite Value in Acceptance Probability:')
      if(exp((lqcurr.split - lqsplit.curr) + (lpr.split - lpr.curr) + (llik.split - llik.curr))==-Inf){
        cat('(-Inf) --> Do NOT accept.')
      } else {
        cat('(Inf) --> ALWAYS accept.')
      }
      cat('\n lqcurr.spl=',lqcurr.split,'\n lqspl.cur=',lqsplit.curr,'\n lpr.spl=',lpr.split,'\n lpr.cur=',lpr.curr,'\n llik.spl=',llik.split, '\n llik.cur=',llik.curr)
      cat('\nbetadf_split=');print(betaDf.split); cat('\nbetadf_LMerge='); print(betaDf.Lmerge)
      cat('\nbD_split=');print(bD.split);cat('\nbD_Lmerge=');print(bD.Lmerge)
      cat('\ns2rd_slpit=');print(sig2rD.split);cat('\nsig2rd_Lmerge=');print(sig2rD.Lmerge)
    }
    
    
    if(is.infinite(exp((lqcurr.split - lqsplit.curr) + (lpr.split - lpr.curr) + (llik.split - llik.curr)))){
	acc.split <- 0
    }else{
	acc.split <- min( 1,  exp((lqcurr.split - lqsplit.curr) + (lpr.split - lpr.curr) + (llik.split - llik.curr)) )
    }
    u <- runif(1)
    if(u < acc.split){
      accepted <- 1
      d.olddraws[S.ij] <- d.split
      # Assign the 2nd component (dj) to the previous dj component and place new di at the end:
      beta.Df[[dj]] <- betaDf.split[[2]]; beta.Df[[length(beta.Df)+1]] <- betaDf.split[[1]]
      b.D[[dj]] <- bD.split[[2]]; b.D[[length(b.D)+1]] <- bD.split[[1]]
      r.D[[dj]] <- rD.split[[2]]; r.D[[length(r.D)+1]] <- rD.split[[1]]
      sig2rD[[dj]] <- sig2rD.split[[2]]; sig2rD <- c(sig2rD, sig2rD.split[[1]])
    } else{
      accepted <- 2
      d.olddraws <- d.olddraws
    }
    
  } else {
    
    ### PROPOSE A MERGE:
    ## Assign all k in S.ij to dj as in launch_merge 
    d.merge <- d.Lmerge
    
    ## Perform one final update for phi_merge
    sig2rD.merge <- rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge)
    # while(sig2rD.merge>500) sig2rD.merge <- rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge)
    vrDmerge <- solve((1/sig2rD.merge) + sum(1/sig2eps[S.ij]))
    meanrDmerge <-  lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                              by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SDf)%*%betaDf.Lmerge - t(t(rep(SD,length(S.ij))))%*%bD.Lmerge , INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]])
    rD.merge <- rnorm(n, mean=vrDmerge%*%rowSums(matrix(unlist(meanrDmerge),nrow=n,byrow=F)), sd=sqrt(vrDmerge))
    
    vbDmerge <- solve(1/sig2bD + sum(apply(matrix((1/sig2eps[S.ij]),ncol=length(S.ij)),2,function(x) x*t(SD)%*%SD)))
    meanbDmerge <- t(SD)%*%rowSums(matrix(unlist(lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                                                           by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SDf)%*%betaDf.Lmerge - matrix(kronecker(rep(1,length(S.ij)),rD.merge),ncol=1) , INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]])),nrow=n,byrow=F))    
    bD.merge <- rnorm(pD, mean=vbDmerge%*%meanbDmerge, sd = sqrt(vbDmerge))
    
    vbetaDfmerge <- solve(pr.Sigma0InvDf + Reduce('+',lapply(1:length(S.ij),FUN=function(k) (1/sig2eps[S.ij[k]])*t(SDf)%*%SDf)))
    meanbetaDfmerge <- pr.Sigma0InvDf%*%pr.beta0Df + t(SDf)%*%
      Reduce('+',lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                           by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SD)%*%bD.merge - matrix(kronecker(rep(1,length(S.ij)),rD.merge),ncol=1) , INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]]))  
    betaDf.merge <- matrix(rmvnorm(1, mean=vbetaDfmerge%*%meanbetaDfmerge, sigma=vbetaDfmerge),ncol=1)
    
    if(length(d.Lsplit)>2){
      currmeans <- lapply(seq_along(1:2),FUN=function(k) matrix(fixed.mean[y.jlabs%in%S],ncol=1 ) + kronecker(rep(1,length(S)),SDf)%*%beta.Df[[d.ij[k]]] + kronecker(rep(1,length(S)),SD)%*%b.D[[d.ij[k]]] + matrix(kronecker(rep(1,length(S)),r.D[[d.ij[k]]]),ncol=1) )
      currsds <- rep(sqrt(sig2eps[S]),each=n)
      lcurrdens <- lapply(seq_along(1:2), FUN=function(k) by(dnorm(Y[y.jlabs%in%S],mean=currmeans[[k]],sd=currsds,log=T),  INDICES = factor(y.jlabs[y.jlabs%in%S]), FUN = sum))
      
      phi <- lapply(seq_along(S), FUN=function(k) -max(lcurrdens[[1]][k],lcurrdens[[2]][k]))
      dens <- lapply(seq_along(1:2), FUN = function(k) exp(lcurrdens[[k]]+unlist(phi)))
      
      probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
      
      probs <- list(probi, (1-probi))     
      if(printSumms==TRUE) cat('probi=', probi)
      lqcurr.merge <- (Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(beta.Df[[d.ij[k]]], mean=betaDf.Lsplit[[k]], sigma=vDfLsplit[[k]], log=T) 
                                          + dnorm(b.D[[d.ij[k]]], mean=bD.Lsplit[[k]], sd=sqrt(vbDLsplit[[k]]), log=T)
                                          + sum(dnorm(r.D[[d.ij[k]]], mean=rD.Lsplit[[k]], sd=sqrt(vrDLsplit[[k]]), log=T) )
                                          + actuar::dinvgamma(sig2rD[[d.ij[k]]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[k]])%*%rD.Lsplit[[k]],log=T ) ))
                       + Reduce(sum, lapply(seq(1,length(S)), FUN = function(k) log(probs[[apply(matrix(d.olddraws[S],nrow=1),2 , function(x) which(d.ij==x))[k]]][[k]])
                       ) ))
    } else {
      lqcurr.merge <- (Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(beta.Df[[d.ij[k]]], mean=betaDf.Lsplit[[k]], sigma=vDfLsplit[[k]], log=T) 
                                          + dnorm(b.D[[d.ij[k]]], mean=bD.Lsplit[[k]], sd=sqrt(vbDLsplit[[k]]), log=T)
                                          + sum(dnorm(r.D[[d.ij[k]]], mean=rD.Lsplit[[k]], sd=sqrt(vrDLsplit[[k]]), log=T) )
                                          + actuar::dinvgamma(sig2rD[[d.ij[k]]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[k]])%*%rD.Lsplit[[k]],log=T ) )))
    }
    
    if(printSumms==TRUE) cat('betadfi=',dmvnorm(beta.Df[[d.ij[1]]], mean=betaDf.Lsplit[[1]], sigma=vDfLsplit[[1]], log=T),'\nbetadfj=',dmvnorm(beta.Df[[d.ij[2]]], mean=betaDf.Lsplit[[1]], sigma=vDfLsplit[[1]], log=T),
                             '\nbdi=', dnorm(b.D[[d.ij[1]]], mean=bD.Lsplit[[1]], sd=sqrt(vbDLsplit[[1]]), log=T),'\nbdj=', dnorm(b.D[[d.ij[2]]], mean=bD.Lsplit[[2]], sd=sqrt(vbDLsplit[[2]]), log=T),
                             '\nrdi=',sum(dnorm(r.D[[d.ij[1]]], mean=rD.Lsplit[[2]], sd=sqrt(vrDLsplit[[1]]), log=T)) ,'\nrdj=' ,sum(dnorm(r.D[[d.ij[2]]], mean=rD.Lsplit[[2]], sd=sqrt(vrDLsplit[[2]]), log=T)) ,
                             '\ns2rdi=', actuar::dinvgamma(sig2rD[[d.ij[1]]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[1]])%*%rD.Lsplit[[1]],log=T ),'\ns2rdi=', actuar::dinvgamma(sig2rD[[d.ij[2]]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[2]])%*%rD.Lsplit[[2]],log=T ))
    
    lqmerge.curr <- (dmvnorm(t(betaDf.merge), mean=t(betaDf.Lmerge), sigma=vDfLmerge, log=T) + dnorm(bD.merge, mean=bD.Lmerge, sd=sqrt(vbD.Lmerge), log=T)
                     + sum(dnorm(rD.merge, mean=rD.Lmerge, sd=sqrt(vrD.Lmerge))) + actuar::dinvgamma(sig2rD.merge, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge,log=T ) )
    
    lpr.merg <- - log(alpha) +log(factorial(length(S.ij) -1) ) + (dmvnorm(t(betaDf.merge), mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf), log=T) 
                                                                  + dnorm(bD.merge, mean=0, sd=sqrt(sig2bD), log=T) + sum(dnorm(rD.merge, mean=rep(0,n), sd=sqrt(sig2rD.merge), log=T)) 
                                                                  + actuar::dinvgamma(sig2rD.merge, shape=pr.A0.rD, scale=pr.B0.rD , log=T) )
    lpr.curr <- (log(factorial(sum(d.olddraws==d.ij[1])-1) + factorial(sum(d.olddraws==d.ij[2])-1) )
                 + Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(beta.Df[[d.ij[k]]],  mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf), log=T)
                                      + dnorm(b.D[[d.ij[k]]], mean=0, sd=sqrt(sig2bD), log=T)
                                      + sum(dnorm(r.D[[d.ij[k]]],mean=rep(0,n), sd=sqrt(sig2rD[[d.ij[k]]]), log=T) )
                                      + actuar::dinvgamma(sig2rD[[d.ij[k]]], shape=pr.A0.rD , scale=pr.B0.rD, log=T ) )) )
    llik.merg <- sum(dnorm(Y[y.jlabs%in%S.ij], mean=matrix(fixed.mean[y.jlabs%in%S.ij],ncol=1 ) + kronecker(rep(1,length(S.ij)),SDf)%*%betaDf.merge 
                           + kronecker(rep(1,length(S.ij)),SD)%*%bD.merge + matrix(kronecker(rep(1,length(S.ij)),rD.merge),ncol=1), 
                           sd = rep(sqrt(sig2eps[S.ij]),each=n) ,log=T))
    
    llik.curr <-  Reduce(sum, lapply(seq(1,length(S.ij)), FUN = function(k)  dnorm(Y[y.jlabs==S.ij[k]], mean=matrix(fixed.mean[y.jlabs==S.ij[k]],ncol=1 ) + SDf%*%matrix(unlist(beta.Df[d.olddraws[S.ij[k]]]),ncol=1)
                                                                                   + SD%*%unlist(b.D[d.olddraws[S.ij[k]]]) + unlist(r.D[d.olddraws[S.ij[k]]]), 
                                                                                   sd = rep(sqrt(sig2eps[S.ij[k]]),each=n) ,log=T)))
    
    if(printSumms==TRUE) cat(' lqcur.mer=',lqcurr.merge,' lqmer.cur=',lqmerge.curr,' lpr.mer=',lpr.merg,' lpr.cur=',lpr.curr,' llik.mer=',llik.merg, ' llik.cur=',llik.curr)
    
    acc.merge <- min( 1,  exp((lqcurr.merge - lqmerge.curr) + (lpr.merg - lpr.curr) + (llik.merg - llik.curr)) )
    
    u <- runif(1)
    if(u < acc.merge){
      accepted <- 3
      d.olddraws[S.ij] <- d.merge
      d.olddraws <- as.numeric(factor(d.olddraws))
      beta.Df <- beta.Df[-d.ij[1]]
      b.D <- b.D[-d.ij[1]]
      r.D <- r.D[-d.ij[1]]
      sig2rD <- sig2rD[-d.ij[1]]
    } else{
      accepted <- 4
      d.olddraws <- d.olddraws
    }
    
  }
  #print(accepted)
  return( list(d.new=d.olddraws, betaDf.new=beta.Df, bD.new=b.D, rD.new=r.D, sig2rD.new=sig2rD, accepted=accepted)) 
}

###########################################################################################################
#####     FUNCTION TO DRAW THE DOMAIN INDICATOR VARIABLE IN NONCONJUGATE MODEL USING GIBBS SAMPLING 
###########################################################################################################  

draw.domA8 <- function( Y, n, J, d.olddraws, y.jlabs, fixed.mean, alpha, SD, SDf, sig2eps, sig2rD, sig2bD, 
                        r.D, b.D, beta.Df, pr.Sigma0InvDf, pr.beta0Df,pr.A0.rD, pr.B0.rD, m, printSumms, t ){
  
  pD = length(b.D[[1]])
  for(j in 1:J){
    # j =1
    K=max(d.olddraws)
    if(printSumms==TRUE) cat('j=',j,': K.i=',K,' d.j=',d.olddraws[j]) 
    ## Tabulate the group sizes
    # j=1
    n.minus.j <- table(d.olddraws[-j])
    Kminus <- length(n.minus.j)
    h <- Kminus + m
    d.new.minus.j <- d.olddraws[-j]
    
    if( Kminus < K ){  
      ## The case where d_j != d_i for any other i:
      d.new.minus.j <- as.numeric(factor(d.olddraws[-j],labels=1:nlevels(factor(d.olddraws[-j]))))
      dj <- Kminus+1
      
      ## Need to draw new phi_c for m-1 possible new states:
      sig2rDnew <- as.list(rinvgamma(m-1, shape=pr.A0.rD , scale=pr.B0.rD ))
      sig2rDall <- c(sig2rD[-d.olddraws[j]], sig2rD[d.olddraws[j]], sig2rDnew )
      r.Dall <- c(r.D[-d.olddraws[j]], r.D[d.olddraws[j]], lapply(seq(1,(m-1)), FUN=function(x) rnorm(n, mean= 0, sd= sqrt(sig2rDnew[[x]]))))
      b.Dall <- c(b.D[-d.olddraws[j]], b.D[d.olddraws[j]], as.list(rnorm(m-1, mean= rep(0,pD), sd = sqrt(sig2bD))))
      beta.Dfall <- c(beta.Df[-d.olddraws[j]],beta.Df[d.olddraws[j]], lapply(seq(1,(m-1)), FUN=function(x) rmvnorm(1, mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf))))
    }  else {
      ## The case where d_j == d_i for some other i:
      
      ## Need to draw new phi_c for m possible new states:
      sig2rDnew <- as.list(rinvgamma(m, shape=pr.A0.rD , scale=pr.B0.rD ))
      sig2rDall <- c(sig2rD, sig2rDnew )
      r.Dall <- c(r.D, lapply(seq(1,m), FUN=function(x) rnorm(n, mean= 0, sd= sqrt(sig2rDnew[[x]]))))
      b.Dall <- c(b.D, as.list(rmvnorm(m, mean= rep(0,pD), sigma= diag(sig2bD,pD))))
      beta.Dfall <- c(beta.Df, lapply(seq(1,m), FUN=function(x) rmvnorm(1, mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf))))
    } 
    
    
    ## Log-likelihoods, and normalizing constant b_j for determining probability of membership for
    ##    d_j=d | all other d_i
    
    currmeans <- lapply(seq_along(1:h),FUN=function(d) matrix(fixed.mean[y.jlabs==j] + SDf%*%matrix(unlist(beta.Dfall[[d]]),pDF,1) + 
                                                                SD%*%matrix(unlist(b.Dall[[d]]),pD,1) + matrix(unlist(r.Dall[[d]],n,1))) )
    
    log.liks <- lapply(seq_along(1:h), FUN=function(d) sum(dnorm(Y[y.jlabs==j],mean=currmeans[[d]],sd=rep(sqrt(sig2eps[j]),n),log=T)))
    #     
    phi= -max(unlist(log.liks))
    e2phi <- lapply(X=1:h, FUN=function(d) exp(phi+log.liks[[d]]))
    weights <- c(as.vector((n.minus.j)),rep(alpha/m,m))/(J - 1 + alpha)
    
    
    ## Normalize using e^(c+phi) / sum_c e^(c+phi) :
    probs.j <- lapply(X=seq(1:h), FUN=function(d) weights[d]*e2phi[[d]]/sum(weights%*%unlist(e2phi)))
    if(printSumms==TRUE) print(unlist(probs.j))
    
    if(sum(is.na(unlist(probs.j)))+sum(is.infinite(unlist(probs.j)))>0){
      cat('NA in prob vec when d.draws=',d.olddraws,'j=',j,
          '\nlogliks=',unlist(log.liks),'\n probs=',unlist(probs.j))
    }
    d.new.j <- sample(1:h,1,prob=probs.j)
    if(printSumms==TRUE) cat(' d.new.j=', d.new.j,'\n')
    
    if( Kminus<K ){ 
      # i.e. we were looking at a singleton to begin with
      if(d.new.j > (Kminus+1)){
        #  i.e. stay a singleton but at a different state
        r.D[[d.olddraws[j]]] <- r.Dall[[d.new.j]]
        b.D[[d.olddraws[j]]] <- b.Dall[[d.new.j]]
        beta.Df[[d.olddraws[j]]] <- beta.Dfall[[d.new.j]]
        sig2rD[[d.olddraws[j]]] <- sig2rDall[[d.new.j]]
        # I don't think we need to reassign the domain a new number 
      }else{
        if(d.new.j <= Kminus){
          # i.e. join an existing domain
          r.D <- r.D[-d.olddraws[j]]
          b.D <- b.D[-d.olddraws[j]]
          beta.Df <- beta.Df[-d.olddraws[j]]
          sig2rD <- sig2rD[-d.olddraws[j]]
          d.olddraws[j] <- unique(d.olddraws[-j][which(d.new.minus.j==d.new.j)])
          d.olddraws <- as.numeric(factor(d.olddraws))
          # I don't think we need to re-order the 
        }
        # else, do nothing because d_j stayed a singleton with it's original state
      } 
    }
    
    if(Kminus==K){
      if(d.new.j > Kminus){
        # i.e. an outcome breaks out to it's own domain
        beta.Df[[length(beta.Df)+1]] <- beta.Dfall[[d.new.j]]
        r.D[[length(r.D)+1]] <- r.Dall[[d.new.j]]
        b.D <- c(b.D,b.Dall[[d.new.j]])
        sig2rD <- c(sig2rD,sig2rDall[[d.new.j]])
        d.olddraws[j] <- Kminus+1
      }else{
        # i.e. an outcome joins a new domain or stays in its original domain
        d.olddraws[j] <- d.new.j
      }
    }
    
    
    if(printSumms==TRUE) cat('d.draws=',d.olddraws,'\n\n')
    ### We can now update the d.olddraws and K
    K <- nlevels(as.factor(d.olddraws))
    
    log.liks <- probs.j <- currmeans <-  NULL
    beta.Dfall <- b.Dall <- r.Dall <- sig2rDall <- NULL
  }
  return(list(d.new=d.olddraws, betaDf.new=beta.Df, bD.new=b.D, rD.new=r.D, sig2rD.new=sig2rD))
}



###########################################################################################################
##### dissMatrix() is a function which takes the matrix of group assignments at each iteration and 
#####         calculates a dissimilarity matrix, which is a summary of the proportion of times the 
#####         observations are NOT grouped together (i.e. high probabilities = very dissimilar)
###########################################################################################################
dissMatrix <- function(foodat,printF=FALSE){
  nvars <- ncol(foodat)
  niter <- nrow(foodat)
  propMat <- matrix(0, nrow=nvars, ncol=nvars)
  for(i in 1:niter){
    if(printF==TRUE) cat('i=',i,'; ')
    distMat <- matrix(0, nrow=nvars, ncol=nvars)
    for(j in 1:nvars){
      distMat[j,which(foodat[i,]==foodat[i,j])] <- 1  
    }
    propMat <- propMat + distMat
  }
  propMat <- propMat/niter
  dissMat <- 1 - propMat
  return(dissMat)
}


###########################################################################################################
##### is.same() is a function which takes the final posterior clustering vector and a row of d.draws and 
#####    returns true if their levels "match". This replaces the identical() function, which will not work
#####    for the labelings since they may not be equal in value but are in fact equal in terms of domains.
###########################################################################################################

is.same <- function(rowd, dom.assns){
  nsame <- 0

  for(i in 1:length(rowd)){
    if(as.numeric(reorder(factor(rowd),dom.assns))[i]==dom.assns[i]){
      nsame <- nsame + 1
     }
  }
  return(nsame)
}

###########################################################################################################
##### summAll() is a function which takes the final posterior draws from the MCMC (after the posterior 
#####     clustering, which may have reduced our number of draws removed after the burnin) and produces key
#####     summary output.
###########################################################################################################


summAll <- function(mcmc.results, dom.assns, full=TRUE){
  varcols <- grep('sig2',colnames(mcmc.results))
  vardraws <- mcmc.results[,varcols]
  coefdraws <- matrix(mcmc.results[,-varcols], nrow=nrow(mcmc.results)) 
  if(full==TRUE){
    ### Summarize parameters with estimate, SE, and posterior intervals:
    resultsMat <- matrix(c(mean(coefdraws[,1]),sd(coefdraws[,1]),quantile(coefdraws[,1],probs=c(0.025,0.975))),nrow=1)
    ## Coefficient parameters
    if(ncol(coefdraws)>1){
      for(j in 2:ncol(coefdraws)){
        resultsMat <- cbind(resultsMat, cbind(mean(coefdraws[,j]),sd(coefdraws[,j]),quantile(coefdraws[,j],probs=c(0.025)),quantile(coefdraws[,j],probs=c(0.975))))
      }
    }
    ## Variance parameters
    for(j in 1:ncol(vardraws)){
      resultsMat <- cbind(resultsMat, cbind(median(vardraws[,j]),sd(vardraws[,j]),quantile(vardraws[,j],probs=c(0.025)),quantile(vardraws[,j],probs=c(0.975))))
    }
    ##### Name the matrix appropriately!
    colnames(resultsMat) <- c(paste(rep(colnames(mcmc.results)[-varcols],each=4),rep(c('est','se','low','up'),ncol(coefdraws)),sep='.'), 
                              paste(rep(colnames(mcmc.results)[varcols],each=4),rep(c('est','se','low','up'),length(varcols)),sep='.'))
  }else{
    ### Summarize parameters with estimate, SE: 
    resultsMat <- matrix(c(mean(coefdraws[,1]),sd(coefdraws[,1])),nrow=1)
    ## Coefficient parameters
    if(ncol(coefdraws)>1){
      for(j in 2:ncol(coefdraws)){
        resultsMat <- cbind(resultsMat, cbind(mean(coefdraws[,j]),sd(coefdraws[,j])))
      }
    }
    
    ## Variance parameters
    for(j in 1:ncol(vardraws)){
      resultsMat <- cbind(resultsMat, cbind(median(vardraws[,j]),sd(vardraws[,j])))
    }
    ##### Name the matrix appropriately!
    resultsMat <- cbind(resultsMat)
    colnames(resultsMat) <- c(paste(rep(colnames(mcmc.results)[-varcols],each=2),rep(c('est','se'),ncol(coefdraws)),sep='.'),
                              paste(rep(colnames(mcmc.results)[varcols],each=2),rep(c('est','se'),length(varcols)),sep='.')) 
    
  }
  
  
  return(resultsMat)
}

###########################################################################################################
##### fixDom.multouts() is a function which takes the starting values from the fitting for group assignments
#####         in the first stage of model-fitting, fixes the outcomes into domains based on LS clustering
#####         allowing for inference on the remaining model parameters.
###########################################################################################################
#####         The function takes the (1) the data vector, Y, stacked by outcomes, 
#####         (2) covariate matrices for each of the ones specified in Luo et al. 2014: Sf, SDf, SD, SO  
#####         (3/4) the outcome labels y.jlabs, and the observation labels, y.ilabs, for the Y vector 
#####         (5/6)  priors as a list, priors, and the initial values as a list, inits
#####         (7/8) min/max number of iterations to perform the mcmc sampler (9) prints iteration num if T
#####         (10) the FIXED domain assignments (11) the min ESS for bO.draws (all J of these)
#####         (12) a T/F to indicate whether to introduce a hierarchical prior on sig2rD
###########################################################################################################

fixDom.multouts <- function(Y, Sf, SDf, SD, SO,  y.jlabs, y.ilabs, priors, inits, min.niter, max.niter, printIter, dom.assns, bOEss ){
  #### y.jlabs = distinguish outcomes
  #### y.ilabs = distinguish subjects
  #### printIter = TRUE to show the iteration for the overall MCMC
  #### printSumms = TRUE to show the clustering process for each outcome at each iteration
  #### sim.seed = random seed used in simulating everything but the covariate matrices in the function multOutData() 
  
  # Define matrices/lists/vectors to store parameters, which are the size of the maximum number of iterations
  niter <- max.niter
  
  ##########################################
  ########## Define the following:
  nobs.vec <- table(y.jlabs)  # number of observations for each outcome (same if all have complete obs on all outcomes)
  pF <- ncol(Sf) 
  pDF<- ncol(SDf) 
  pD <- ncol(SD) 
  pO <- ncol(SO)
  J <- nlevels(factor(y.jlabs)) # number of outcomes
  n <- nlevels(factor(y.ilabs))# number of subjects/outcome
  y.mislabs <- apply(matrix(Y,ncol=1), 1, is.na) 
  nmis <- sum(y.mislabs)
  
  
  ##########################################
  ########## Matrices to Store Parameter Draws: 
  ## Domain assignments:
  d.draws <- dom.assns
  
  ## Determine Number of Domains for Stage 2:
  K <- nlevels(factor(d.draws))
  ## Determine whether we will have >1 domains (else we won't need several domain-specific parameters)
  if(K>1){
    multdomain <- TRUE
  }else{
    multdomain <- FALSE
  }
  ## Determine whether we have singletons in our data (also will not need several domain-specific parameters):
  Ktemp <- K - sum(table(d.draws)==1)
  if(Ktemp<K){
    ## identify numeric value of the singleton domain(s)
    singleton <- as.numeric(which(table(d.draws)==1))
  }else{
    singleton <- NULL
  }
  
  ## Non-domain specific parameters:
  betaf.draws <- matrix(NA, nrow=niter, ncol=pF)
  bO.draws <- matrix(NA, nrow=niter, ncol=pO*J)
  sig2bO.draws <- matrix(NA, nrow=niter, ncol=pO)
  r.draws <- matrix(NA, nrow=niter, ncol=n)
  sig2r.draws <- rep(NA, niter)
  sig2eps.draws <- matrix(NA, nrow=niter, ncol=J)
  sig2rD.draws <- matrix(NA, nrow=niter, ncol=K)
  sig2bD.draws <- matrix(NA, nrow=niter, ncol=pD)
  
  ## Domain-specific parameters: 
  betaDf.draws <- list()
  bD.draws <- list()
  rD.draws <- list()
  
  
  ##########################################
  ########## Initial Values (entered from inits=list()) :
  betaf.draws[1,] <- inits$betaF
  bO.draws[1,] <- inits$bO
  colnames(bO.draws) <- rep(1:J,each=pO)
  sig2bO.draws[1,] <- inits$sig2bO
  r.draws[1,] <-  inits$r
  sig2r.draws[1] <- inits$sig2r
  sig2eps.draws[1,] <- inits$sig2eps
  sig2rD.draws[1,] <- inits$sig2rD
  sig2bD.draws[1,] <-inits$sig2bD
  
  betaDf.draws[[1]] <- inits$betaDf
  bD.draws[[1]] <- inits$bD
  rD.draws <- inits$rD
  if(!is.null(singleton)){
    bD.draws[[1]][singleton] <- 0
  }

  ##########################################
  ########## Prior Values (entered from inits=list()) :
  pr.beta0f <- priors$beta0f
  pr.Sigma0Invf <- priors$SigInv.betaf
  pr.beta0Df <- priors$beta0Df
  pr.Sigma0InvDf <- priors$SigInv.betaDf 
  pr.A0.bO <- priors$sig2.bO[1]
  pr.B0.bO <- priors$sig2.bO[2]
  pr.A0.bD <- priors$sig2.bD[1]
  pr.B0.bD <- priors$sig2.bD[2] 
  pr.A0.eps <- priors$sig2.eps[1]
  pr.B0.eps <- priors$sig2.eps[2]
  pr.A0.r <- priors$sig2.r[1]
  pr.B0.r <- priors$sig2.r[2]
  pr.A0.rD <- priors$sig2.rD[1]
  pr.B0.rD <- priors$sig2.rD[2] 
  pr.alpha <- priors$alpha
  
  Rvec <- t(t(kronecker(c(rep(1,J)),r.draws[1,])))
 
  if( K > 1){
    Lambda <- cbind(1*(d.draws==1),model.matrix( ~ as.factor(d.draws)  ,contrasts='contr.sum')[,-1])
  }else{
    Lambda <- cbind(rep(1,J))
  }
  RDvec <- kronecker(Lambda,diag(1,n))%*%unlist(rD.draws)
  Zmat.D <- kronecker(Lambda, SD)
  Fmat.DF <- kronecker(Lambda,SDf) 
  Fmat.F <- kronecker(rep(1, J),Sf)  
  Zmat.O <- kronecker(diag(1,J), SO) 
  
  ### To Determine convergence
  MeHg.coef <- matrix(NA, nrow=niter,ncol=J)
  MeHg.coef[1,] <- betaf.draws[1,] + bO.draws[1,] + unlist(bD.draws[[1]][d.draws] )
  MeHg.ess <- 2000

  ## label outcomes according to domain assignment:
  d.labs <- rep(d.draws,nobs.vec)  
  ## number of outcomes in each domain:
  njD <- table(d.draws)  
  
  t <- 2
  while(t%in%2:niter){
    
    if(printIter==TRUE){
      if(t%%10000==0) if(printIter==TRUE) cat('\n t=',t,' K=',K,';\n')
    }
    
    # print(t)
    ########################
    ### ---> Draw missing values
    means.Y <- Fmat.F%*%(betaf.draws[t-1,]) + Fmat.DF%*%unlist(betaDf.draws[[t-1]]) + Zmat.D%*%unlist(bD.draws[[t-1]]) + Zmat.O%*%(bO.draws[t-1,]) + RDvec + Rvec
    sds.Y <- rep(sqrt(sig2eps.draws[t-1,]),each=n)
    Ymis <- rnorm(nmis, mean=means.Y[y.mislabs], sd=sds.Y[y.mislabs])
    Y[y.mislabs] <- Ymis    
    
    ################################################
    #### ---> Gibbs Steps for Domain-Specific Parameters:
    means.Df <- lapply(seq_along(1:K), FUN=function(k) pr.Sigma0InvDf%*%pr.beta0Df + 
                         t(kronecker(rep(1,njD[k]),SDf))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%
                         (as.matrix(Y[d.labs==k],ncol=1) - kronecker(rep(1,njD[k]),(Sf%*%betaf.draws[(t-1),])) - kronecker(rep(1,njD[k]),(SD%*%(bD.draws[[t-1]][[k]]))) 
                          - kronecker(diag(1,njD[k]),SO)%*%bO.draws[(t-1),(colnames(bO.draws) %in% y.jlabs[d.labs==k])] 
                          - as.matrix(kronecker(rep(1,njD[k]),r.draws[(t-1),]) - kronecker(rep(1,njD[k]),unlist(rD.draws[[k]])),ncol=1)) )
    var.betaDf <- lapply(seq_along(1:K), FUN=function(k) solve(pr.Sigma0InvDf + sum(1/1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])])*t(SDf)%*%SDf) )#t(kronecker(rep(1,njD[k]),SDf))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%kronecker(rep(1,njD[k]),SDf)) )
    betaDf.draws[[t]] <- lapply(X=seq_along(1:K), FUN = function(k) mvrnorm(1, mu=unlist(var.betaDf[[k]])%*%unlist(means.Df[[k]]) , Sigma =unlist(var.betaDf[[k]])))  
    
    if(multdomain){
      means.bD <- lapply(seq_along(1:K), FUN=function(k) t(kronecker(rep(1,njD[k]),SD))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%
                           (as.matrix(Y[d.labs==k],ncol=1) - kronecker(rep(1,njD[k]),(Sf%*%betaf.draws[(t-1),])) - kronecker(rep(1,njD[k]),(SDf%*%betaDf.draws[[t]][[k]])) 
                            - kronecker(diag(1,njD[k]),SO)%*%bO.draws[(t-1),(colnames(bO.draws) %in% y.jlabs[d.labs==k])] - as.matrix(kronecker(rep(1,njD[k]),r.draws[(t-1),]) 
                                                                                                                                      - kronecker(rep(1,njD[k]),(rD.draws[[k]])),ncol=1)) )
      var.bD <- lapply(seq_along(1:K), FUN=function(k) solve(diag(1/sig2bD.draws[(t-1),],pD) + t(kronecker(rep(1,njD[k]),SD))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%kronecker(rep(1,njD[k]),SD)) )
      bD.draws[[t]] <- lapply(X=seq_along(1:K), FUN = function(k) mvrnorm(1, mu=unlist(var.bD[[k]])%*%unlist(means.bD[[k]]) , Sigma =unlist(var.bD[[k]])))  
      if(!is.null(singleton)) bD.draws[[t]][singleton] <- 0
      
      bDd.sums <- apply(matrix(unlist(bD.draws[[t]]), nrow=pD),1,function(x) sum(x^2))
      if(!is.null(singleton)){
        sig2bD.draws[t,] <- rinvgamma(pD, shape=pr.A0.bD + (K-length(singleton))/2, scale=pr.B0.bD + .5*bDd.sums)
      }else{
        sig2bD.draws[t,] <- rinvgamma(pD, shape=pr.A0.bD + K/2, scale=pr.B0.bD + .5*bDd.sums)
      }
      
      var.rD <- lapply(seq_along(1:K), FUN=function(k) solve((1/sig2rD.draws[t-1,k]) + sum(1/sig2eps.draws[t-1,d.draws==k])) )
      separ.rD <- by((Y - Fmat.F%*%(betaf.draws[t-1,]) - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]])  - Zmat.O%*%(bO.draws[t-1,]) - Rvec),
                     INDICES=factor(y.jlabs),FUN=function(x) 1*x)
      means.rDbyout <- lapply(seq(1:J), FUN=function(j) (1/sig2eps.draws[t-1,j])*unlist(separ.rD[[j]]) )
      means.rDmat <- matrix(unlist(means.rDbyout),ncol=J,byrow=F)
      means.rD <- lapply(seq_along(1:K), FUN=function(k) 
        if( sum(d.draws==k)==1){
          means.rDmat[,k] 
        }else{
          rowSums(means.rDmat[,d.draws==k])
        }        )
      rD.draws <- lapply(seq_along(1:K), FUN=function(k) rnorm(n, mean=unlist(var.rD[[k]])%*%unlist(means.rD[[k]]) , sd=sqrt(unlist(var.rD[[k]])))) 
      if(!is.null(singleton)) for(i in 1:length(singleton)) rD.draws[[singleton[i]]] <- rep(0,n)
      RDvec <- kronecker(Lambda,diag(1,n))%*%unlist(rD.draws)
      
#       
#       if(!is.null(singleton)){
     sig2rDsums <- matrix(apply(matrix(unlist(rD.draws), nrow=K, byrow=T),1,function(x) sum(x^2)),nrow=1)
#         s2rDtemp <- c(3,3)
#         while(sum(s2rDtemp>2)>0){
#           s2rDtemp <- apply(sig2rDsums, 2, function(x) rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*x))
#         }
#       }
#       
      sig2rD.draws[t,] <- apply(sig2rDsums, 2, function(x) rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*x))
      if(sum(sig2rD.draws[t,]>3)>0) cat('\nt=',t,' sig2rD=',sig2rD.draws[t,])
    }else{
      ### Assign redundant parameters to be 0:
      rD.draws <- list(rep(0,n))
      bD.draws[[t]] <- 0
      sig2rD.draws[t,] <- rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD)
      sig2bD.draws[t,] <- rinvgamma(pD, shape=pr.A0.bD + K/2, scale=pr.B0.bD)
    }
        # print(rD.draws)
    
    ###################
    #### ---> Gibbs Steps for non-domain specific parameters:
    ## Changed to sum over the n-subjects
    VbetafInv <- solve(pr.Sigma0Invf + sum(1/sig2eps.draws[t-1,])*t(Sf)%*%Sf)
    means.betaf1 <- VbetafInv%*%(pr.Sigma0Invf%*%pr.beta0f +(rep(1/sig2eps.draws[t-1,],each=n)*rep(Sf,J))%*%(Y - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]]) - Zmat.O%*%(bO.draws[(t-1),]) - Rvec - RDvec))
    betaf.draws[t,] <- rnorm(1, mean=means.betaf1,  sd=sqrt(VbetafInv)) 
    
    VbOInv <- by(sig2eps.draws[t-1,], INDICES=factor(1:J), FUN=function(x) solve(diag(1/sig2bO.draws[t-1,],pO)+(1/x)*t(SO)%*%SO)) 
    means.bO <- lapply(seq_along(1:J), FUN=function(j) unlist(VbOInv[[j]])%*%(1/sig2eps.draws[t-1,j])*t(SO)%*%
                         (Y[y.jlabs==j]- Sf%*%(betaf.draws[t,]) - SDf%*%betaDf.draws[[t]][[d.draws[j]]] - SD%*%bD.draws[[t]][[d.draws[j]]] - r.draws[t-1,] - rD.draws[[d.draws[j]]]))                     
    bO.draws[t,] <- unlist(lapply(seq_along(1:J), FUN=function(j) mvrnorm(1, mu=unlist(means.bO[[j]]), Sigma=unlist(VbOInv[[j]])) ))
    
    
    ## EDIT: r.draws <- rnorm ( ... aggregate(mean.diff, by=list(outcomes=y.ilabs), FUN=function(x) sum(x/c(1:J)))[,'V1'] )
    VrInv <- solve(1/sig2r.draws[t-1] + sum(1/sig2eps.draws[t-1,]))
    mean.diffr <- (Y - Fmat.F%*%(betaf.draws[t,]) - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]])  - Zmat.O%*%(bO.draws[t,]) - RDvec)
    mean.diffrlist <- by(mean.diffr, INDICES=factor(y.jlabs), FUN=function(x) 1*x)
    means.rbyout <-  lapply(seq(1:J), FUN=function(j) (1/sig2eps.draws[t-1,j])*unlist(mean.diffrlist[[j]]))
    means.r <- apply(matrix(unlist(means.rbyout),nrow=J,byrow=T),2,function(x) VrInv*sum(x))
    r.draws[t,] <- rnorm(n, mean=means.r, sd=sqrt(VrInv))
    
    
    Rvec <- t(t(kronecker(c(rep(1,J)),r.draws[t,]))) 
    
    
    ### ---> Variance parameters:
    bOj.sums <- apply(matrix(bO.draws[t,], nrow=pO),1,function(x) sum(x^2))
    sig2bO.draws[t,] <- rinvgamma(pO, shape=pr.A0.bO+J/2, scale=pr.B0.bO+.5*bOj.sums)
    
    
    
    sig2r.draws[t] <- rinvgamma(1, shape=pr.A0.r + n/2, scale=pr.B0.r + .5*t(r.draws[t,])%*%r.draws[t,])
    
    
    mean.diff2 <- (Y - Fmat.F%*%(betaf.draws[t,]) - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]])  - Zmat.O%*%(bO.draws[t,]) - RDvec - Rvec)
    sig2eps.draws[t,] <- rinvgamma(J, shape=pr.A0.eps + n/2, scale=pr.B0.eps + .5*aggregate(mean.diff2, by=list(outcomes=y.jlabs), FUN=function(x) sum(x^2))[,'V1'])
    

    ###### After minimum number of iterations performed, check ESS for bO.draws
    ## niterT = the total number of iterations performed in this second stage (for recording purposes)
    MeHg.coef[t,] <- betaf.draws[t,] + bO.draws[t,] + unlist(bD.draws[[t]][d.draws])
    
    niterT <- t
    
    if(t>min.niter){
      if(t%%200==0){
        ess.test <- multiESS(MeHg.coef[1:t,])
        cat('\nt=',t,' MeHg.ess=',ess.test)
        # cat('\nESS>ess.test:',sum(ess.test>=bOESS))
        if(ess.test>=MeHg.ess){
          cat('\nConvergence based on',MeHg.ess,' is reached at t=',t)  
	        t <- niter           
        }
      }
    }
    t <- t + 1 
  }
  return( list( niterT.2=niterT, ess.test=ess.test, betaf.draws=betaf.draws,  sig2bO.draws=sig2bO.draws, sig2bD.draws=sig2bD.draws, 
                sig2r.draws=sig2r.draws, sig2rD.draws=sig2rD.draws, sig2eps.draws=sig2eps.draws, betaDf.draws=betaDf.draws ))
}


###########################################################################################################
##### renumber() is a function which takes the matrix of group assignments at each iteration and 
#####         relabels them in increasing order to make the groups more comparable across iterations
###########################################################################################################

renumber <- function(domdraws, Kpost=NULL){
  assn1st <- apply(domdraws, 1, function(x) seq_along(x)[!duplicated(x)] ) # index of 1st unique outcomes
  domorder <- apply(domdraws, 1, function(x) x[seq_along(x)[!duplicated(x)]] ) # value of unique outcomes in order
  relabds <- matrix(NA, nrow=nrow(domdraws), ncol=ncol(domdraws))
  J <- ncol(domdraws)
  
  for(i in 1:nrow(domdraws)){
    # i = 1
    if(is.null(Kpost)){
      Kpost <- length(domorder[[i]])
      relabds[i, ] <- factor(domdraws[i,], levels = domorder[[i]], labels = 1:Kpost)
      Kpost <- NULL
    } else {
      # Kpost <- length(domorder[,i])
      relabds[i, ] <- factor(domdraws[i,], levels = domorder[,i], labels = 1:Kpost)
    }
  }
  
  p=1
  relabs2gps <- relabds
  patterns <- list()
  groups <- list()
  k=1
  while(nrow(relabs2gps)>0){
    p1 <- relabs2gps[1,] #assn1st[,p]
    patterns[[k]] <- p1
    relabs2gps <- matrix(relabs2gps[-1,],ncol=J)
    gp <- which(apply(relabs2gps,1,function(x) sum(x==p1)==J))
    groups[[k]] <- c(p,gp+p)
    if(length(gp)>0)  relabs2gps <- matrix(relabs2gps[-gp,],ncol=J)
    p <- p + length(gp) + 1
    k <- k + 1
  }
  
  return(list(newds=relabds,gpindices=groups,gppatterns=patterns))
}

###########################################################################################################
##### dahlLSclust() is a function which takes the matrix of group assignments at each iteration and 
#####         computes least squared error between the association matrix at each iteration, t, and
#####         the estimated similarity matrix, which is a summary of the proportion of times the 
#####         observations are grouped together (i.e. 1 - dissimilarity matrix)
###########################################################################################################

dahlLSclust <- function(d.draws, S.hat ){
  nouts <- ncol(d.draws)
  niter <- nrow(d.draws)
  sums.t <- rep(NA, niter)
  for(t in 1:niter){
    ## Create association matrix at each iteration:
    assocMat <- matrix(0, nrow=nouts, ncol=nouts)
    for(j in 1:nouts){
      assocMat[j,which(d.draws[t,]==d.draws[t,j])] <- 1  
    }
    sumij <- 0
    for(i in 1:nouts){
      sumij <- sumij + sum((assocMat[i,] - S.hat[i,])^2  )
    }
    sums.t[t] <- sumij
  } 
  return(sums.t)
}

