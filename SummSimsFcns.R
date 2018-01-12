library(xtable)


summSimsMO <- function(fullres, domAssn,  True, vardigits=4, J=20, summParDtrue=TRUE){
  
  nsim=nrow(fullres)
  algMat <- matrix(NA, nrow=7, ncol=1)
  rownames(algMat) <- c('Stage 1: Iterations','Stage 1: ESS','Stage 2: Iterations','Stage 2: ESS',
                        'Perc. Times Domains=D','No. Times Domains<D','No. Times Domains>D')
  colnames(algMat) <- c('Simulation Results')
  
  algMat['Stage 1: Iterations',] <- round(median(fullres[,'niterT.1']))
  algMat['Stage 1: ESS',] <- round(mean(fullres[,'ess.test1']))
  algMat['Stage 2: Iterations',] <- round(median(fullres[,'niterT.2']))
  algMat['Stage 2: ESS',] <- round(mean(fullres[,'ess.test2']))
  algMat['Perc. Times Domains=D',] <- 100*sum(fullres[,'Kpost']==max(domAssn))/nrow(fullres)
  algMat['No. Times Domains<D',] <- sum(fullres[,'Kpost']<max(domAssn))
  algMat['No. Times Domains>D',] <- sum(fullres[,'Kpost']>max(domAssn))
  
  clustMat <- matrix(NA, nrow=length(unique(fullres$Kpost)), ncol=4)
  colnames(clustMat) <- c('No. Domains','No. Simulations','No. Patterns','No. Correct - True D')
  
  clustMat[,'No. Domains'] <- paste('D =',sort(unique(fullres$Kpost)))
  clustMat[,'No. Simulations'] <- table(fullres$Kpost)
  clustMat[,'No. Patterns'] <- unlist(lapply(seq_along(table(fullres$Kpost)), function(x) nrow(unique(fullres[fullres$Kpost==names(table(fullres$Kpost))[x],2:(J+1)]))))
  clustMat[names(table(fullres$Kpost))==True$D,'No. Correct - True D'] <- sum(apply(fullres[,2:(J+1)], 1, function(x) identical(domAssn, as.numeric(x))))
  clustMat[names(table(fullres$Kpost))!=True$D,'No. Correct - True D'] <- rep('',sum(names(table(fullres$Kpost))!=True$D))
  
  
  all4Mat <- matrix(NA, nrow=3, ncol=5)
  rownames(all4Mat) <- c("$\\beta_f$", "$\\sigma_{bO}^2$", "$\\sigma_{bD}^2$")
  # colnames(all4Mat) <- c('True','Estimate','SE', '2.5%','97.5%')
  colnames(all4Mat) <- c('True','Estimate', '2.5%','97.5%','MSE')
  
  mse <- function(x, true, nsims){
    sum((x-true)^2)/nsims
  }
  if(is.na(True$s2bD)){
    True$s2bD = 0
  }
  # all4Mat[1,] <- c(True$betaf, round(apply(fullres[,grep('betaf',colnames(fullres))], 2, mean),digits=vardigits))
  # all4Mat[2,] <- c(round(True$s2bO,digits=vardigits), round(apply(fullres[,grep('sig2bO',colnames(fullres))], 2, mean),digits=vardigits))
  # all4Mat[3,] <- c(round(True$s2bD,digits=vardigits), round(apply(fullres[,grep('sig2bD',colnames(fullres))], 2, mean),digits=vardigits))
  if(summParDtrue){
    fullres <- fullres[fullres$Kpost==True$D,]
  }
  all4Mat[1,] <- c(True$betaf, round(c(mean(fullres[,'betaf.est']), mean(fullres[,'betaf.low']), mean(fullres[,'betaf.up'])
                                       ,mse(fullres[,'betaf.est'],true=True$betaf,nsims)),digits=vardigits) )
  all4Mat[2,] <- c(round(True$s2bO,digits=vardigits), round(c(mean(fullres[,'sig2bO.est']), mean(fullres[,'sig2bO.low']), mean(fullres[,'sig2bO.up'])
                                                              ,mse(fullres[,'sig2bO.est'],true=True$s2bO,nsims)),digits=vardigits) )
  all4Mat[3,] <- c(round(True$s2bD,digits=vardigits), round(c(mean(fullres[,'sig2bD.est']), mean(fullres[,'sig2bD.low']), mean(fullres[,'sig2bD.up'])
                                                              ,mse(fullres[,'sig2bD.est'],true=True$s2bD,nsims)),digits=vardigits) )
  
  domVarsMat <- matrix(NA, nrow=(6*True$D + True$D), ncol=3)
  colnames(domVarsMat) <- c('True','Estimate','SE')
  rownames(domVarsMat) <- c(paste0(paste0('$\\beta_{Df',seq(1,6)),',',rep(1:True$D,each=6),'}$'), paste0('$\\sigma^2_{rD,',1:True$D,'}$'))
  
  ## Pick only the covariate estimates/SEs for the matching D 
  #  This does not exclude simulations which did not pick D in the posterior
  betaDfcols.est <- colnames(fullres)[grep('betaDf', colnames(fullres))][grep('.est',colnames(fullres)[grep('betaDf', colnames(fullres))])][1:(6*True$D)]
  betaDfcols.se <- colnames(fullres)[grep('betaDf', colnames(fullres))][grep('.se',colnames(fullres)[grep('betaDf', colnames(fullres))])][1:(6*True$D)]
  s2rDcols.est <- colnames(fullres)[grep('sig2rD', colnames(fullres))][grep('.est',colnames(fullres)[grep('sig2rD', colnames(fullres))])][1:(True$D)]
  s2rDcols.se <- colnames(fullres)[grep('sig2rD', colnames(fullres))][grep('.se',colnames(fullres)[grep('sig2rD', colnames(fullres))])][1:(True$D)]
  
  domVarsMat[1:(6*True$D),'True'] <- c(True$betaDf)
  domVarsMat[1:(6*True$D),'Estimate'] <- round(apply(fullres[,betaDfcols.est], 2, function(x) mean(na.omit(x))) ,3)
  domVarsMat[1:(6*True$D),'SE'] <- round(apply(fullres[,betaDfcols.se], 2, function(x) mean(na.omit(x))) ,3)
  
  domVarsMat[(6*True$D+1):nrow(domVarsMat),'True'] <- True$s2rD
  domVarsMat[(6*True$D+1):nrow(domVarsMat),'Estimate'] <- round(apply(as.matrix(fullres[,s2rDcols.est]), 2, function(x) median(na.omit(x))) ,digits=vardigits)
  domVarsMat[(6*True$D+1):nrow(domVarsMat),'SE'] <- round(apply(as.matrix(fullres[,s2rDcols.se]), 2, function(x) median(na.omit(x))) ,digits=vardigits)
  
  epsMat <- matrix(NA, nrow=J, ncol=2)
  colnames(epsMat) <- c('Estimate','SE')
  rownames(epsMat) <- paste0('$\\sigma^2_{\\epsilon,',1:J,'}$')
    
  s2epscols.est <- colnames(fullres)[grep('sig2eps', colnames(fullres))][grep('.est',colnames(fullres)[grep('sig2eps', colnames(fullres))])]
  s2epscols.se <- colnames(fullres)[grep('sig2eps', colnames(fullres))][grep('.se',colnames(fullres)[grep('sig2eps', colnames(fullres))])]
  
  epsMat[, 'Estimate'] <- round(apply(fullres[,s2epscols.est], 2, function(x) median(na.omit(x))) ,3)
  epsMat[, 'SE'] <- round(apply(fullres[,s2epscols.se], 2, function(x) median(na.omit(x))) ,3)
  
  return(list(algMat=algMat, domVarsMat=domVarsMat, all4Mat=all4Mat, epsMat=epsMat, clustMat=clustMat))
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
