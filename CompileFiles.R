

compileResults <- function(nsims, J, filedir, simname ){
  
  #########################
  ### Create vector with all final number of domains:
  Kpost <- rep(NA, nsims)
  succ.sims <- rep(NA,nsims)
  for(i in 1:nsims){
    sim.i <- paste(simname,i,'.txt',sep='')
    if(file.exists(paste0('V:/Documents/THESIS/SIMULATIONS-MO/',filedir,sim.i)) ){
      readtab <- read.table(file=paste('V:/Documents/THESIS/SIMULATIONS-MO/',filedir,sim.i,sep=''), header = T)
      Kpost[i] <- readtab[,'Kpost']
      succ.sims[i] <- 1
    }else{
      succ.sims[i] <- Kpost[i] <- 0
    }
  }
  
  #########################
  ### Create matrix which can hold all parameters in the model with largest number of domains:
  fixed.dim <- 1 + J + 8 + 12  # njob + labels + sim.stats + (3params*4stats)
  tabResults <- matrix(NA, nrow=1, ncol = (fixed.dim + 2*(1 + J + 7*max(Kpost))))
  # --> column names for this matrix will be the same as the largest model:
  colnames(tabResults) <- colnames(read.table(file=paste0('V:/Documents/THESIS/SIMULATIONS-MO/',filedir,paste(simname,which.max(Kpost),'.txt',sep='')), header = T))
  
  #########################
  ### Read in tables and adjust the betaDF and sig2rD entries in Kpost < max(Kpost)
  for(i in 2:(nsims+1)){
    # cats('\ni=',i-1,'\nKpost=',Kpost[i-1])
    destfile=paste0('V:/Documents/THESIS/SIMULATIONS-MO/',filedir,paste(simname,i-1,'.txt',sep=''))
    if(file.exists(destfile)){
     readtab <- read.table(file=destfile, header = T)
      if(Kpost[(i-1)]==max(Kpost)){
        tabResults <- rbind(tabResults, as.numeric(readtab[1,]))
      }else{
        firstFixed <- readtab[,c(1:fixed.dim)]
        tabMinFixed <- readtab[,-c(1:fixed.dim)]
        fillmiss.betaDF <-  cbind((tabMinFixed[1,-grep('sig2',colnames(tabMinFixed))]), matrix(NA,nrow=1, ncol=(max(Kpost)-Kpost[i-1])*12))
        
        tabMinbetaDF <- tabMinFixed[,-grep('betaD',colnames(tabMinFixed))]
        fillmiss.sig2rD <- cbind((tabMinbetaDF[1,grep('sig2rD',colnames(tabMinbetaDF))]), matrix(NA,nrow=1, ncol=(max(Kpost)-Kpost[i-1])*2))
        tabMinS2rD <- tabMinbetaDF[,-grep('sig2rD',colnames(tabMinbetaDF))]
        tempTab <- cbind(firstFixed,fillmiss.betaDF, tabMinS2rD, fillmiss.sig2rD)
        colnames(tempTab) <- colnames(tabResults)
        tabResults<- rbind(tabResults,tempTab)
      }
    }
  }
  return(tabResults)
}


## compileResults(nsims, J, filedir, simname )
fullres <- compileResults(nsims = 100, J = 20, filedir = 'INSERT_FOLDER_NAME_WHERE_RESULTS_STORED', simname = 'J.20.n.500.D.7.a.ß')
write.table(fullres[-1,], file=newfilename,  col.names = T)

source('SummSimsFcns.R')
source('SetParam.R')
## List of true parameters
True1 <- list(betaf = betaf, s2bD=var(bD), s2rD=s2rD, s2r=s2r, s2bO=s2bO, betaDf=betaDf, D=D, nsims=100)
## Function to summarize simulations
summTabs1 <- summSimsMO(fullres = fullres1, domAssn = domAssn, True = True1, J=J)

## Print tables of algorithm performance and parameter estimation
print(xtable(summTabs1$algMat, digits=0, caption = 'Summary of overall algorithm performance.'))
print(xtable(summTabs1$clustMat, digits=0, caption = 'Summary of posterior clusters found across all simulations.'), include.rownames = F)
print(xtable(summTabs1$all4Mat, digits=4, caption = 'Summary of overall mercury effect and associated random effect variance parameter estimates.'),sanitize.rownames.function = identity)
print(xtable(summTabs1$epsMat, digits=3, caption = 'Summary of residual variance parameter estimates.'),sanitize.rownames.function = identity)
