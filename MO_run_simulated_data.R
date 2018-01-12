library(plyr)
library(dplyr)

source('clustMultOutwImp_forsims.R')
source('SimulateData.R')

source('MO_simulated_data.R')
njob = sample(1:100,1) 
seedfile <- read.table('Seedfile.txt')
  
simData <- multOutData(betaf=betaf, bD=bD, betaDf=betaDf, s2bO=s2bO, s2rD=s2rD, s2r=s2r, n=n, letter=letter,
                       domAssn=domAssn, seed=njob, seedfile=seedfile, Dstart=10)

print('Starting MCMC:\n\n')
fullresults <- cluster.multouts(Y=simData$Y, Sf=simData$covs$Sf, SDf=simData$covs$SDf, SD=simData$covs$SD, SO=simData$covs$Sf,  
                                y.jlabs=simData$y.jlabs, y.ilabs=simData$y.ilabs, priors=simData$priors, inits=simData$inits, 
                                domAssn=domAssn, m=10, tune=c(5,5), nburn=500, min.niter=5000, min.niter2=5000, max.niter=30000, 
                                printIter=TRUE, printSumms=FALSE, sim.seed=njob )

## Domain Assignments:
domRes <- data.frame('True'=domAssn, 'Est'=t(fullresults[1,paste('d',1:length(domAssn),sep='.')])) 
domRes <- domRes %>% tbl_df() %>%
  mutate(Match=True==X2.5., Est=X2.5.) %>%
  select(True, Est, Match)
domRes

##### Fixed Effects Parameter Estimates:
# overall fixed effect
betafRes <- matrix(fullresults[1,grep('betaf',colnames(fullresults))], ncol=4, byrow = T)
colnames(betafRes) <- c('Est', 'SD','2.5%','97.5%')
betafRes <- cbind(True=betaf, betafRes)
betafRes
# domain-specific fixed effects
betaDFnames <- unique(apply(matrix(colnames(fullresults)[grep('betaDf',colnames(fullresults))],ncol=1), 1, function(x) substr(x,1,9) ))
betaDfRes <- matrix(fullresults[1,grep('betaDf',colnames(fullresults))], ncol=2, byrow = T)
rownames(betaDfRes) <- betaDFnames
colnames(betaDfRes) <- c('Est', 'SD')
betaDfRes <- cbind('True'=c(betaDf), betaDfRes)
betaDfRes

##### Random Effects/Variance Components Estimates:
# variance - domain-specific random subject effect
sig2rDRes <- matrix(fullresults[1,grep('sig2rD',colnames(fullresults))], ncol=2, byrow = T)
colnames(sig2rDRes) <- c('Est', 'SD')
sig2rDRes <- cbind(True=s2rD, sig2rDRes)
sig2rDRes

# variance - outcome-specific random effects
sig2bORes <- matrix(fullresults[1,grep('sig2bO',colnames(fullresults))], ncol=4, byrow = T)
colnames(sig2bORes) <- c('Est', 'SD','2.5%','97.5%')
sig2bORes <- cbind(True=s2bO, sig2bORes)
sig2bORes

# variance - domain-specific random effects
sig2bDRes <- matrix(fullresults[1,grep('sig2bD',colnames(fullresults))], ncol=4, byrow = T)
colnames(sig2bDRes) <- c('Est', 'SD','2.5%','97.5%')
sig2bDRes <- cbind(True=var(bD), sig2bDRes)
sig2bDRes

# variance - residual variance
sig2epsRes <- matrix(fullresults[1,grep('sig2eps',colnames(fullresults))], ncol=2, byrow = T)
rownames(sig2epsRes) <- paste('sig2eps',1:J,sep='.')
colnames(sig2epsRes) <- c('Est', 'SD')
sig2epsRes
