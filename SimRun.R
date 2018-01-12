source('clustMultOutwImp_forsims.R')
source('SimulateData.R')

args <- commandArgs(trailingOnly=TRUE)
njob <- as.numeric(args[1])
set.param.script <- args[2]
stored.seeds <- args[3]

cat('Job number: ', njob, '\n')
cat('Set param script: ', set.param.script, '\n')

source(set.param.script)
seedfile <- read.table(stored.seeds)

simData <- multOutData(betaf=betaf, bD=bD, betaDf=betaDf, s2bO=s2bO, s2rD=s2rD, s2r=s2r, n=n, letter=letter,
                       domAssn=domAssn, seed=njob, seedfile=seedfile, Dstart=10)
print('Starting Algorithm:\n\n')
fullresults <- cluster.multouts(Y=simData$Y, Sf=simData$covs$Sf, SDf=simData$covs$SDf, SD=simData$covs$SD, SO=simData$covs$Sf,  
                                y.jlabs=simData$y.jlabs, y.ilabs=simData$y.ilabs, priors=simData$priors, inits=simData$inits, 
                                domAssn=domAssn, m=10, tune=c(5,5), nburn=500, min.niter=5000, min.niter2=5000, max.niter=30000, 
                                printIter=TRUE, printSumms=FALSE, sim.seed=njob )

resultsclust <- data.frame(njob, fullresults )
print(simData$fileName)
write.table(resultsclust, file=simData$fileName)