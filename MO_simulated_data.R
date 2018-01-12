## Parameter values similar to our results
letter = 'a'

J = 20   # outcomes
n = 500  # subjects
pDF = 6  # domain-specific covariates treated as fixed effects

domAssn <- c(1,1,1,2,2,2,2,2,1,3,3,4,4,5,5,6,6,7,7,6)  # domain assigment
D <- length(unique(domAssn))                           # unique domains

betaf <-  0.05   # overall mercury effect

## Domain specific Fixed effects: 
betaDf = array(0,c(pDF,D))
betaDf[1,] <- c( -0.02,-0.10, 0.10, 0.15, 0.25,-0.15,-0.05 ) # child's sex
betaDf[2,] <- c(  0.10, 0.01, 0.05, 0.02, 0.05, 0.10,-0.15 ) # mother's age
betaDf[3,] <- c(  0.25, 0.10, 0.10, 0.15, 0.20, 0.15,-0.10 ) # HOME
betaDf[4,] <- c(  0.15, 0.02, 0.05, 0.10, 0.02, 0.15,-0.05 ) # KBIT
betaDf[5,] <- c(  0.15, 0.05, 0.10, 0.05, 0.08, 0.02,-0.05 ) # Hollingshead
betaDf[6,] <- c( -0.02,-0.15,-0.05, 0.15, 0.20, 0.25, 0.05 ) # child's age

##  Domain-specific random effects for mercury exposure:
bD1 <- -.02
bD2 <- 0.05
bD3 <- 0.03
bD4 <- -.02
bD5 <- -.01
bD6 <- -.02
bD7 <- -.01
bD <- c(bD1, bD2, bD3, bD4, bD5, bD6, bD7)
var(bD); sd(bD); sum(bD)

##  Variance of domain-specific random subject effects:
s2rD = c(.3,.5,.5,.6,.4,.3,.7) 

## Variance of subject-specific random effects: 
s2r=.4^2

## Variance of outcome-specific random effects for mercury exposure: 
s2bO = .02^2


