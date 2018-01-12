## Parameter values similar to our results

letter = 'a'

J = 20
n = 500
pDF = 6

domAssn <- c(1,1,1,2,2,2,2,2,1,3,3,4,4,5,5,6,6,7,7,6)
D <- length(unique(domAssn))

betaf <-  0.05

## Domain specific Fixed effects: 
# large discrepancies
ncovs = pDF
betaDf = array(0,c(ncovs,D))
betaDf[1,] <- c( -0.02,-0.10, 0.10, 0.15, 0.25,-0.15,-0.05 ) # child's sex
betaDf[2,] <- c(  0.10, 0.01, 0.05, 0.02, 0.05, 0.10,-0.15 ) # mother's age
betaDf[3,] <- c(  0.25, 0.10, 0.10, 0.15, 0.20, 0.15,-0.10 ) # HOME
betaDf[4,] <- c(  0.15, 0.02, 0.05, 0.10, 0.02, 0.15,-0.05 ) # KBIT
betaDf[5,] <- c(  0.15, 0.05, 0.10, 0.05, 0.08, 0.02,-0.05 ) # Hollingshead
betaDf[6,] <- c( -0.02,-0.15,-0.05, 0.15, 0.20, 0.25, 0.05 ) # child's age

##  Domain specific Random effects:
bD1 <- -.02
bD2 <- 0.05
bD3 <- 0.03
bD4 <- -.02
bD5 <- -.01
bD6 <- -.02
bD7 <- -.01
bD <- c(bD1, bD2, bD3, bD4, bD5, bD6, bD7)
var(bD); sd(bD); sum(bD)

s2rD = c(.3,.5,.5,.6,.4,.3,.7) 

## Subject-specific Random effects: 
s2r=.4^2
s2bO = .02^2
