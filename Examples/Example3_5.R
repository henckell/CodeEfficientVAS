library("xtable")

#### getting necessary functions

source("~/CodeEfficientVAS/Examples/helpFunctions.R")

set.seed(1)

## all edge weights one
SavedResults<-matrix(numeric(64), ncol=8)
SavedCoeff<-matrix(numeric(48), ncol=8)
SavedEps<-matrix(numeric(48), ncol=8)


for(i in 1:8){
  # random edge coefficients and error variances
  Coeff <- rnorm(6)
  Eps <- runif(6)

  # setting upt adjacency matrix
  B <- t(matrix(c(0,0,0,0,0,0,
                  0,0,0,0,0,0,
                  0,0,0,0,0,0,
                  Coeff[1],Coeff[2],0,0,0,0,
                  0,Coeff[3],Coeff[4],Coeff[5],0,0,
                  0,0,0,Coeff[6],0,0), ncol=6))

  # computing corresponding covariance matrix
  Cov <- Covariance(B,Eps)

  # computing asymptotic variances
  Results <- matrix(numeric(8), ncol=1)

  Results[1,1] <- aVar(c(4),c(5),c(1,2),Cov) 
  Results[2,1] <- aVar(c(4),c(5),c(1,2,3),Cov) 
  Results[3,1] <- aVar(c(4),c(5),c(2),Cov) 
  Results[4,1] <- aVar(c(4),c(5),c(2,3),Cov)
  Results[5,1] <- aVar(c(4),c(5),c(1,2,6),Cov) 
  Results[6,1] <- aVar(c(4),c(5),c(1,2,3,6),Cov)
  Results[7,1] <- aVar(c(4),c(5),c(2,6),Cov) 
  Results[8,1] <- aVar(c(4),c(5),c(2,3,6),Cov)

  # saving results
  SavedResults[,i]<-Results
  SavedCoeff[,i]<-Coeff
  SavedEps[,i]<-Eps
}

# print(SavedResults)

# table of variances for Example 3.5
xtable(SavedResults[,1:6])
