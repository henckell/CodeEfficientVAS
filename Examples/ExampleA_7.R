library(xtable)

setwd("~/Dokumente/Papers/Rcode/CodeEfficientVAS/Examples")

#### getting necessary functions

source("~/Dokumente/Papers/Rcode/CodeEfficientVAS/Examples/helpFunctions.R")

RNGkind(sample.kind = "Rounding")

# setting sample size
n <- 1000

# Setting up adjacency matrix
# Node order: a1,a2,b2,b1,x,y 
B <- matrix(c(0,0,0,0,0,0,
              1,0,0,0,0,0,
              0,0,0,0,0,0,
              0,0,1,0,0,0,
              1,0,0,1,0,0,
              0,1,20,0,1,0),nrow=6)

# Computing scaled asymptotic varainces
Cov <- Covariance(t(B),rep(1/3,6))

avar <- c()
avar[1] <- aVar(5,6,c(2,4),Cov)/n
avar[2] <- aVar(2,6,c(4,5),Cov)/n
avar[3] <- aVar(4,6,c(2,5),Cov)/n

# Simulating data to compute empirical variances
set.seed(60)

k <- 1000
beta <- matrix(numeric(k*3),nrow=3)
for(i in 1:k){
n <- 1000
a1 <- runif(n,-1,1)
a2 <- runif(n,-1,1) +a1
b2 <- runif(n,-1,1) 
b1 <- runif(n,-1,1) +b2
x <- runif(n,-1,1) + a1 + b1
y <- runif(n,-1,1) + a2 + 20* b2 + x

fit <- lm(y~x+a2+b1-1)
beta[,i] <- fit$coefficients
}

# empirical mean and variance
apply(beta,1,stats::var)
apply(beta,1,mean)

# produce plots for Example A.7, with k=1000 OLS regression
# residuals vs fitted values
pdf("residVSfitted.pdf")
plot(fitted(fit), resid(fit), 
          ylab="Residuals", xlab="Fitted values", 
          main="Residual plot") 
abline(0, 0)  
dev.off()

# residuals vs X
pdf("residVSx.pdf")
plot(x, resid(fit), 
     ylab="Residuals", xlab="X", 
     main="Residual plot") 
abline(0, 0) 
dev.off()


# produce table for Example A.7 
res <- rbind(apply(beta,1,mean),apply(beta,1,stats::var),avar)
xtable(t(res))
