# given adjacency matrix and error variances for a causal linear model
# computes the covariance matrix
Covariance <- function(B,eps){
  n <- nrow(B)
  Id <- diag(n)
  Cov <- solve(Id - B)%*% diag(eps) %*% t(solve(Id - B))
  return(Cov)
}

# computed the asymptotic variance for the 
# OLS regression coefficient corresponding to X in the regression Y~X+W
aVar <- function(x,y,W,CovG){
  if(length(W)==0){
    sigma_xx <- CovG[x,x]
    sigma_yy <- CovG[y,y]
    sigma_yx <- CovG[y,x]
    
    aVarW <- (sigma_yy - sigma_yx %*% solve(sigma_xx) %*% sigma_yx) / sigma_xx
    return(aVarW)
  }
  else{
    ## struggling with the transpose (probably as Cov(x,W) is of numeric class)
    ## as.matrix seems to do an automatic transpose (seems a bit sketchy though)
    Sigma_xw <- CovG[x,W]
    conVarX <- CovG[x,x] -Sigma_xw %*% solve(CovG[W,W]) %*% as.matrix(Sigma_xw)
    
    
    ## same isse as above with Cov(y,W2)
    W2 <- append(x,W)
    Sigma_yw2 <- CovG[y,W2]
    conVarY <- CovG[y,y] -Sigma_yw2 %*% solve(CovG[W2,W2]) %*% as.matrix(Sigma_yw2)
    aVarW <- conVarY / conVarX
    return(aVarW)
  }
}
