## function that computes partial causal effect from B matrix x, y and the conditioning set z 
causalEffect <- function (B, y, x, z) {
  # requires topolically ordered B
  p <- ncol(B)
  vec <- matrix(0, p, 1)
  vec[x] <- 1
  if (y - x > 1) {
    for (i in (x + 1):y){
      if(sum(i==z)==0){vec[i] <- B[i, ] %*% vec} 
      else{vec[i]=0}
    }
    return(vec[y])
  }
  else {
    return(B[y, x])
  }
}

## function that computes the joint total effect of x on y given B
causalEffectLeo <- function (B, y, x) {
  # requires topolically ordered B
  k <- length(x)
  effect <- numeric(k)
  for(i in 1:k){
    effect[i] <- causalEffect(B,y,x[i],x[-i])
  }
  return(effect)
}

## function that computes possiple parents
possparents <- function(x, amat){
  posspaX <- which(amat[x,]==1)
  return(posspaX)
}

## function that computes parents 
parents <- function(x, amat){
  paX <- setdiff(possparents(x,amat),which((amat[,x] + t(amat[x,]))==2))
  return(paX)
}

## redone samp function that fixes issue when pluggin in a vector of length 1
resamp <- function(x,...){if(length(x)==1) x else sample(x,...)} 


## function that computes O given x,y and the dag.amat matrix of a graph, works for amenable CPDAGs as well
ComputeO <- function(x,y,amat){
  
  ## obtaining forb(x,y,G)
  Forbidden <- pcalg:::forbiddenNodes(amat,x,y)
  
  ## obtaining an(y)
  AnY <- pcalg:::possibleAnProper(amat,y)
  
  ## intersecting the two to get causal node set (not including x as intended)
  CausalNodes <- intersect(Forbidden, AnY)
  
  ###### constructing best.set
  ## initiating best set
  best.set <-  c()
  
  ## adding parents of all causal nodes iteratively, can probably be done more elegantly
  for(i in CausalNodes){ 
    pa_i <- parents(i,amat)
    best.set <- union(best.set, pa_i)
    ## changed to union
  }
  
  ## removing forbidden nodes that are parents
  best.set <- setdiff(best.set, Forbidden)
  
  ## removing x as oddly it is not classified as a forbidden node
  best.set <- setdiff(best.set, x)
  
  return(best.set)
}


## help lm function
lmLeo <- function(data,x,y,set,X,Y){
  k<- length(x)
  if(length(set)==0 & k==1){
    fit <- lm(Y~X-1)
  }  else{
    Data<-data[,union(x,set)]
    fit <- lm(Y~.-1,data=Data)
  }
  Beta <- fit$coefficients[1:k]
  return(Beta)
}

## function that given B and epsilon computes covariance matrix for the causal linear models
Covariance <- function(B,eps){
  n <- nrow(B)
  Id <- diag(n)
  Cov <- solve(Id - B)%*% diag(eps) %*% t(solve(Id - B))
  return(Cov)
}