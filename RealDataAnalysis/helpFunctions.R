#function to compute estimated variances given for our two estimates of interest
# for any given data, adjacency matrix and pair (X,Y)
CIeval <- function(data,B,x,y){
   
  O <- ComputeO(x,y,B)
  paX <- parents(x,B)
  
  Beta <- numeric(3)
  CI <- list()
  SE <- list()
  
  data.frame <- as.data.frame(data)
  
  X<-data.frame[,x]
  Y<-data.frame[,y]
  
  ## fit with empty conditioning set
  fit.Triv <- lm(Y~X-1)
  Beta[1] <- fit.Triv$coefficients[1]
  CI[[1]] <- confint(fit.Triv,1)
  SE[[1]] <- vcov(fit.Triv)[1,1]
  
  ## DAG
  if(length(paX)==0){
    fit.paX <- lm(Y~X-1)
  }  else{
    Data.paX<-data.frame[,union(x,paX)]
    fit.paX <- lm(Y~.-1,data=Data.paX)
  }
  Beta[2] <- fit.paX$coefficients[1]
  CI[[2]] <- confint(fit.paX,1)
  SE[[2]] <- vcov(fit.paX)[1,1]
  
  if(length(O)==0){
    fit.O <- lm(Y~X-1)
  } else{
    Data.O<-data.frame[,union(x,O)]
    fit.O <- lm(Y~.-1,data=Data.O)
  }
  Beta[3] <- fit.O$coefficients[1]
  CI[[3]] <- confint(fit.O,1)
  SE[[3]] <- vcov(fit.O)[1,1]
  
  results <- list(x,y,paX,O,Beta,CI,SE)
  return(results)
}

# center with 'apply()'
center_apply <- function(x) {
  apply(x, 2, function(y) y - mean(y))
}

#sample function that works with single integer
resamp <- function(x,...){if(length(x)==1) x else sample(x,...)} 

#possible parent recovering function
possparents <- function(x, amat){
  posspaX <- which(amat[x,]==1)
  return(posspaX)
}

# parents recovering function
parents <- function(x, amat){
  paX <- setdiff(possparents(x,amat),which((amat[,x] + t(amat[x,]))==2))
  return(paX)
}

# descendants recovering function
de <- function(x, amat){
  deX <- setdiff(pcalg::possDe(amat,x,type="dag"),x)
  return(deX)
}

# O-set recovering function
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

