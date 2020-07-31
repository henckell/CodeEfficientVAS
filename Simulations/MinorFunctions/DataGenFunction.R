#### Data genaration functions: take an adjacency matrix and error parameters and return set of data 

## uniform errors
DataGenUnif<-function(sample.size,B,limit){
  p <- length(B[1,])
  X<-matrix(numeric(sample.size*p),ncol=p)
  X[,1] = runif(sample.size,-limit[1],limit[1])
  X[,2] = B[2,1] %*% X[,1] + runif(sample.size,-limit[2],limit[2])
  if(p>=3){
    for(i in 3:p){
      X[,i] <-X[,1:(i-1)]  %*% t(t(B[i,1:(i-1)])) + runif(sample.size,-limit[i],limit[i])
    }
  }
  return(X)
}

## student-t distributed errors
DataGenStu<-function(sample.size,B,df,v){
  p <- length(B[1,])
  X<-matrix(numeric(sample.size*p),ncol=p)
  X[,1] = rt(sample.size,df[1]) * sqrt(v[1] * (df[1]-2)/df[1])
  X[,2] = B[2,1] %*% X[,1] + rt(sample.size,df[2]) * sqrt(v[2] * (df[2]-2)/df[2])
  if(p>=3){
    for(i in 3:p){
      X[,i] <-X[,1:(i-1)]  %*% t(t(B[i,1:(i-1)])) + rt(sample.size,df[i]) * sqrt(v[i] * (df[i]-2)/df[i])
    }
  }
  return(X)
}

## logistically dsitributed errors
DataGenLog<-function(sample.size,B,scale){
  p <- length(B[1,])
  X<-matrix(numeric(sample.size*p),ncol=p)
  X[,1] = rlogis(sample.size,scale=scale[1])
  X[,2] = B[2,1] %*% X[,1] + rlogis(sample.size,scale=scale[2])
  if(p>=3){
    for(i in 3:p){
      X[,i] <-X[,1:(i-1)]  %*% t(t(B[i,1:(i-1)])) + rlogis(sample.size,scale=scale[i])
    }
  }
  return(X)
}
