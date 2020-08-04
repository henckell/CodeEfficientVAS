library(pcalg)
library(MASS)
library(gRbase)
library(ggm)

## this file contains 3 functions: 
##data.eval to compute the the total effect estimates for a given data set and graphs
##generateB: to help generating appropriate adjacency matrices
##MainFunction: to sample linear causal models, generate data and collect estimates to compute empirical MSEs

## function that given data and adjustment sets of interest computes the corresponding regression coefficients
data.eval <- function(x,y,data,cpdag.permitted,cpdag.desc,
                      paX.dag,adjust.dag,all.dag,all.nforb.dag,O.dag,
                      paX.cpdag,adjust.cpdag,all.cpdag,all.nforb.cpdag,O.cpdag,
                      allowed=0){
  
  k <- length(x)
  n<-6
  
  SingleBeta <- matrix(numeric(2*n*k),ncol=2*n)
  
  Data.Frame <- as.data.frame(data)
  
  X<-Data.Frame[,x]
  Y<-Data.Frame[,y]
  
  ## fit with empty conditioning set
  if(k==1){
    fit.Triv <- lm(Y~X-1)
  }  else{
    fit.Triv <- lm(Y~.-1,data=X)
  }
  SingleBeta[,1] <- fit.Triv$coefficients
  
  ## DAG
  SingleBeta[,2] <- lmLeo(Data.Frame,x,y,paX.dag,X,Y)
  
  SingleBeta[,3] <- lmLeo(Data.Frame,x,y,adjust.dag,X,Y)
  
  SingleBeta[,4] <- lmLeo(Data.Frame,x,y,all.dag,X,Y)
  
  SingleBeta[,5] <- lmLeo(Data.Frame,x,y,all.nforb.dag,X,Y)
  
  SingleBeta[,6] <- lmLeo(Data.Frame,x,y,O.dag,X,Y)
  
  
  ## CPDAG
  if(cpdag.permitted == 1 | allowed != 0){
    SingleBeta[,(n+1):(2*n)] <- NA
  } else{
    if(cpdag.desc==1){
      SingleBeta[,(n+1):(2*n)] <- 0
    } else {
      
      SingleBeta[,(n+1)] <- fit.Triv$coefficients[1:k]  
      SingleBeta[,(n+2)] <- lmLeo(Data.Frame,x,y,paX.cpdag,X,Y)
      SingleBeta[,(n+3)] <- lmLeo(Data.Frame,x,y,adjust.cpdag,X,Y)
      SingleBeta[,(n+4)] <- lmLeo(Data.Frame,x,y,all.cpdag,X,Y)
      SingleBeta[,(n+5)] <- lmLeo(Data.Frame,x,y,all.nforb.cpdag,X,Y)
      SingleBeta[,(n+6)] <- lmLeo(Data.Frame,x,y,O.cpdag,X,Y)
    }
  }
  return(SingleBeta)
}


# function to generate appropriate adjacency matrices
generateB <- function(size, neigh, graph.type){
  G <- randDAG(size, neigh, graph.type)
  
  ## compute edge-weight matrix, which is topologically ordered 
  B.dag <- t(as(G,"matrix"))
  
  subset<-which(B.dag!=0)
  B.dag[subset]=runif(length(subset),0.1,2)
  
  randomMatrix<-matrix(rbinom(size^2,1,0.5),ncol=size)
  subset<-which(randomMatrix==0)
  randomMatrix[subset]<--1
  B.dag=B.dag*randomMatrix
  
  ## topologically order B
  topo.order <- topoSort(G)
  topo.order <- unlist(topo.order)
  
  B.dag <- B.dag[topo.order,topo.order]
  colnames(B.dag) <- 1:size
  rownames(B.dag) <- 1:size
  
  ## topologically order G :)  
  #G<-graph.adjacency(t(B.dag),weighted=TRUE)
  #G <- igraph:::as_graphnel(G)
  
  ## compute adjacency matrix
  return(B.dag)
}

############ 

MainFunction <- function(size, neigh, sample.size, reps, graph.type="er",
                         error.type="normal",x.size=1){
  
  ## setting number of sets because less work this way in the future to change
  n <- 6
  
  #### generate DAG and compute matters of interest from it 
  set<- 1:size
  DescendantsSize <- 0
  allowed <- 1
  tcpdag.permitted <- 0
  
  B.dag <- generateB(size, neigh, graph.type)

  amat.dag <- B.dag
  amat.dag[B.dag != 0] <- 1
  
  ## loop that draws x till de(x) \neq \emptyset
  while(DescendantsSize==0 | allowed!=0 | tcpdag.permitted == 0){
    #if(length(set)==0) break
    
    if(length(set)<=x.size-1){
      set<-1:size
      DescendantsSize <- 0
      allowed <- 1
      
      B.dag <- generateB(size, neigh, graph.type)
      
      amat.dag <- B.dag
      amat.dag[B.dag != 0] <- 1
    }
  
    ##ema: changed to sample x
    x <- resamp(set,x.size)

    ##get desc of x
    deX.dag <- c()  
    deX.dag <- lapply(x,2,FUN=pcalg::possDe,m=B.dag,type="dag")
    deX.dag <- unique(Reduce(intersect, deX.dag))
    
    
    deX.dag <- setdiff(deX.dag,x)
    DescendantsSize <- length(deX.dag)
    tcpdag.permitted <- 0
    
    if(DescendantsSize==0){set = setdiff(set,x)}
    else{
      ## drawing a y that is a descendant of x
      y <- resamp(deX.dag,1)
      
      true.cpdag <- dag2essgraph(B.dag)
      amat.truecpdag <- t(as(true.cpdag,"matrix"))
      amat.truecpdag[amat.truecpdag != 0] <- 1
      
      tcpdag.permitted <- pcalg:::isAmenable(amat.truecpdag,x,y,type="cpdag") 
      
      Forbidden <- pcalg:::forbiddenNodes(amat.dag,x,y)
      allowed <- length(intersect(Forbidden,c(x)))
      }
    
    if(tcpdag.permitted == 0 | allowed >= 1){
      DescendantsSize <- 0
      set <- setdiff(set, x) 
    }
    

  }
  
  Forbidden <- pcalg:::forbiddenNodes(amat.dag,x,y)
  
  ## computing sets 4 of interest for DAG
  paX.dag <- lapply(x,parents,amat=amat.dag)
  paX.dag <- unique(Reduce(union, paX.dag))
  paX.dag <- setdiff(paX.dag,x)
  
  #Forbidden <- pcalg:::forbiddenNodes(amat.dag,x,y)

  adjust.dag <- lapply(union(x,y),pcalg:::possAn,m=amat.dag,type="dag",possible=FALSE)
  adjust.dag <- unique(Reduce(union, adjust.dag))
  adjust.dag <- setdiff(adjust.dag,Forbidden)
  adjust.dag <- setdiff(adjust.dag,union(x,y))
  
  
  all.dag <- setdiff(1:size,union(x,y))

  all.nforb.dag <- setdiff(all.dag, Forbidden) 
  
  O.dag <- ComputeO(x,y,amat.dag)


  ## computing causal effect
  CE <- causalEffectLeo(B.dag,y,x)
  ####
  
  #### true CPDAG

  
  ## store sizes of sets of interest
  sizes <- numeric(n-1)
  sizes[1] <- length(paX.dag)
  sizes[2] <- length(adjust.dag)
  sizes[3] <- length(all.dag)
  sizes[4] <- length(all.nforb.dag)
  sizes[5] <- length(O.dag)
  
  ## check whether sets are VAS

  
  #print(adj.dag)
  
  names(sizes) <- c("paX DAG","adjust DAG","all DAG","all n-forb DAG","O DAG")
  
  ### checking whether the sets of interest are actually valid
  #valid.truecpdag <- numeric(3)
  #valid.truecpdag[1] <- gac(amat.truecpdag,x,y,posspaX.truecpdag, type= "cpdag")$gac
  #valid.truecpdag[2] <- gac(amat.truecpdag,x,y,paX.truecpdag, type= "cpdag")$gac
  #valid.truecpdag[3] <- gac(amat.truecpdag,x,y,O.truecpdag, type= "cpdag")$gac

  
  #### generate data and estimate parameters of interest
    
  descend.issues <-numeric(reps) 
  amen.issues <- numeric(reps)
  
  #Beta <- matrix(numeric(7*reps),ncol=7)

  Beta <- lapply(1:x.size,function(x) matrix(numeric(2*n*reps),ncol=2*n))
  size.cpdag <- matrix(numeric((n-1)*reps),ncol=(n-1))
  allowed.cpdag <- numeric(reps)
   
  #### preparing data simulation ####
  
  ## normal
  if(error.type == "normal"){
    ## randomly drawing error variances  
    epsilon <- runif(size,0.5,1.5)
    
    ## compute covariance matrix for G given epsilon
    CovG <- Covariance(B.dag,epsilon)
  }
  
  ## uniform
  if(error.type == "unif"){
    ## randomly drawing error variances  
    limit <- runif(size,1.2,2.1)
  }
  
  ## student
  if(error.type == "stu"){
    ## randomly drawing error variances  
    df <- rep(5,size)
    v <- runif(size,0.5,1.5)
  }
  
  ## logistic 
  if(error.type == "log"){
    ## randomly drawing error variances  
    scale <- runif(size,0.4,0.7)
  }
  
  #### simulating and evaluating data ####
  
  for(i in 1:reps){
    
    if(error.type == "unif"){
      Data <- DataGenUnif(sample.size,B.dag,limit)
    }
    
    if(error.type == "stu"){
      Data <- DataGenStu(sample.size,B.dag,df,v)
    }
    
    if(error.type == "log"){
      Data <- DataGenLog(sample.size,B.dag,scale)
    }

    if(error.type=="normal"){
    ## generated data given covariance matrix of the DAG
    Data <- mvrnorm(sample.size,rep(0,size),CovG)
    
    ## tranform data to format required by ges
    score <- new("GaussL0penObsScore", Data)
    
    ## estimating the CPDAG C
    summary <- ges(score)
    cpdag <- summary$essgraph
    
    ## compute the amat of C
    amat.cpdag <- t(as(cpdag,"matrix"))
    amat.cpdag[amat.cpdag != 0] <- 1
    } else{
      amat.cpdag <- lingam(Data)$Bpruned
       amat.cpdag[amat.cpdag != 0] <- 1
    }
      
    
    ## store whether C is amenable
    cpdag.permitted <- 1 - pcalg:::isAmenable(amat.cpdag,x,y,type="cpdag")
    amen.issues[i] <- cpdag.permitted


    ## store whether y is an ancestor of x in the estimated CPDAG
    AnY <- pcalg:::possibleAnProper(amat.cpdag,y)
    y.desc.x <- length(intersect(x,AnY))
    cpdag.desc <- c(0)
    if(y.desc.x==0) cpdag.desc <- 1
    descend.issues[i] <- cpdag.desc
    
    ## compute sets of interest for estimated CPDAG
    O.cpdag <- ComputeO(x,y,amat.cpdag)
    
    paX.cpdag <- lapply(x,parents,amat=amat.cpdag)
    paX.cpdag <- unique(Reduce(union, paX.cpdag))
    paX.cpdag <- setdiff(paX.cpdag,union(x,y))

    # store whether a valid adjustment set exists 
    Forbidden.cpdag <- pcalg:::forbiddenNodes(amat.cpdag,x,y)
    allowed.cpdag[i] <- length(intersect(Forbidden.cpdag,c(x)))
    
    adjust.cpdag <- lapply(union(x,y),pcalg:::possAn,m=amat.cpdag,type="cpdag",
                           possible=TRUE)
    adjust.cpdag <- unique(Reduce(union, adjust.cpdag))
    adjust.cpdag <- setdiff(adjust.cpdag,Forbidden.cpdag)
    adjust.cpdag <- setdiff(adjust.cpdag,union(x,y))
    
    all.cpdag <- setdiff(1:size,union(x,y))
    
    all.nforb.cpdag <- setdiff(all.cpdag,Forbidden.cpdag)
    
    ## stores sets of interest sizes
    if(cpdag.permitted==1 | allowed.cpdag[i] != 0){
      size.cpdag[i,1] <- NA
      size.cpdag[i,2] <- NA
      size.cpdag[i,3] <- NA
      size.cpdag[i,4] <- NA
      size.cpdag[i,5] <- NA
    }else{
      size.cpdag[i,1] <- length(paX.cpdag)
      size.cpdag[i,2] <- length(adjust.cpdag)
      size.cpdag[i,3] <- length(all.cpdag)
      size.cpdag[i,4] <- length(all.nforb.cpdag)
      size.cpdag[i,5] <- length(O.cpdag)
      }

    
    ## estimating parameters of interest
    for(k in 1:x.size){
    Beta[[k]][i,]<- data.eval(x,y,Data,amen.issues[i],descend.issues[i],  
                         paX.dag,
                         adjust.dag,
                         all.dag,
                         all.nforb.dag,
                         O.dag,
                         paX.cpdag,
                         adjust.cpdag,
                         all.cpdag,
                         all.nforb.cpdag,
                         O.cpdag,allowed.cpdag[i])[k,]
    }
  }
  
  colnames(size.cpdag) <- c("paX","adjust","all","all n-forb","O")
  
  ####  evaluating the computed regression coefficients 
  for(k in 1:x.size){
  colnames(Beta[[k]]) <- c("empty","paX DAG","adjust DAG","all DAG","all n-forb DAG","O DAG",
                           "paX est. CPDAG","adjust est. CPDAG","all est. CPDAG","all n-forb est. CPDAG",
                      "O est. CPDAG","empty est. CPDAG")
  }
  
  Bias.l <- lapply(1:x.size, function(k) (apply(Beta[[k]],2,mean,na.rm = TRUE) - CE[k]))
  Var.l <- lapply(1:x.size, function(k) (apply(Beta[[k]],2,stats::var,na.rm = TRUE)))
  MSE.l <- lapply(1:x.size, function(k) (Bias.l[[k]]^2 + Var.l[[k]]))
  
  Bias <- matrix(numeric(x.size*2*n),nrow=x.size)
  Var <- matrix(numeric(x.size*2*n),nrow=x.size)
  MSE <- matrix(numeric(x.size*2*n),nrow=x.size)
  
  for(i in 1:x.size){
    Bias[i,] <- Bias.l[[i]]
    Var[i,] <- Var.l[[i]]
    MSE[i,] <- MSE.l[[i]]
  }
  
  colnames(Bias) <- c("empty","paX DAG","adjust DAG","all DAG","all n-forb DAG","O DAG",
                      "empty est. CPDAG",
                      "paX est. CPDAG","adjust est. CPDAG","all est. CPDAG","all n-forb est. CPDAG",
                      "O est. CPDAG")
  colnames(Var) <- c("empty","paX DAG","adjust DAG","all DAG","all n-forb DAG","O DAG",
                     "empty est. CPDAG",
                     "paX est. CPDAG","adjust est. CPDAG","all est. CPDAG","all n-forb est. CPDAG",
                     "O est. CPDAG")
  colnames(MSE) <- c("empty","paX DAG","adjust DAG","all DAG","all n-forb DAG","O DAG",
                     "empty est. CPDAG",
                     "paX est. CPDAG","adjust est. CPDAG","all est. CPDAG","all n-forb est. CPDAG",
                     "O est. CPDAG")
  
  
  #Results <- cbind(Bias,Var,MSE)
  #colnames(Results) <- c("Bias","Var","MSE")
  
  Issues <- cbind(amen.issues,descend.issues)
  colnames(Issues) <- c("not amen","not desc")
  
  return(list(graph.type=graph.type,error.type=error.type,x.size=x.size,x=x,y=y,
              Bias=Bias, Var=Var, MSE=MSE,
              size=size,neigh=neigh, sample.size=sample.size,
              CE=CE,Beta=Beta,
              Sizes=sizes,
              size.cpdag = size.cpdag,Issues=Issues,
              allowed.cpdag = allowed.cpdag
              ))
}

## Example run
# MainFunction(10,2,1000,100)
