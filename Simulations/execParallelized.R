setwd(".../CodeEfficientVAS")
source("Simulations/MainFunctions.R")
source("Simulations/MinorFunctions/helpFunctions.R")
source("Simulations/MinorFunctions/DataGenFunction.R")

library("pcalg")
library("MASS")
library("foreach")
library("doParallel")


#setting up the parameters
nc <- 30
registerDoParallel(nc) ## Request nc processors in parallel
possible_size<-c(10,20,50,100)
possible_neigh<-c(2,3,4,5)
possible_sample.size<-c(125,500,2000,8000)
reps<-100
possible_graph.type<-c("er","power")
possible.x.size <- c(2,3) # <- c(1) ## change for singleton X simulation

num_setings  <-  5000
seed1<-1

# |X| either two or three simulation

finalRes <- foreach(i=seed1: num_setings,.packages = 'pcalg') %dopar%{
  set.seed(i)
  
  reps<-100
  
  size <- sample(possible_size,1)
  neigh <- sample(possible_neigh,1)
  graph.type <- sample(possible_graph.type,1)
  sample.size <- sample(possible_sample.size,1) 
  x.size <- sample(possible.x.size,1)
  
  error.type <- sample(c("normal","not-normal"),1)
  
  if(error.type=="not-normal") {
    error.type <- sample(c("unif","stu","log"),1)
  }
  
  Res<-MainFunction(size,neigh,sample.size,reps,graph.type,error.type,x.size=x.size)
  

  
  Output <- list(seed=i,size=size,neigh=neigh,sample.size=sample.size,x.size=x.size,
                 graph.type=Res$graph.type,error.type=Res$error.type,
                 x=Res$x,y=Res$y,CE=Res$CE,
                 Beta=Res$Beta,Bias=Res$Bias, Var=Res$Var, MSE=Res$MSE,
                 Issues=Res$Issues,adj.dag = Res$adj.dag, adj.cpdag = Res$adj.cpdag, allowed.cpdag = Res$allowed.cpdag,
                 Sizes=Res$Sizes,Sizes.cpdag=Res$size.cpdag)
  
  save(file=paste("data/resJointX",i,".Rdata",sep = ""),Output)
}

# |X|=1 simulation

seed1<-1 # reset seed

finalRes <- foreach(i=seed1: num_setings,.packages = 'pcalg') %dopar%{
  set.seed(i)
  
  reps<-100
  
  size <- sample(possible_size,1)
  neigh <- sample(possible_neigh,1)
  graph.type <- sample(possible_graph.type,1)
  sample.size <- sample(possible_sample.size,1) 
  x.size <- 1
  
  error.type <- sample(c("normal","not-normal"),1)
  
  if(error.type=="not-normal") {
    error.type <- sample(c("unif","stu","log"),1)
  }
  
  Res<-MainFunction(size,neigh,sample.size,reps,graph.type,error.type,x.size=x.size)
  
  
  
  Output <- list(seed=i,size=size,neigh=neigh,sample.size=sample.size,x.size=x.size,
                 graph.type=Res$graph.type,error.type=Res$error.type,
                 x=Res$x,y=Res$y,CE=Res$CE,
                 Beta=Res$Beta,Bias=Res$Bias, Var=Res$Var, MSE=Res$MSE,
                 Issues=Res$Issues,adj.dag = Res$adj.dag, adj.cpdag = Res$adj.cpdag, allowed.cpdag = Res$allowed.cpdag,
                 Sizes=Res$Sizes,Sizes.cpdag=Res$size.cpdag)
  
  save(file=paste("data/resSingleX",i,".Rdata",sep = ""),Output)
}





