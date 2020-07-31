
#source("/u/henckell/Dokumente/Papers/research/new/DraftPapers/Leo/efficiency-adjustment/Rcode/simulation-paper.R")
setwd("~/Dokumente/Papers/Rcode/CodeEfficientVAS/Simulations")
source("~/Dokumente/Papers/Rcode/CodeEfficientVAS/MainFunctions")

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

# ptm <- proc.time()

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
  # save(file="/userdata/henckell/MSEquicktest.Rdata", Output)
  
  save(file=paste("/userdata/henckell/joint/resJointX",i,".Rdata",sep = ""),Output)
 #  save(file=paste("/Users/henckell/Dokumente/res",i,".Rdata",sep = ""),Output)
}

possible.x.size <- c(1) 

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
  # save(file="/userdata/henckell/MSEquicktest.Rdata", Output)
  
  save(file=paste("/userdata/henckell/joint/resSingleX",i,".Rdata",sep = ""),Output)
  #  save(file=paste("/Users/henckell/Dokumente/res",i,".Rdata",sep = ""),Output)
}





