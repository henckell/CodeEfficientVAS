library(readxl)
library(pcalg)
library(latex2exp)
library(ggplot2)
library(EnvStats)
library(gridExtra)
library(grid)

#### getting necessary functions

source("~/CodeEfficientVAS/RealDataAnalysis/helpFunctions.R")

#### preparing data
baseData1 <- read_xls("Data/1. cd3cd28.xls")
baseData2 <- read_xls("Data/2. cd3cd28icam2.xls") ## ignored by Mooji
baseData3 <- read_xls("Data/3. cd3cd28+aktinhib.xls")
baseData4 <- read_xls("Data/4. cd3cd28+g0076.xls")
baseData5 <- read_xls("Data/5. cd3cd28+psitect.xls") ## changes B according to Mooji
baseData6 <- read_xls("Data/6. cd3cd28+u0126.xls")
baseData7 <- read_xls("Data/7. cd3cd28+ly.xls") ## changes B according to Mooji
baseData8 <- read_xls("Data/8. pma.xls")
baseData9 <- read_xls("Data/9. b2camp.xls")

data.transform <- function(data){
  data <- as.matrix(data)
  data <- log(data)
  data <- center_apply(data)
  return(data)
}

baseData<- data.transform(baseData)

# splom(baseData)

# colnames(baseData) <- c()

#### adjacency matrices for our three graphs, as well as their tweaked versions
# Consensus graph
B <- t(matrix(c(0,0,0,0,0,0,0,1,1,0,0,
               1,0,0,0,0,0,0,1,1,0,0,
               0,0,0,0,1,0,0,0,0,0,0,
               0,0,1,0,1,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,
               0,1,0,0,0,0,0,1,0,0,0,
               0,0,0,0,1,0,0,1,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,
               0,0,1,1,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,1,1,0,0,
               0,0,0,0,0,0,0,1,1,0,0), ncol=11))

B.1 <- t(matrix(c(0,0,0,0,0,0,0,1,1,0,0,
               1,0,0,0,0,0,0,1,1,0,0,
               0,0,0,0,1,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,
               0,1,0,0,0,0,0,1,0,0,0,
               0,0,0,0,1,0,0,1,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,
               0,0,1,1,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,1,1,0,0,
               0,0,0,0,0,0,0,1,1,0,0), ncol=11))

B.2 <- t(matrix(c(0,0,0,0,0,0,0,1,1,0,0,
               1,0,0,0,0,0,0,1,1,0,0,
               0,0,0,0,1,0,0,0,0,0,0,
               0,0,1,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,
               0,1,0,0,0,0,0,1,0,0,0,
               0,0,0,0,1,0,0,1,0,0,0,
               0,0,0,0,0,0,0,0,0,0,0,
               0,0,1,1,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,1,1,0,0,
               0,0,0,0,0,0,0,1,1,0,0), ncol=11))

names<-c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")
 
rownames(B) <- names
colnames(B) <- names

rownames(B.1) <- names
colnames(B.1) <- names

rownames(B.2) <- names
colnames(B.2) <- names

#Sachs graph
C <- t(matrix(c(0,0,0,0,0,0,0,1,1,0,0,
                1,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,1,0,1,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,0,1,0,0,0,
                0,0,0,0,0,1,0,1,0,0,0,
                0,0,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0), ncol=11))
#names<-c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")

C.1 <- t(matrix(c(0,0,0,0,0,0,0,1,1,0,0,
                1,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,0,1,0,0,0,
                0,0,0,0,0,1,0,1,0,0,0,
                0,0,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0), ncol=11))

C.2 <- t(matrix(c(0,0,0,0,0,0,0,1,1,0,0,
                1,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,1,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,0,1,0,0,0,
                0,0,0,0,0,1,0,1,0,0,0,
                0,0,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0), ncol=11))

names<-c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")
rownames(C) <- names
colnames(C) <- names


#Mooiji acyclic graph
D <- t(matrix(c(0,1,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,1,0,0,0,0,1,0,0,
                0,0,0,0,1,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,1,0,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0), ncol=11))

D.1 <- t(matrix(c(0,1,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,1,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,1,0,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0), ncol=11))

D.2 <- t(matrix(c(0,1,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,1,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,1,0,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0), ncol=11))

names<-c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")
rownames(D) <- names
colnames(D) <- names

## main function to compute estimated variances for all appropriate pairs (X,Y)
MainFunction <- function(Data,B){
  baseData<- data.transform(Data)
  
  x<-1:11

  descendants <- lapply(x,de,B)
  length.de <- sapply(descendants,length)
  no.de <- which(length.de==0)
  res<-list()
  loi<-c()

  # removing all nodes that have no descendant from the consideration for X
  x.clean <- setdiff(x,no.de)
  
  # computing and saving all descendants for each node in the graph
  for(i in 1:length(x.clean)){
    j<-which(x==x.clean[i])
    soi <- descendants[[j]]
    loi[i]<-length.de[j]
   }

  pa <- X <- Y <- O <- list()
  betas <-matrix(numeric(sum(loi)*3),ncol=3)
  SE.e <- SE.pa <-  SE.O <- numeric(sum(loi))

  n1<-0
  n2<-0
  
  res <- list()
  
  # computing estimated variances for all apropriate pairs
  for(i in 1:length(x.clean)){
    j<-which(x==x.clean[i])
    soi <- descendants[[j]]
    loi[i]<-length.de[j]
  
    n2 <- n1 +n2
    n1<- 0
    for(k in 1:loi[i]){ 
      res[[k + n2]] <- CIeval(baseData,B,j,soi[k])
      tot <- CIeval(baseData,B,j,soi[k])
      betas[k + n2,] <- tot[[5]]
      # eCI[k + n2,] <- tot[[6]][[1]]
      X[[k + n2]] <- tot[[1]]
      Y[[k + n2]] <- tot[[2]]
      pa[[k + n2]] <- tot[[3]]
      O[[k + n2]] <- tot[[4]]
      SE.pa[k+n2] <- tot[[7]][[2]]
      SE.O[k+n2] <- tot[[7]][[3]]
      n1<-n1+1
    }
  }
  
  ratio <- SE.O/SE.pa
  
  return(list(var.pa=SE.pa,var.O=SE.O,ratio=ratio,X=X,Y=Y,pa=pa,O=O))
}


# computing all estimated variances for consensus graph in the 8 exp. conditions
B1 <- MainFunction(baseData1,B)
B3 <- MainFunction(baseData3,B)
B4 <- MainFunction(baseData4,B)
B5 <- MainFunction(baseData5,B.1)
B6 <- MainFunction(baseData6,B)
B7 <- MainFunction(baseData7,B.2)
B8 <- MainFunction(baseData8,B)
B9 <- MainFunction(baseData9,B)

# removing cases where O=P
b1 <- B1$ratio[-which(B1$ratio==1)]
b3 <- B3$ratio[-which(B3$ratio==1)]
b4 <- B4$ratio[-which(B4$ratio==1)]
b5 <- B5$ratio[-which(B5$ratio==1)]
b6 <- B6$ratio[-which(B6$ratio==1)]
b7 <- B7$ratio[-which(B7$ratio==1)]
b8 <- B8$ratio[-which(B8$ratio==1)]
b9 <- B9$ratio[-which(B9$ratio==1)]

# computing all estimated variances for Sachs et al. graph in the 8 exp. conditions
C1 <- MainFunction(baseData1,C)
C3 <- MainFunction(baseData3,C)
C4 <- MainFunction(baseData4,C)
C5 <- MainFunction(baseData5,C.1)
C6 <- MainFunction(baseData6,C)
C7 <- MainFunction(baseData7,C.2)
C8 <- MainFunction(baseData8,C)
C9 <- MainFunction(baseData9,C)

# removing cases where O=P
c1 <- C1$ratio[-which(C1$ratio==1)]
c3 <- C3$ratio[-which(C3$ratio==1)]
c4 <- C4$ratio[-which(C4$ratio==1)]
c5 <- C5$ratio[-which(C5$ratio==1)]
c6 <- C6$ratio[-which(C6$ratio==1)]
c7 <- C7$ratio[-which(C7$ratio==1)]
c8 <- C8$ratio[-which(C8$ratio==1)]
c9 <- C9$ratio[-which(C9$ratio==1)]


# computing all estimated variances for Mooij and Heskes graph in the 8 exp. conditions
D1 <- MainFunction(baseData1,D)
D3 <- MainFunction(baseData3,D)
D4 <- MainFunction(baseData4,D)
D5 <- MainFunction(baseData5,D.1)
D6 <- MainFunction(baseData6,D)
D7 <- MainFunction(baseData7,D.2)
D8 <- MainFunction(baseData8,D)
D9 <- MainFunction(baseData9,D)

# removing cases where O=P
d1 <- D1$ratio[-which(D1$ratio==1)]
d3 <- D3$ratio[-which(D3$ratio==1)]
d4 <- D4$ratio[-which(D4$ratio==1)]
d5 <- D5$ratio[-which(D5$ratio==1)]
d6 <- D6$ratio[-which(D6$ratio==1)]
d7 <- D7$ratio[-which(D7$ratio==1)]
d8 <- D8$ratio[-which(D8$ratio==1)]
d9 <- D9$ratio[-which(D9$ratio==1)]


## pooled violin plots with geomeans added as red dots and medians as black dot
y1 <- c(b1,b3,b4,b5,b6,b7,b8,b9)
y2 <- c(c1,c3,c4,c5,c6,c7,c8,c9)
y3 <- c(d1,d3,d4,d5,d6,d7,d8,d9)
Y <- c(y1,y2,y3)

dat <- data.frame(y=c(y1,y2,y3),x=c(rep(1,length(y1)),rep(2,length(y2)),rep(3,length(y3))))
df <- data.frame(y=c(geoMean(y1),geoMean(y2),geoMean(y3)),x=c(1,2,3))
df2  <- data.frame(y=c(median(y1),median(y2),median(y3)),x=c(1,2,3))

g <- ggplot(data=dat,aes(x=factor(x),y=y),fill=factor(x)) + geom_violin() + 
  ylim(0,2) + 
  geom_abline(slope=0,intercept=1,color = "blue") + 
  ggtitle("Ratio of estimated variances") +
  ylab(TeX("$\\widehat{var}(\\hat{\\beta}_{yx.\\mathbf{o}})/\\widehat{var}(\\hat{\\beta}_{yx.\\mathbf{p}})$")) + 
  # ylab("ytitle") + 
  # $\\widehat{var}(\\hat{\\beta}_{yx.\\mathbf{o}})/ \\widehat{var}(\\hat{\\beta}_{yx.\\mathbf{p}})$')) + 
  xlab("") + scale_color_brewer(palette="Dark2") +
  scale_x_discrete(labels=c(parse(text = TeX('$G_{C}$')),parse(text = TeX('$G_{S}$')),
                                   parse(text = TeX('$G_{M}$')))) +
                     # c("Consensus","Sachs et al.",TeX('$mode(L_{ij})$'))) +
  # "G0076","Psitectorigenin","U0126","LY204002","PMA","2CAMP")) + 
  theme(axis.text.x=element_text(size=12)) + 
  geom_point(shape=15,data=df,aes(x=x,y=y), colour = 'red',inherit.aes=FALSE) +
  geom_point(shape=15,data=df2,aes(x=x,y=y), colour = 'black',inherit.aes=FALSE) +
  theme(axis.title.x = element_text(size = 15, margin = margin(10,0,0,0))) + 
  theme(axis.title.y = element_text(size = 15, margin = margin(10,0,0,0)))

  pdf(file="JointViolin.pdf")
  g
  dev.off()
  