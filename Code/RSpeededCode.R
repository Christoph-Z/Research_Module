library(profvis)


##--------------------------------------------------------------------
##Generating data

##Sample size
n <- 50000
##Number of regressors
p <- 5
##Number of regressors unequal to zero
p.True <-2
##True Beta Vector
beta <- c(1.5,3,0,0,0)

##Generate data for X
x_1 <- rep(1,n)         
x_2 <- rnorm(n,2,1)
x_3 <- rnorm(n,0,1)
x_4 <- rnorm(n,5,2)
x_5 <- rnorm(n,4,1)

X.all <- cbind(x_1,x_2,x_3,x_4,x_5)
y.mean <- X.all%*%beta
y <- X.all%*%beta+rnorm(n,0,1)
X <- X.all

##--------------------------------------------------------------------
##Basic speed comparison

##Calculate Pmatrix 
Rprof("Pmatrix default")
X%*%solve(t(X)%*%X)%*%t(X)%*%y
Rprof(NULL)

Rprof("Pmatrix")
X%*%(solve(t(X)%*%X)%*%(t(X)%*%y))
Rprof(NULL)

##Print Summary
summaryRprof("Pmatrix default")
summaryRprof("Pmatrix")

##Sample Generating for Bootstrap
Rprof("Sample")
sample(seq(1,n,1),n-10)
Rprof(NULL)

Rprof("Sample.int")
seq(1,n,1)[sample.int(n,n-10)]
Rprof(NULL)

##Print Summary
summaryRprof("Sample")
summaryRprof("Sample.int")

##--------------------------------------------------------------------
##Set of all possible models

Models <- function(p.True,p){
  ##For the CV we need to compute the set A of all possible models 
  #Creats the set {1,...,p} from which we want to generate the Powerset
  Index <- seq(1,p,1)                        
  
  #Denotes the number of Possible Models out of {1,...,p} Regressors 
  col.A <- 2^p-1                                           
  
  #Denote A as Powerset of {1,...,p}
  A <- matrix(0L,nrow = p, ncol = col.A)      
  k <- choose(p,1)
  l <- 1
  for (i in 1:p) {
    #combn spits out all combinations of i elements in Index
    A[1:i,l:k] <- combn(Index,i)             
    k <-k + choose(p,i+1)
    l <- l + choose(p,i)
  }
  return(A)
}

A <- Models(p.True,5)

##Set of Cat II Models
Partion <- function(p.True,p){
  ##For the CV we need to compute the set A of all possible models 
  #Creats the set {1,...,p} from which we want to generate the Powerset
  Index <- seq(1,p,1)                        
  
  #Denotes the number of Possible Models out of {1,...,p} Regressors 
  col.A <- 2^p-1                                           
  
  #Denote A as Powerset of {1,...,p}
  A <- matrix(0L,nrow = p, ncol = col.A)      
  k <- choose(p,1)
  l <- 1
  for (i in 1:p) {
    #combn spits out all combinations of i elements in Index
    A[1:i,l:k] <- combn(Index,i)             
    k <-k + choose(p,i+1)
    l <- l + choose(p,i)
  }
  
  
  ##Split A into two dijoint subsets, ie, the set of Category I Models and Category II Models 
  ##
  ##To Compute the Set of CI and CII we assume a certain strurcture on the Data. We need X to be s.t all
  ##"TRUE" regrossers are in the first column of the X Matrix, i.e, if we have P "true" regressors, they
  ##are given throght X[1],...,X[P]. And evrey X[P+] is trash...
  ##
  ##We may need this sets later n for Simulation study's
  coln <- c()
  for (i in p.True:col.A) {
    if(all(seq(1,p.True,1) %in% A[,i])){    
      coln <- c(coln,i)
    }
  }
  C2 <- A[,coln]           #Set of all Category II Models
  C1 <-A[,-coln]           #Set of all Category I Models
  return(C2)
}


##--------------------------------------------------------------------
##AIC and BIC

Pmatrix <- list()
for (i in 1:length(A[1,])) {
  X.train <- X[,A[,i]]
  X.train.transpos <- t(X.train)
  Pmatrix[[paste0("element",i)]] <- X.train%*%(solve(X.train.transpos%*%X.train)%*%X.train.transpos)
}



##AIC/BIC
InfoCrit <- function(y,X,I,A,Criterion = "AIC"){
  ##Calculate Information Criterion
  
  ##For AIC
  if(Criterion == "AIC"){
    
    j <- 0
    InfoCriterion <- c()
    for (i in I) {
      j <- j+1
      
      ##Number of regressors
      k <- sum(A[,j] != 0)
      
      ##Expectation observation j
      mu <- i%*%y

      ##Varianz
      var <- 1/n*sum((y-mu)^2)
      
      ##Calculate Log Likelihood value
      LogL <- -n/2*(log(2*pi)+log(var)+1)
      
      InfoCriterion[j] <- 2*k-2*LogL
    }
  }
  ##For BIC
  if(Criterion == "BIC"){
    j <- 0
    InfoCriterion <- c()
    for (i in I) {
      j <- j+1
      
      ##Number of regressors
      k <- sum(A[,j] != 0)
      
      ##Expectation observation j
      mu <- i%*%y
      
      ##Varianz
      var <- 1/n*sum((y-mu)^2)
      
      ##Calculate Log Likelihood value
      LogL <- -n/2*(log(2*pi)+log(var)+1)
      
      InfoCriterion[j] <- log(n)*k - 2*LogL
    }
  }
  
  #Choose the Model which minimze the Information Criterion
  TheChosenOne <- which.min(InfoCriterion)
  return(A[,TheChosenOne])
}

Rprof("AIC")
InfoCrit(y,X,Pmatrix,A)
Rprof(NULL)
summaryRprof("AIC")


##--------------------------------------------------------------------
##CV1

#Combination of Sample Splits
train <- combn(seq(1,n,1),n-1)

Pmatrix <- list()
k <- 0
for (i in 1:length(A[1,])) {
  A.i <- A[,i]
  for (j in 1:length(train[1,])) {
    k <- k+1
    train.j <- train[,j]
    X.train <- X[train.j,A.i] 
    Pmatrix[[paste0("element",k)]] <- X[-train.j,A.i]%*%solve(t(X.train)%*%X.train)%*%t(X.train)
  }
}

CV1 <- function(n,y,X,I,A,train){
  #Compute Prediction Errors for all sets in A
  MeanPred.Error <- c()
  k <- 0
  for (i in 1:length(A[1,])) {
    Pred.Error <- c()
    for (j in 1:n) {
      k <- k+1
      #To make the Code a bit shorter
      train.j <- train[,j]                 
      #Prediction Error for a given alpha and a given subset
      for(v in I[k]){
        Error.temp <- y[-train.j]-v%*%y[train.j]
        Pred.Error[j] <- t(Error.temp)%*%Error.temp
      }
    }
    MeanPred.Error[i] <- 1/length(n)*sum(Pred.Error)
  }
  TheChosenOne <- which.min(MeanPred.Error)
  return(A[,TheChosenOne])
}

Rprof("CV1")
CV1(n,y,X,Pmatrix,A,train)
Rprof(NULL)
summaryRprof("CV1")

##--------------------------------------------------------------------
##CV1 Shao
col.A <- 2^p-1
gram <- t(X)%*%X

proj_diag <- matrix(0L,nrow = n, ncol = col.A)
for(j in 1:col.A){
  for(i in 1:n){
    alpha <- A[A[,j]!=0,j]
    proj_diag[i,j] <- t(X[i,alpha])%*%solve(gram[alpha,alpha])%*%X[i,alpha]
  }
}

LOOCV <- function(y,B,C,D, Alpha = A){
  #this functions computes the regressors chosen by LOOCV
  #
  #y            Data for dependent variable
  #B            Data for Regressors
  #C            gram matrix of regressors
  #D            diagonal elements of the projection matrix
  #
  #Alpha        Set of possible Modelvaritaions for which the CV should be calculated 
  #             (By Default use all possible Models)  
  
  p <- length(X[1,])#Number of Possible Regressors
  n <- length(y)#Number of Observations
  col.A <- length(A[1,])
  
  #compute the OLS estimates for model alpha
  b <- matrix(0L, nrow = p, ncol = col.A)
  for(j in 1:col.A){
    b[A[A[,j]!=0,j],j]<- solve(C[A[A[,j]!=0,j],A[A[,j]!=0,j]])%*%t(B[,A[A[,j]!=0,j]])%*%y
  }
  
  #Compute Prediction Errors for all sets in A
  MeanPred.Error <- numeric(col.A)
  for (i in 1:col.A) {
    alpha <- A[A[,i]!=0,i]
    if(length(alpha)==1){
      MeanPred.Error[i] <- t((rep(1,n)-D[,i])^(-2))%*%((y - B[,alpha]*b[alpha,i])^2)/n
    }else{
      MeanPred.Error[i] <- t((rep(1,n)-D[,i])^(-2))%*%((y - B[,alpha]%*%b[alpha,i])^2)/n 
    }
  }
  TheChosenOne <- which.min(MeanPred.Error)
  return(A[,TheChosenOne])
}

Rprof("LOOCV")
y <- X%*%beta + rnorm(n,0,1)
LOOCV(y,X,gram,proj_diag,A)
Rprof(NULL)
summaryRprof("LOOCV")

##--------------------------------------------------------------------
##MCCV

MCCV <- function(n,n_v, y, X, A, b, Replacement = FALSE){  #n_v = #leaved out data points
  ##Combinations of Sample partions for fitting the model
  
  ##For the Mone Carlo CV with b subsets of {1,...,n}

  train <- matrix(ncol = b, nrow = n-n_v)
  for (i in 1:b) {
    train[,i] <- seq(1,n,1)[sample.int(n,n-n_v,replace = Replacement)] 
  }
  
  #Compute Prediction Errors for all sets in A
  MeanPred.Error <- c()
  for (i in 1:length(A[1,])) {
    Pred.Error <- c()
    for (j in 1:length(train[1,])) {
      #To make the Code a bit shorter
      train.j <- train[,j]                 
      X.train <- X[train.j,A[,i]]    
      #Prediction Error for a given alpha and a given subset
      Error.temp <- y[-train.j]-X[-train.j,A[,i]]%*%(solve(t(X.train)%*%X.train)%*%(t(X.train)%*%y[train.j]))
      Pred.Error[j] <- t(Error.temp)%*%Error.temp
    }
    MeanPred.Error[i] <- 1/length(train[1,])*sum(Pred.Error)
  }
  TheChosenOne <- which.min(MeanPred.Error)
  return(A[,TheChosenOne])
}


##Combination of Samplesplits
train.MCCV <- matrix(ncol = b, nrow = n-n_v)
for (i in 1:b) {
  train[,i] <- sample(seq(1,n,1),n-n_v,replace = FALSE)
}

##Pmatrix
Pmatrix_MCCV <- list()
k <- 0
for (i in 1:length(A[1,])) {
  A.i <- A[,i]
  for (j in 1:length(train.MCCV[1,])) {
    k <- k+1
    train.j <- train.MCCV[,j]
    X.train <- X[train.j,A.i] 
    Pmatrix_MCCV[[paste0("element",k)]] <- X[-train.j,A.i]%*%(solve(t(X.train)%*%X.train)%*%t(X.train))
  }
}

MCCV.v2 <- function(n,n_v, y, X,I, A, b,train, Replacement = FALSE){  #n_v = #leaved out data points
  #Compute Prediction Errors for all sets in A
  MeanPred.Error <- c()
  k <- 0
  for (i in 1:length(A[1,])) {
    Pred.Error <- c()
    for (j in 1:length(train[1,])) {
      k <- k+1
      #To make the Code a bit shorter
      train.j <- train[,j]                 
      #Prediction Error for a given alpha and a given subset
      for(v in I[k]){
        Error.temp <- y[-train.j]-v%*%y[train.j]
        Pred.Error[j] <- t(Error.temp)%*%Error.temp
      }
    }
    MeanPred.Error[i] <- 1/length(n_v)*sum(Pred.Error)
  }
  TheChosenOne <- which.min(MeanPred.Error)
  return(A[,TheChosenOne])
}


##--------------------------------------------------------------------
##Simulation 1

##Vector of diffrent Sample sizes
N <- seq(15,300,2)

##Number of Repetitions of the Experiment
m <- 2000

##Set of posssible models
A <- Models(2,5)

##Vector of probabilties of frequencys of choosing a Cat II Model for each samplezize by m repetations 
BIC <- rep(0,length(N))
AIC <- rep(0,length(N))
CV <- rep(0,length(N))
MCV <- rep(0,length(N))

for (q in 1:length(N)) {
  n <- N[q]
  n_v <- floor(n-n^(3/4))
  
  ##Number of Samplepartions
  b <- 5*n
  X <- X.all[1:n,]
  
  ##Pmatrix AIC/BIC
  Pmatrix_INF <- list()
  for (i in 1:length(A[1,])) {
    X.train <- X[,A[,i]]
    X.train.transpos <- t(X.train)
    Pmatrix_INF[[paste0("element",i)]] <- X.train%*%(solve(X.train.transpos%*%X.train)%*%X.train.transpos)
  }
  
  ##Pmatrix CV1
  train <- combn(seq(1,n,1),n-1)
  
  Pmatrix_CV <- list()
  k <- 0
  for (i in 1:length(A[1,])) {
    for (j in 1:length(train[1,])) {
      k <- k+1
      train.j <- train[,j]
      X.train <- X[train.j,A[,i]] 
      Pmatrix_CV[[paste0("element",k)]] <- X[-train.j,A[,i]]%*%solve(t(X.train)%*%X.train)%*%t(X.train)
    }
  }
  
  ##Combination of Samplesplits
  train.MCCV <- matrix(ncol = b, nrow = n-n_v)
  for (i in 1:b) {
    train.MCCV[,i] <- sample(seq(1,n,1),n-n_v,replace = FALSE)
  }
  
  ##Pmatrix
  Pmatrix_MCCV <- list()
  k <- 0
  for (i in 1:length(A[1,])) {
    A.i <- A[,i]
    for (j in 1:length(train.MCCV[1,])) {
      k <- k+1
      train.j <- train.MCCV[,j]
      X.train <- X[train.j,A.i] 
      Pmatrix_MCCV[[paste0("element",k)]] <- X[-train.j,A.i]%*%(solve(t(X.train)%*%X.train)%*%t(X.train))
    }
  }
  
  for (j in 1:m) {
    ##Generating Data
    y <- y.mean[1:n]+rnorm(n,0,1)
    
    #Define some temporary Variables, i.e, define for a given Method the in iteration i choosed Model
    temp.BIC <- InfoCrit(y,X,Pmatrix_INF,A,Criterion = "BIC")
    temp.AIC <- InfoCrit(y,X,Pmatrix_INF,A) 
    temp.CV1 <- CV1(n,y,X,Pmatrix_CV,A,train)
    temp.MCV <- MCCV.v2(n,n_v,y,X,Pmatrix_MCCV,A,b,train.MCCV)
    
    ##Counting how many times a Cat II Model is picked for n= N[i] in m interations
    BIC[q] <- BIC[q] + (sum( temp.BIC < (p.True+1) & temp.BIC >0) == p.True)
    AIC[q] <- AIC[q] + (sum( temp.AIC < (p.True+1) & temp.AIC > 0) == p.True)
    CV[q] <- CV[q] + (sum( temp.CV1 < (p.True+1) & temp.CV1 > 0) == p.True)
    MCV[q] <- MCV[q] + (sum( temp.MCV < (p.True+1) & temp.MCV > 0) == p.True)
  }
  
}


##--------------------------------------------------------------------
##Simulation 2


##Samplesize
n <- 500

##Number of Repetitions of the Experiment
m <- 2000

##Number of True regressors 
p.True<- 2

##Number of Regressors
p <- 5

##Calculate the set of Category II Models
C2 <- Partion(p.True,p)

##Gives the Dimenson of the choosen Model for each iteration
BIC <- c()
AIC <- c()
CV <- c()
MCV <- c()
BICV <- c()

##Number of leaved out data and number of samplepartions
n_v <- 394
b <- 5*n

##Calculate BICD for BICV with ibd package
BIBD <- ibd(n,b,n-n_v,pbar = TRUE)$N

##Calculationg Inverses
X <- X.all[1:n,]

##Pmatrix AIC/BIC
Pmatrix_INF <- list()
for (i in 1:length(C2[1,])) {
  X.train <- X[,C2[,i]]
  X.train.transpos <- t(X.train)
  Pmatrix_INF[[paste0("element",i)]] <- X.train%*%(solve(X.train.transpos%*%X.train)%*%X.train.transpos)
}

##Pmatrix CV1
train <- combn(seq(1,n,1),n-1)

Pmatrix_CV <- list()
k <- 0
for (i in 1:length(C2[1,])) {
  for (j in 1:length(train[1,])) {
    k <- k+1
    train.j <- train[,j]
    X.train <- X[train.j,C2[,i]] 
    Pmatrix_CV[[paste0("element",k)]] <- X[-train.j,C2[,i]]%*%solve(t(X.train)%*%X.train)%*%t(X.train)
  }
}

##Combination of Samplesplits
train.MCCV <- matrix(ncol = b, nrow = n-n_v)
for (i in 1:b) {
  train.MCCV[,i] <- sample(seq(1,n,1),n-n_v,replace = FALSE)
}

##Pmatrix
Pmatrix_MCCV <- list()
k <- 0
for (i in 1:length(C2[1,])) {
  A.i <- C2[,i]
  for (j in 1:length(train[1,])) {
    k <- k+1
    train.j <- train.MCCV[,j]
    X.train <- X[train.j,A.i] 
    Pmatrix_MCCV[[paste0("element",k)]] <- X[-train,A.i]%*%(solve(t(X.train)%*%X.train)%*%t(X.train))
  }
}

for (i in 1:m) {
  ##Generating Data
  y <- y.mean + rnorm(n,0,1)
  
  BIC[i] <- sum( InfoCrit(y,X,Pmatrix_INF,C2,Criterion = "BIC") != 0)
  AIC[i] <- sum( InfoCrit(y,X,Pmatrix_INF,C2,Criterion = "AIC") != 0)
  CV[i] <- sum( CV1(n,y,X,Pmatrix_CV,C2,train) != 0)
  MCV[i] <- sum( MCCV.v2(n,n_v,y,X,Pmatrix_MCCV,C2,b,train.MCCV) != 0)
  #BICV[i] <- sum( CV(n_v,y,X, BICV = BIBD) !=0)
}

##Computes the probability for a Criterion to Chooses a Cat II Model of a given size.
Probabilities <- matrix(0L,4,p-p.True+1)
for(i in 0:(p-p.True)){
  Probabilities[,i+1] <- c(sum( BIC == (p.True + i)), sum( AIC == (p.True + i)),  sum( CV == (p.True + i)) , sum( MCV == (p.True + i)) )/m
}
Probabilities






