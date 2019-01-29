##Data for True Regression Model
#Intercept
x_1 <- rep(1,n)         
x_2 <- rnorm(n,2,1)
x_3 <- rnorm(n,0,1)
x_4 <- rnorm(n,5,2)
x_5 <- rnorm(n,4,1)

##Possible values for Beta and X
beta <- c(1.5,3,0,0,0)
X <- cbind(x_1,x_2,x_3,x_4,x_5)

eps <- rnorm(n,0,1)     #Errorterm
y <- X%*%beta+eps

##Number of Possible Regressors
p <- length(X[1,])                          

#Creats the set {1,...,p} from which we want to generate the Powerset
Index <- seq(1,p,1)    

#Denotes the number of Possible Models out of {1,...,p}
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

##Sample Splitting CV(1)
n <- length(X[,1])





##List Inverses CV(1) for fixed n

#Combination of Sample Splits
train <- combn(seq(1,n,1),n-1)

Inverses_CV1 <- list()
k <- 0
for (i in 1:length(A[1,])) {
  for (j in 1:length(train[1,])) {
    k <- k+1
    X.train <- X[train[,j],A[,i]] 
    Inverses_CV1[[paste0("element",k)]] <- solve(t(X.train)%*%X.train)
  }
}





##List Inverses MCCV for fixed n

##Combination of Samplesplits
train <- matrix(ncol = b, nrow = n-n_v)
for (i in 1:b) {
  train[,i] <- sample(seq(1,n,1),n-n_v,replace = FALSE)
}

Inverses_MCCV <- list()
k <- 0
for (i in 1:length(A[1,])) {
  for (j in 1:length(train[1,])) {
    k <- k+1
    X.train <- X[train[,j],A[,i]] 
    Inverses_MCCV[[paste0("element",k)]] <- solve(t(X.train)%*%X.train)
  }
}





##List Inveses AIC/BIC for fixed n

Inverses_IC <- list()
for (i in 1:length(A[1,])) {
  X.train <- X[,A[,i]]
  Inverses_IC[[paste0("element",i)]] <- solve(t(X.train)%*%X.train)
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
        
        ##Data for the choosen Model
        X.model <- X[,A[,j]] 
        
        ##Expectation observation j
        mu <- X.model%*%i%*%t(X.model)%*%y
        
        error <- sum((y-mu)^2)
        
        ##Varianz
        var <- 1/n*error
        
        ##Calculate Log Likelihood value
        LogL <- -n/2*(log(2*pi)-log(var)-1)
      
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
        
        ##Data for the choosen Model
        X.model <- X[,A[,j]] 
        
        ##Expectation observation j
        mu <- X.model%*%i%*%t(X.model)%*%y
        
        error <- sum((y-mu)^2)
        
        ##Varianz
        var <- 1/n*error
        
        ##Calculate Log Likelihood value
        LogL <- -n/2*(log(2*pi)-log(var)-1)
      
      InfoCriterion[j] <- log(n)*k - 2*LogL
    }
    }
  
  #Choose the Model which minimze the Information Criterion
  TheChosenOne <- which.min(InfoCriterion)
  return(A[,TheChosenOne])
}

Rprof("AIC")
InfoCrit(y,X,Inverses_IC,A)
Rprof(NULL)
summaryRprof("AIC")



##CV
CV <- function(n_v,y,X,I,A,train){
  #Compute Prediction Errors for all sets in A
  MeanPred.Error <- c()
  k <- 0
  for (i in 1:length(A[1,])) {
    Pred.Error <- c()
    for (j in 1:length(train[1,])) {
      k <- k+1
      #To make the Code a bit shorter
      train.j <- train[,j]                 
      X.train <- X[train.j,A[,i]]    
      #Prediction Error for a given alpha and a given subset
      for(v in I[k]){
        Pred.Error[j] <- norm(as.matrix(y[-train.j]-X[-train.j,A[,i]]%*%v%*%t(X.train)%*%y[train.j]),"2")^2
      }
    }
    MeanPred.Error[i] <- 1/length(train[1,])*sum(Pred.Error)
  }
  TheChosenOne <- which.min(MeanPred.Error)
  return(A[,TheChosenOne])
}

Rprof("CV")
CV(395,y,X,Inverses_MCCV,A,train)
Rprof(NULL)
summaryRprof("CV")

Rprof("CV2")
CV(n_v,y,X,MonteCarlo = 5*500)
Rprof(NULL)
summaryRprof("CV2")
