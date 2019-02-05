#faster AIC and BIC 

#We try to fasten our earlier AIC and BIC selection by
# - outsourcing the OLS estimation (this shouldn't be that much faster)
# - using algebra to create a leaner criterion function

#generate the data and further inputs
n <- 500
p <- 5
beta <- c(1.5,3,0,0,0)
X <- DataGen_X(n)
y <- X%*%beta + rnorm(n,0,1)
gram <- t(X)%*%X

#Generate the set of all possible Alphas

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

InfoCrit_fast <- function(y,X,B, Alpha = A, Criterion = "AIC"){
  #this function chooses a model alpha by AIC or BIC respectively (under the assumption 
  # of a Gaussian error term)
  #
  #y            data for the dependent variable
  #X            data fot the regressors
  #B            the empirical gram matrix of the regressors
  #Alpha        Set of possible Modelvaritaions for which the CV should be calculated 
  #             (By Default use all possible Models)
  #Criterion    Decides if we want to use the Akaike or Bayesian information criterion
  
  ##Change some Variable names to keep the code shorter
  A <- Alpha
  
  ##Number of Observations
  n <- length(y)
  
  ##Set of possible Models
  if(is.null(A)){
    
    ##Number of Possible Regressors
    p <- length(X[1,])
    
    ##Creats the set {1,...,p} from which we want to generate the Powerset
    Index <- seq(1,p,1)    
    
    ##Denotes the number of Possible Models out of {1,...,p} 
    col.A <- 2^p-1     
    
    ##Denote A as Powerset of {1,...,p}
    A <- matrix(0L,nrow = p, ncol = col.A)      
    k <- choose(p,1)
    l <- 1
    for (i in 1:p) {
      ##combn spits out all combinations of i elements in Index
      A[1:i,l:k] <- combn(Index,i)               
      k <-k + choose(p,i+1)
      l <- l + choose(p,i)
    }
  }
  
  col.A <- length(A[1,])
  
  
  ##Vector of Likelihood values for different Models
  InfoCriterion <- numeric(col.A)
  if(Criterion == "AIC"){
    for (i in 1:length(A[1,])) {
      
      ##Number of regressors
      k <- sum(A[,i] != 0)
      
      ##Data for the choosen Model
      X.model <- X[,A[,i]] 
      
      ##Expectation observation i
      mu <- X.model%*%solve(B[A[,i],A[,i]])%*%t(X.model)%*%y
      
      error <- sum((y-mu)^2)/n
      
      ##Calculate Log Likelihood value
      LogL <- -n/2*(log(2*pi)-log(error)-1)
      
      ##Calculate Information Criterion
      InfoCriterion[i] <- 2*k-2*LogL
      
    }
  }
  
  if(Criterion == "BIC"){
    for (i in 1:length(A[1,])) {
      
      ##Number of regressors
      k <- sum(A[,i] != 0)
      
      ##Data for the choosen Model
      X.model <- X[,A[,i]] 
      
      ##Expectation observation i
      mu <- X.model%*%solve(t(X.model)%*%X.model)%*%t(X.model)%*%y
      
      error <- sum((y-mu)^2)/n
      
      ##Calculate Log Likelihood value
      LogL <- -n/2*(log(2*pi)-log(error)-1)
      
      ##Calculate Information Criterion
      InfoCriterion[i] <- log(n)*k - 2*LogL
    }
  }
  
  #Choose the Model which minimze the Information Criterion
  TheChosenOne <- which.min(InfoCriterion)
  return(A[,TheChosenOne])
}

InfoCrit <- function(y,X,Alpha = NULL, Criterion = "AIC"){
  #Alpha        Set of possible Modelvaritaions for which the CV should be calculated 
  #             (By Default use all possible Models)
  #Criterion    Decides if we want to use the Akaike or Bayesian information criterion
  
  ##Change some Variable names to keep the code shorter
  A <- Alpha
  
  ##Number of Observations
  n <- length(y)
  
  ##Set of possible Models
  if(is.null(A)){
    
    ##Number of Possible Regressors
    p <- length(X[1,])
    
    ##Creats the set {1,...,p} from which we want to generate the Powerset
    Index <- seq(1,p,1)    
    
    ##Denotes the number of Possible Models out of {1,...,p} 
    col.A <- 2^p-1     
    
    ##Denote A as Powerset of {1,...,p}
    A <- matrix(0L,nrow = p, ncol = col.A)      
    k <- choose(p,1)
    l <- 1
    for (i in 1:p) {
      ##combn spits out all combinations of i elements in Index
      A[1:i,l:k] <- combn(Index,i)               
      k <-k + choose(p,i+1)
      l <- l + choose(p,i)
    }
  }
  
  ##Vector of Likelihood values for different Models
  InfoCriterion <- c()
  for (i in 1:length(A[1,])) {
    
    ##Number of regressors
    k <- sum(A[,i] != 0)
    
    ##Data for the choosen Model
    X.model <- X[,A[,i]] 
    
    ##Expectation observation i
    mu <- X.model%*%solve(t(X.model)%*%X.model)%*%t(X.model)%*%y
    
    error <- sum((y-mu)^2)
    
    ##Varianz
    var <- 1/n*error
    
    ##Calculate Log Likelihood value
    LogL <- -n/2*log(2*pi)-n/2*log(var)-1/(2*var)*error
    
    ##Calculate Information Criterion
    
    ##For AIC
    if(Criterion == "AIC"){
      InfoCriterion[i] <- 2*k-2*LogL
    }
    
    ##For BIC
    if(Criterion == "BIC"){
      InfoCriterion[i] <- log(n)*k - 2*LogL
    }
    
  }
  
  #Choose the Model which minimze the Information Criterion
  TheChosenOne <- which.min(InfoCriterion)
  return(A[,TheChosenOne])
}

#-------------------------------------------------------------------------------
#Inpsecting whether there are speed improvements
#-------------------------------------------------------------------------------

Rprof("faster")
y <- X%*%beta + rnorm(n,0,1)
InfoCrit_fast(y,X,gram,A)
InfoCrit_fast(y,X,gram,A,"BIC")
Rprof(NULL)
summaryRprof("faster")
#$`by.self`
#self.time self.pct total.time total.pct
#"%*%"       0.1      100        0.1       100
#
#$by.total
#total.time total.pct self.time self.pct
#"%*%"                  0.1       100       0.1      100
#"InfoCrit_fast"        0.1       100       0.0        0
#
#$sample.interval
#[1] 0.02
#
#$sampling.time
#[1] 0.1

Rprof("incumbent")
y <- X%*%beta + rnorm(n,0,1)
InfoCrit(y,X,A)
InfoCrit(y,X,A,"BIC")
Rprof(NULL)
summaryRprof("incumbent")
#$`by.self`
#self.time self.pct total.time total.pct
#"%*%"                  0.04       50       0.04        50
#"as.list.default"      0.02       25       0.02        25
#"getInlineInfo"        0.02       25       0.02        25
#
#$by.total
#total.time total.pct self.time self.pct
#"InfoCrit"                   0.08       100      0.00        0
#"%*%"                        0.04        50      0.04       50
#"cmpfun"                     0.04        50      0.00        0
#"compiler:::tryCmpfun"       0.04        50      0.00        0
#"doTryCatch"                 0.04        50      0.00        0
#"tryCatch"                   0.04        50      0.00        0
#"tryCatchList"               0.04        50      0.00        0
#"tryCatchOne"                0.04        50      0.00        0
#"as.list.default"            0.02        25      0.02       25
#"getInlineInfo"              0.02        25      0.02       25
#"as.list"                    0.02        25      0.00        0
#"cmp"                        0.02        25      0.00        0
#"cmpBuiltinArgs"             0.02        25      0.00        0
#"cmpCall"                    0.02        25      0.00        0
#"cmpForBody"                 0.02        25      0.00        0
#"cmpPrim2"                   0.02        25      0.00        0
#"cmpSymbolAssign"            0.02        25      0.00        0
#"findLocalsList"             0.02        25      0.00        0
#"findLocalsList1"            0.02        25      0.00        0
#"FUN"                        0.02        25      0.00        0
#"funEnv"                     0.02        25      0.00        0
#"genCode"                    0.02        25      0.00        0
#"h"                          0.02        25      0.00        0
#"lapply"                     0.02        25      0.00        0
#"make.functionContext"       0.02        25      0.00        0
#"tryInline"                  0.02        25      0.00        0

$sample.interval
[1] 0.02

$sampling.time
[1] 0.08