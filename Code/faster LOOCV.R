#A faster LOOCV function

#Since in our simulation design the regressors are fixed, some parts of the OLS estimate 
#doesn't change in the different simulation steps and hence could be pulled out of the
#CV function to reduce redundant computations.
#
#In each new simulation step we only draw new error terms and hence only the y values differ.
#Thus we can give the matrix inverses and the set of Alphas directly to the function instead
#recomputing it over and over.
#
#Another point where we want to save computing power is by using an algebraic reformulation 
#of the average squared prediction error as in Shao(93).

#---------------------------------------------------------------------------------------
#Generate the Data
n <- 500
p <- 5
p.True <- 2
beta <- c(1.5,3,0,0,0)
DataGen_X <- function(n){
  x_1 <- rep(1,n)         
  x_2 <- rnorm(n,2,1)
  x_3 <- rnorm(n,0,1)
  x_4 <- rnorm(n,5,2)
  x_5 <- rnorm(n,4,1)
  return(X<- cbind(x_1,x_2,x_3,x_4,x_5))
}
X <- DataGen_X(n)
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
#generate the diagonal elements of the projection matrices corresponding to model alpha
proj_diag <- matrix(0L,nrow = n, ncol = col.A)
for(j in 1:col.A){
  for(i in 1:n){
    alpha <- A[A[,j]!=0,j]
    proj_diag[i,j] <- t(X[i,alpha])%*%solve(gram[alpha,alpha])%*%X[i,alpha]
  }
}

y <- X%*%beta + rnorm(n,0,1)
#b <- matrix(0L, nrow = p, ncol = col.A)
#for(j in 1:col.A){
#  b[A[A[,j]!=0,j],j]<- solve(gram[A[A[,j]!=0,j],A[A[,j]!=0,j]])%*%t(X[,A[A[,j]!=0,j]])%*%y
#}
#MeanPred.Error <- numeric(col.A)
#for (i in 1:col.A) {
#  alpha <- A[A[,i]!=0,i]
#  if(length(alpha)==1){
#    MeanPred.Error[i] <- t((rep(1,n)-proj_diag[,i])^(-2))%*%((y - X[,A[A[,i]!=0,i]]*b[A[A[,i]!=0,i],i])^2)/n
#  }else{
#    MeanPred.Error[i] <- t((rep(1,n)-proj_diag[,i])^(-2))%*%((y - X[,A[A[,i]!=0,i]]%*%b[A[A[,i]!=0,i],i])^2)/n
#  }
#}
#A[,which.min(MeanPred.Error)]

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

CV <- function(n_v, y, X, Alpha = NULL, MonteCarlo = NULL, Replacement = FALSE, BICV = NULL ){  #n_v = #leaved out data points
  #n_v          Number of leaved out data points
  #y,X          Data for Regression
  #Alpha        Set of possible Modelvaritaions for which the CV should be calculated 
  #             (By Default use all possible Models)
  #MonteCarlo   Number of Subsets of {1,...,n} which is randomly drawn for a Monte Carlo CV
  #             (By Default do K-Fold CV)
  #Replacement  Replacement for Monte Carlo (Default False)
  #
  ##BICV        Incidence Matrix for a BICV
  
  ##Change some Variable names to keep the code shorter
  A <- Alpha
  b <- MonteCarlo
  
  ##Number of Possible Regressors
  p <- length(X[1,])                          
  
  ##Set of possible Models
  if(is.null(A)){
    
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
  }
  
  
  ##Number of Observations
  n <- length(y)
  
  
  ##Combinations of Sample partions for fitting the model
  
  ##For the Mone Carlo CV with b subsets of {1,...,n}
  if(!is.null(b)){
    train <- matrix(ncol = b, nrow = n-n_v)
    for (i in 1:b) {
      train[,i] <- sample(seq(1,n,1),n-n_v,replace = Replacement)
    }
  }else if(!is.null(BICV)){
    #Convert the Incedence matrix of BIBD in our Notation
    train <- BICV * seq(1,n,1)
  }else{
    ##For the general case with all subsets of {1,...,n} 
    train <- combn(seq(1,n,1),n-n_v)
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
      Pred.Error[j] <- norm(as.matrix(y[-train.j]-X[-train.j,A[,i]]%*%solve(t(X.train)%*%X.train)%*%t(X.train)%*%y[train.j]),"2")^2
      
    }
    MeanPred.Error[i] <- 1/length(train[1,])*sum(Pred.Error)
  }
  TheChosenOne <- which.min(MeanPred.Error)
  return(A[,TheChosenOne])
}

#test for speed
system.time(LOOCV(y,X,gram,proj_diag,A))
#user  system elapsed 
#0       0       0 
Rprof("LOOCV")
y <- X%*%beta + rnorm(n,0,1)
LOOCV(y,X,gram,proj_diag,A)
Rprof(NULL)
summaryRprof("LOOCV")
#$`by.self`
#self.time self.pct total.time total.pct
#"tryCatchOne"      0.02      100       0.02       100
#
#$by.total
#total.time total.pct self.time self.pct
#"tryCatchOne"               0.02       100      0.02      100
#".rs.callAs"                0.02       100      0.00        0
#"Rprof"                     0.02       100      0.00        0
#"tryCatch"                  0.02       100      0.00        0
#"tryCatchList"              0.02       100      0.00        0
#"withCallingHandlers"       0.02       100      0.00        0
#
#$sample.interval
#[1] 0.02
#
#$sampling.time
#[1] 0.02

#Compare this with the other CV method
Rprof("CV")
y <- X%*%beta + rnorm(n,0,1)
CV(1,y,X,A)
Rprof(NULL)
summaryRprof("CV")
#$`by.self`
#self.time self.pct total.time total.pct
#"as.matrix"                 0.50    31.65       1.16     73.42
#"%*%"                       0.18    11.39       0.18     11.39
#"t"                         0.16    10.13       0.20     12.66
#"CV"                        0.12     7.59       1.56     98.73
#"La.svd"                    0.12     7.59       0.22     13.92
#"solve"                     0.10     6.33       0.48     30.38
#"diag"                      0.08     5.06       0.08      5.06
#"norm"                      0.04     2.53       1.44     91.14
#"svd"                       0.04     2.53       1.40     88.61
#"solve.default"             0.04     2.53       0.20     12.66
#"colnames<-"                0.04     2.53       0.04      2.53
#"t.default"                 0.04     2.53       0.04      2.53
#".rs.enqueClientEvent"      0.02     1.27       0.02      1.27
#"is.atomic"                 0.02     1.27       0.02      1.27
#"is.finite"                 0.02     1.27       0.02      1.27
#"min"                       0.02     1.27       0.02      1.27
#"ncol"                      0.02     1.27       0.02      1.27
#"rownames"                  0.02     1.27       0.02      1.27
#
#$by.total
#total.time total.pct self.time self.pct
#"CV"                         1.56     98.73      0.12     7.59
#"norm"                       1.44     91.14      0.04     2.53
#"svd"                        1.40     88.61      0.04     2.53
#"as.matrix"                  1.16     73.42      0.50    31.65
#"solve"                      0.48     30.38      0.10     6.33
#"La.svd"                     0.22     13.92      0.12     7.59
#"t"                          0.20     12.66      0.16    10.13
#"solve.default"              0.20     12.66      0.04     2.53
#"%*%"                        0.18     11.39      0.18    11.39
#"diag"                       0.08      5.06      0.08     5.06
#"colnames<-"                 0.04      2.53      0.04     2.53
#"t.default"                  0.04      2.53      0.04     2.53
#".rs.enqueClientEvent"       0.02      1.27      0.02     1.27
#"is.atomic"                  0.02      1.27      0.02     1.27
#"is.finite"                  0.02      1.27      0.02     1.27
#"min"                        0.02      1.27      0.02     1.27
#"ncol"                       0.02      1.27      0.02     1.27
#"rownames"                   0.02      1.27      0.02     1.27
#"double"                     0.02      1.27      0.00     0.00
#"hook"                       0.02      1.27      0.00     0.00
#"matrix"                     0.02      1.27      0.00     0.00
#"Rprof"                      0.02      1.27      0.00     0.00
#
#$sample.interval
#[1] 0.02
#
#$sampling.time
#[1] 1.58

