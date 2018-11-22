##Generate some data for regression
n <- 50


##Data for True Regression Model
p.True<- 4                   #Number of True regressors 
x_1 <- rep(1,n)         #Intercept
x_2 <- rnorm(n,2,1)
x_3 <- rnorm(n,0,1)
x_4 <- rnorm(n,5,2)

eps <- rnorm(n,0,1)     #Errorterm


##Generate some unecessary extra data 
x_5 <- rnorm(n,4,1)
x_6 <- rnorm(n,2,2)
x_7 <- rnorm(n,1,3)
x_8 <- rnorm(n,3,7)
x_9 <- rnorm(n,0,1)
x_10 <- rnorm(n,5,5)


beta.vec <- c(1.5,0.5,3,2,rep(0,6)) 
X <- cbind(x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,x_10)


##Genereta dependent variables acording to the Above statet Model
y <- X%*%beta.vec + eps    #True Model

##Number of Regressors
p <- length(X[1,])

##For the CV we need to compute the set A of all possible models 
Index <- seq(1,p,1)                         #Creats the set {1,...,p} from which we want to generate the Powerset

col.A <- 0                                  #Denotes the number of Possible Models out of {1,...,p} Regressors           
for (i in 1:p) {
  col.A <- col.A + choose(p,i)   
}

A <- matrix(0L,nrow = p, ncol = col.A)      #Denote A as Powerset of {1,...,p}
k <- choose(p,1)
l <- 1
for (i in 1:length(beta.vec)) {
  A[1:i,l:k] <- combn(Index,i)             #combn spits out all combinations of i elements in Index 
  k <-k + choose(p,i+1)
  l <- l + choose(p,i)
}


##Split A into two dijoint subsets, ie, the set of Category I Models and Category II Models 
##We may need this sets later n for Simulation study's
coln <- c()
for (i in p.True:col.A) {
  if(all(seq(1,p.True,1) %in% A[,i])){    
    coln <- c(coln,i)
  }
}
C1 <- A[,coln]           #Set of all Category I Models
C2 <-A[,-coln]           #Set of all Category II Models


##Leave one out Crossvalidation

CV <- function(n_v, y, X, Alpha = NULL, MonteCarlo = NULL, Replacement = FALSE ){  #n_v = #leaved out data points
  #n_v          Number of leaved out data points
  #y,X          Data for Regression
  #Alpha        Set of possible Modelvaritaions for which the CV should be calculated 
  #             (By Default use all possible Models)
  #MonteCarlo   Number of Subsets of {1,...,n} which is randomly drawn for a Monte Carlo CV
  #             (By Default do K-Fold CV)
  #Replacement  Replacement for Monte Carlo (Default False)
  
  ##Change some Variable names to keep the code shorter
  A <- Alpha
  b <- MonteCarlo
  
  ##Number of Possible Regressors
  p <- length(X[1,])                          
  
  ##Set of possible Models
  if(is.null(A)){
    
    Index <- seq(1,p,1)    #Creats the set {1,...,p} from which we want to generate the Powerset
    
    col.A <- 0                                  #Denotes the number of Possible Models out of {1,...,p}             
    for (i in 1:p) {
      col.A <- col.A + choose(p,i)   
    }
  
    A <- matrix(0L,nrow = p, ncol = col.A)      #Denote A as Powerset of {1,...,p}
    k <- choose(p,1)
    l <- 1
    for (i in 1:p) {
      A[1:i,l:k] <- combn(Index,i)              #combn spits out all combinations of i elements in Index 
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
  }else{
    
    ##For the general case with all subsets of {1,...,n} 
    train <- combn(seq(1,n,1),n-n_v)
  }
  
  
  #Compute Prediction Errors f??r all sets in A
  MeanPred.Error <- c()
  for (i in 1:length(A[1,])) {
    
    Pred.Error <- c()
    for (j in 1:length(train[1,])) {
      
      train.j <- train[,j]                #To make the Code a bit shorter 
      X.train <- X[train.j,A[,i]]         
      Pred.Error[j] <- norm(as.matrix(y[-train.j]-X[-train.j,A[,i]]%*%solve(t(X.train)%*%X.train)%*%t(X.train)%*%y[train.j]),"2")^2
    
    }
    MeanPred.Error[i] <- 1/length(train[1,])*sum(Pred.Error)
  }
  TheChosenOne <- which.min(MeanPred.Error)
  return(A[,TheChosenOne])
}

CV(1,y,X,MonteCarlo = 50)




##Balance Incomplete Set
v <- 7        #lenght of Superstet , i.e, number of Observations
b <- 7        #Number of Blocks, i.e, subsets
k <- 3        #Number of elements in the blocks
r <- 3        #Observation i apears in r different blocks
lambda <- 1   #The pair i,j apears in lamda diffrent blocks

A <- matrix(0L,ncol = v, nrow = b)    #matrix of all subsets 

##Defining the constraints a BIBD should fufill
##Counts number of 1 and 0 in each row (Zeile), i.e, how often i occours (counts r)
rho_1 <- function(A){
  v <- length(A[,1])
  rho_1 <- c()
  for (i in 1:v) {
    rho_1[i] <- sum(A[i,])
  }
  return(rho_1)
}

rho_0 <- function(A){
  v <- length(A[,1])
  b <- length(A[1,])
  rho_0 <- rep(b,v)-rho_1(A)
  return(rho_0)
}

##Count number of 1 and 0 in each column (Spale), i.e, how many elements are in each block (Counts 3)
xi_1 <- function(A){
  b <- length(A[1,])
  xi_1 <- c()
  for (i in 1:b) {
    xi_1[i] <- sum(A[,i])
  }
  return(xi_1)
}

xi_0 <- function(A){
  b <- length(A[1,])
  v <- length(A[,1])
  xi_0 <- rep(v,b)-xi_1(A)
  return(xi_0)
}

##Scalar Product between each two rows, i.e , counts of often pairs ocoure
pi_1 <- function(A){
  v <- length(A[,1])
  comb <- combn(v,2)
  pi_1 <- c()
  for (i in 1:length(comb[1,])) {
    comb_ij <- comb[,i]
    pi_1[i] <- A[comb_ij[1],]%*%A[comb_ij[2],]
  }
  return(pi_1)
}

pi_0 <- function(A){
  pi_0 <- rep(length(A[1,]),choose(length(A[,1]),2)) - pi_1(A)
  return(pi_0)
}

##Later on a BIBD should hold the following constraints
# rho_1 <= r      rho_0 <= b-r
# xi_1 <= k       xi_0 <= v-k
# pi_1 <= lambda  pi_0 <= b-lambda

##Constructing BIBD by forward checking the constraints

A.vec <- c(A)

##Run the following algorithm as long as all constraints are satisfied 
while(
  !(all(rho_1(A) <= r) & all(rho_0(A) <= b-r) & 
  all(xi_1(A) <= k) & all(xi_0(A) <= v-k) & 
  all(pi_1(A) <= lambda) & all(pi_0(A) <= b-lambda))
  ){
  
  ##Change randome one entry in A
  m <- sample(length(A.vec),1)
  #Either Change a 1 to 0 
  if(A.vec[m] == 1){
    A.vec[m] <- 0
    A <- matrix(A.vec, ncol = v, nrow = b)
    #check wether or not the change will hurt a constraint. If this happens undo he change
    if(!all(!rho_0(A)>b-r) | !all(!xi_0(A)>v-k) | !all(!pi_0(A)>b-lambda)){
      A.vec[m] <- 1
      A <- matrix(A.vec, ncol = v, nrow = b)
    }else{}
  #or change a 0 to 1
  }else{
    A.vec[m] <- 1
    A <- matrix(A.vec, ncol = v, nrow = b)
    #check wether or not the change will hurt a constraint. If this happens undo he change
    if(!all(!rho_1(A)>r) | !all(!xi_1(A)>k) | !all(!pi_1(A)>lambda)){
      A.vec[m] <- 0
      A <- matrix(A.vec, ncol = v, nrow = b)
    }else{}
  }
  A <- matrix(A.vec, ncol = v, nrow = b)
}










