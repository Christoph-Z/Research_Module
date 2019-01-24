set.seed(100)
library(ibd)
##Generate some data for regression
n <- 50

##Function that Generates Data
#The only purpose of this function is to save some lines of code. Hence that the function
#is not very general
DataGen <- function(n,p,p.True){
  ##Data for True Regression Model
  #Intercept
  x_1 <- rep(1,n)         
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
  
  ##Possible values for Beta and X
  beta.pos <- c(1.5,3,2,5,3,7,4,6,5,7)
  X.pos <- cbind(x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,x_10)
  
  beta.vec <- c(beta.pos[1:p.True],rep(0,p-p.True)) 
  X<- X.pos[,1:p]
  
  
  ##Genereta dependent variables acording to the Above statet Model
  y <- X%*%beta.vec + eps    #True Model
  return(cbind(y,X))
}

#Number of True regressors 
p.True<- 2

##Number of Regressors
p <- 5

#Data
Data <- DataGen(n,p,p.True)
y <- Data[,1]
X <- Data[,-1]

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

##Leave one out Crossvalidation

#Spits out a vector of Integers, i.e, alpha
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



##Akaike and Schwarz Information Criterion
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


####################SIMULATION STUDIES#####################

##Simulation I
##Simulation of P(M_CV is in Cat II) 

##Number of True regressors 
p.True<- 2

##Number of Regressors
p <- 5

##Vector of diffrent Sample sizes
N <- seq(15,300,2)

##Number of Repetitions of the Experiment
m <- 500

##Vector of probabilties of frequencys of choosing a Cat II Model for each samplezize by m repetations 
BIC <- rep(0,length(N))
AIC <- rep(0,length(N))
CV1 <- rep(0,length(N))
MCV <- rep(0,length(N))

for (i in 1:length(N)) {
  n <- N[i]
  n_v <- floor(n-n^(3/4))
  
  ##Number of Samplepartions
  b <- 5*n
  
  for (j in 1:m) {
    ##Generating Data
    Data <- DataGen(n,p,p.True)
    y <- Data[,1]
    X <- Data[,-1]
    
    #Define some temporary Variables, i.e, define for a given Method the in iteration i choosed Model
    temp.BIC <- InfoCrit(y,X,Criterion = "BIC")
    temp.AIC <- InfoCrit(y,X) 
    temp.CV1 <- CV(1,y,X)
    temp.MCV <- CV(n_v,y,X,MonteCarlo = b)
    
    ##Counting how many times a Cat II Model is picked for n= N[i] in m interations
    BIC[i] <- BIC[i] + (sum( temp.BIC < (p.True+1) & temp.BIC >0) == p.True)
    AIC[i] <- AIC[i] + (sum( temp.AIC < (p.True+1) & temp.AIC > 0) == p.True)
    CV1[i] <- CV1[i] + (sum( temp.CV1 < (p.True+1) & temp.CV1 > 0) == p.True)
    MCV[i] <- MCV[i] + (sum( temp.MCV < (p.True+1) & temp.MCV > 0) == p.True)
  }
  
}


##Simulation II 
##based on the Simulations in Shao

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
CV1 <- c()
MCV <- c()
BICV <- c()

##Number of leaved out data and number of samplepartions
n_v <- 394
b <- 5*n

##Calculate BICD for BICV with ibd package
BIBD <- ibd(n,b,n-n_v,pbar = TRUE)$N

 for (i in 1:m) {
  ##Generating Data
  Data <- DataGen(n,p,p.True)
  y <- Data[,1]
  X <- Data[,-1]
  
  BIC[i] <- sum( InfoCrit(y,X,C2,Criterion = "BIC") != 0)
  AIC[i] <- sum( InfoCrit(y,X,C2,Criterion = "AIC") != 0)
  CV1[i] <- sum( CV(1,y,X,C2) != 0)
  MCV[i] <- sum( CV(n_v,y,X, MonteCarlo = b) != 0)
  #BICV[i] <- sum( CV(n_v,y,X, BICV = BIBD) !=0)
}

##Computes the probability for a Criterion to Chooses a Cat II Model of a given size.
Probabilities <- matrix(0L,4,p-p.True+1)
for(i in 0:(p-p.True)){
  Probabilities[,i+1] <- c(sum( BIC == (p.True + i)), sum( AIC == (p.True + i)),  sum( CV1 == (p.True + i)) , sum( MCV == (p.True + i)) )/m
}
Probabilities


##Simulation III 
##Wahl n_v

##Vector of different Sample sizes
N <- seq(15,300,2)

##Number of True regressors 
p.True<- 2

##Number of Regressors
p <- 5

##Number of Iteration of Modelfitting
m <- 1000

##Criterion Function
##Is minimal iff CV chooses the true model and prefers less conservative Choices
L <- function(model){
   + sum( model > p.True) -sum( (model <= p.True & model > 0) )
}

##Vector of optimal n_v given n
M <- rep(0,length(N))

##Simulation
for (i in 1:length(N)) {
  n <- N[i]
  ##Define the number b of Samplepartions as function of n
  ##Shao 93 claimed that b=O(n)
  b <- n*10
    
  ##Grid of possible n_v values for a given Samplezize n in N
  N_v <- seq(2,n-p-1,2)
    
  ##Vector of all Criterionfunction values given n_v
  N_v.CritValue <- rep(0,length(N_v))
  for (k in 1:length(N_v)) {
    n_v <- N_v[k]
      
    ##For statistical Influence fit the Model m times
    for (l in 1:m) {
        
      ##Generate Data
      Data <- DataGen(n,p,p.True)
      y <- Data[,1]
      X <- Data[,-1]
      
      ##Select a Model
      model <- CV(n_v,y,X, MonteCarlo = b)
      N_v.CritValue[k] <-  N_v.CritValue[k] + L(model)
    }
  TheChoosenOne <- which.min(N_v.CritValue)
  M[i] <- N_v[TheChoosenOne]
  }
 }







