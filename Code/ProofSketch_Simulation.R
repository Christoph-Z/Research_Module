n <- 500
p <- 5
p.True <- 2
beta <- c(4,7,0,0,0)
DataGen_X <- function(n){
  x_1 <- rep(1,n)         
  x_2 <- rnorm(n,0,1)
  x_3 <- rnorm(n,0,1)
  x_4 <- rnorm(n,0,1)
  x_5 <- rnorm(n,4,1)
  return(X<- cbind(x_1,x_2,x_3,x_4,x_5))
}
X <- DataGen_X(n)
sigma <- 1
y <- X%*%beta+rnorm(n,0,sigma)

##--------------------------------------------------------------------
##ESPE for CV(1) and CV(n_v) for Models in Category II

ESPE_CV1 <- function(d,n){
  return((1+d/n)*sigma^2)
}
ESP_CVn.v <- function(d,n){
  return(1+d/n^(3/4))
}

##--------------------------------------------------------------------
##ESPE for CV(1) and CV(n_v) for Models in Category II

CV <- function(n,y,X,Method = "CV1"){
  if(Method =="CV1"){
    n_v <- 1
  }else{n_v <- n-floor(n^{3/4})}
  
  train <- combn(seq(1,n,1),n-n_v)
  Pred.Error <- c()
  for (i in 1:length(train[1,])) {
    #To make the Code a bit shorter
    train.i <- train[,i]                 
    X.train <- X[train.i,]    
    #Prediction Error for a given alpha and a given subset
    Error <- y[-train.i]-X[-train.i,]%*%solve(t(X.train)%*%X.train)%*%t(X.train)%*%y[train.i]
    Pred.Error[i] <- t(Error)%*%Error
  }
  MeanPred.Error <- 1/n*sum(Pred.Error)
  return(MeanPred.Error)
}


N.grid <-seq(1,50,1)

CV1_d2 <- c()
for (i in 1:length(N.grid)) {
  CV1_d2[i] <- CV(N.grid[i],y[1:N.grid[i]],X[1:N.grid[i],1:2])
}

CV1_d3 <- c()
for (i in 1:length(N.grid)) {
  CV1_d3[i] <- CV(N.grid[i],y[1:N.grid[i]],X[1:N.grid[i],1:3])
}

d <- data.frame("Samplesize"=N.grid,"ESPE1_2"=ESPE_CV1(2,N.grid),"ESPE1_3"=ESPE_CV1(3,N.grid),"ASPE1_2"=CV1_d2,"ASPE1_3"=CV1_d3)

ggplot(d,aes(x=Samplesize)) +
  geom_line(aes(y=ESPE1_2,color = "M*"),size =1,linetype=5) +
  geom_line(aes(y=ASPE1_2,color = "M*"),size =1.2) +
  geom_line(aes(y=ESPE1_3,color = "M*+1"),size =1,linetype=5) +
  geom_line(aes(y=ASPE1_3,color = "M*+1"),size =1.2) +
  ylab("Expected Squared Prediction Error") + ggtitle("Asymptotic Properties CV(1)") +
  theme_bw() +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20,hjust = 0.5,face="bold")) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(plot.title = element_text(margin = margin(t = 0, r = 0 , b = 20, l = 0))) +
  theme(legend.title=element_blank(),legend.text = element_text(size = 13)) +
  theme(legend.position = c(0.85,0.8))



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



MCCV1_d2 <- c()
for (i in 1:length(N.grid)) {
  MCCV1_d2[i] <- CV(N.grid[i]-floor(N.grid[i]^{3/4}),y[1:N.grid[i]],X[1:N.grid[i],1:2], MonteCarlo = 5*N.grid[i])
}

MCCV1_d3 <- c()
for (i in 1:length(N.grid)) {
  MCCV1_d3[i] <- CV(N.grid[i]-floor(N.grid[i]^{3/4}),y[1:N.grid[i]],X[1:N.grid[i],1:3],MonteCarlo = 5*N.grid[i])
}



d <- data.frame("Samplesize"=N.grid,"ESPE1_2"=ESPE_CVn.v(2,N.grid),"ESPE1_3"=ESPE_CVn.v(3,N.grid),"ASPE1_2"=MCCV1_d2,"ASPE1_3"=MCCV1_d3)

