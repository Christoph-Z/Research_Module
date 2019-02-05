#alternative MCCV

#just want to test whether the different approach to calculating the ASPE yields any speed 
#improvements

n <- 500
n_v <- n-floor(n^0.75)
b <- 5*n
p <- 5
beta <- c(1.5,3,0,0,0)
X <- DataGen_X(n)
y <- X%*%beta + rnorm(n,0,1)
gram <- t(X)%*%X

MCCV <- function(n_v, y, X, B, Alpha = NULL, MonteCarlo = NULL, Replacement = FALSE){  #n_v = #leaved out data points
  #n_v          Number of leaved out data points
  #y,X          Data for Regression
  #B            gram matrix
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
    train <- matrix(ncol = b, nrow = n-n_v)
    for (i in 1:b) {
      train[,i] <- sample(seq(1,n,1),n-n_v,replace = Replacement)
    }
  
    col.A <- length(A[1,])
  
  #Compute Prediction Errors for all sets in A
  MeanPred.Error <- numeric(col.A)
  for (i in 1:col.A) {
    Pred.Error <- numeric(b)
    ginv <- solve(B[A[,i],A[,i]])
    b_alpha <- ginv%*%X[,A[,i]]%*%y
    
    for (j in 1:length(train[1,])) {
      #To make the Code a bit shorter
      train.j <- train[,j]                 
      X.train <- X[train.j,A[,i]]    
      #Prediction Error for a given alpha and a given subset
      Pred.Error[j] <- sum((solve(diag(n_v)-X[-train.j,A[,i]]%*%ginv%*%t(X[-train.j,A[,i]]))%*%(y[-train.j]-X[-train.j,A[,i]]%*%b_alpha))^2)/n_v
      
    }
    MeanPred.Error[i] <- 1/length(train[1,])*sum(Pred.Error)
  }
  TheChosenOne <- which.min(MeanPred.Error)
  return(A[,TheChosenOne])
}

CV <- function(n_v, y, X, B, Alpha = NULL, MonteCarlo = NULL, Replacement = FALSE, BICV = NULL ){  #n_v = #leaved out data points
  #n_v          Number of leaved out data points
  #y,X          Data for Regression
  #B            gram matrix
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

#--------------------------------------------------------------------------------------
#testing performance
#--------------------------------------------------------------------------------------
Rprof("new")
MCCV(n_v,y,X,gram,A,MonteCarlo=b)
Rprof(NULL)
summaryRprof("new")
#$`by.self`
#self.time self.pct total.time total.pct
#"solve.default"          496.98    96.40     502.80     97.53
#"diag"                     8.40     1.63       8.68      1.68
#"%*%"                      5.48     1.06       5.48      1.06
#"solve"                    2.62     0.51     512.56     99.42
#"t"                        0.38     0.07       0.46      0.09
#"MCCV"                     0.34     0.07     515.56    100.00
#"as.matrix"                0.20     0.04       0.22      0.04
#"ls"                       0.18     0.03       0.18      0.03
#"is.data.frame"            0.16     0.03       0.24      0.05
#"list.dirs"                0.10     0.02       0.10      0.02
#"rownames"                 0.08     0.02       0.16      0.03
#"nrow"                     0.08     0.02       0.08      0.02
#"sum"                      0.08     0.02       0.08      0.02
#"NROW"                     0.06     0.01       0.08      0.02
#"t.default"                0.06     0.01       0.06      0.01
#"colnames<-"               0.04     0.01       0.28      0.05
#"seq.default"              0.04     0.01       0.06      0.01
#"tryInline"                0.02     0.00       0.12      0.02
#"parent.env"               0.02     0.00       0.08      0.02
#"seq"                      0.02     0.00       0.08      0.02
#"$"                        0.02     0.00       0.02      0.00
#".Call"                    0.02     0.00       0.02      0.00
#".rs.makeCompletions"      0.02     0.00       0.02      0.00
#"as.matrix.default"        0.02     0.00       0.02      0.00
#"c"                        0.02     0.00       0.02      0.00
#"is.array"                 0.02     0.00       0.02      0.00
#"is.matrix"                0.02     0.00       0.02      0.00
#"isBlankLineRd"            0.02     0.00       0.02      0.00
#"pmin"                     0.02     0.00       0.02      0.00
#"psub"                     0.02     0.00       0.02      0.00
#"unique.default"           0.02     0.00       0.02      0.00
#
#$by.total
#total.time total.pct self.time self.pct
#"MCCV"                             515.56    100.00      0.34     0.07
#"solve"                            512.56     99.42      2.62     0.51
#"solve.default"                    502.80     97.53    496.98    96.40
#"diag"                               8.68      1.68      8.40     1.63
#"%*%"                                5.48      1.06      5.48     1.06
#"t"                                  0.46      0.09      0.38     0.07
#"<Anonymous>"                        0.46      0.09      0.00     0.00
#"colnames<-"                         0.28      0.05      0.04     0.01
#"FUN"                                0.28      0.05      0.00     0.00
#"lapply"                             0.28      0.05      0.00     0.00
#"Reduce"                             0.28      0.05      0.00     0.00
#"is.data.frame"                      0.24      0.05      0.16     0.03
#"as.matrix"                          0.22      0.04      0.20     0.04
#"ls"                                 0.18      0.03      0.18     0.03
#".rs.getCompletionsSearchPath"       0.18      0.03      0.00     0.00
#".rs.objectsOnSearchPath"            0.18      0.03      0.00     0.00
#"rownames"                           0.16      0.03      0.08     0.02
#"tryInline"                          0.12      0.02      0.02     0.00
#".rs.getCompletionsFile"             0.12      0.02      0.00     0.00
#"cmp"                                0.12      0.02      0.00     0.00
#"cmpCall"                            0.12      0.02      0.00     0.00
#"cmpfun"                             0.12      0.02      0.00     0.00
#"compiler:::tryCmpfun"               0.12      0.02      0.00     0.00
#"doTryCatch"                         0.12      0.02      0.00     0.00
#"genCode"                            0.12      0.02      0.00     0.00
#"h"                                  0.12      0.02      0.00     0.00
#"tryCatch"                           0.12      0.02      0.00     0.00
#"tryCatchList"                       0.12      0.02      0.00     0.00
#"tryCatchOne"                        0.12      0.02      0.00     0.00
#"list.dirs"                          0.10      0.02      0.10     0.02
#".rs.getCompletionsPackages"         0.10      0.02      0.00     0.00
#".rs.listDirs"                       0.10      0.02      0.00     0.00
#"cmpSymbolAssign"                    0.10      0.02      0.00     0.00
#"nrow"                               0.08      0.02      0.08     0.02
#"sum"                                0.08      0.02      0.08     0.02
#"NROW"                               0.08      0.02      0.06     0.01
#"parent.env"                         0.08      0.02      0.02     0.00
#"seq"                                0.08      0.02      0.02     0.00
#"findCenvVar"                        0.08      0.02      0.00     0.00
#"findFunDef"                         0.08      0.02      0.00     0.00
#"sample"                             0.08      0.02      0.00     0.00
#"suppressWarnings"                   0.08      0.02      0.00     0.00
#"withCallingHandlers"                0.08      0.02      0.00     0.00
#"t.default"                          0.06      0.01      0.06     0.01
#"seq.default"                        0.06      0.01      0.04     0.01
#".rs.getHelp"                        0.06      0.01      0.00     0.00
#".rs.getHelpFromObject"              0.06      0.01      0.00     0.00
#".rs.getHelpFunction"                0.06      0.01      0.00     0.00
#"Rd2HTML"                            0.06      0.01      0.00     0.00
#"tools:::httpd"                      0.06      0.01      0.00     0.00
#"writeBlock"                         0.06      0.01      0.00     0.00
#"writeContent"                       0.06      0.01      0.00     0.00
#"writeSection"                       0.06      0.01      0.00     0.00
#"addParaBreaks"                      0.04      0.01      0.00     0.00
#"of1"                                0.04      0.01      0.00     0.00
#"unique"                             0.04      0.01      0.00     0.00
#"writeLines"                         0.04      0.01      0.00     0.00
#"writeLinesUTF8"                     0.04      0.01      0.00     0.00
#"$"                                  0.02      0.00      0.02     0.00
#".Call"                              0.02      0.00      0.02     0.00
#".rs.makeCompletions"                0.02      0.00      0.02     0.00
#"as.matrix.default"                  0.02      0.00      0.02     0.00
#"c"                                  0.02      0.00      0.02     0.00
#"is.array"                           0.02      0.00      0.02     0.00
#"is.matrix"                          0.02      0.00      0.02     0.00
#"isBlankLineRd"                      0.02      0.00      0.02     0.00
#"pmin"                               0.02      0.00      0.02     0.00
#"psub"                               0.02      0.00      0.02     0.00
#"unique.default"                     0.02      0.00      0.02     0.00
#".rs.appendCompletions"              0.02      0.00      0.00     0.00
#".rs.emptyCompletions"               0.02      0.00      0.00     0.00
#".rs.getCompletionsArgument"         0.02      0.00      0.00     0.00
#".rs.getCompletionsFunction"         0.02      0.00      0.00     0.00
#".rs.getKnitParamsForDocument"       0.02      0.00      0.00     0.00
#".rs.getRCompletions"                0.02      0.00      0.00     0.00
#".rs.injectKnitrParamsObject"        0.02      0.00      0.00     0.00
#".rs.listFilesFuzzy"                 0.02      0.00      0.00     0.00
#".rs.normalizePath"                  0.02      0.00      0.00     0.00
#".signalSimpleWarning"               0.02      0.00      0.00     0.00
#"as.vector"                          0.02      0.00      0.00     0.00
#"cb$putconst"                        0.02      0.00      0.00     0.00
#"cmpBuiltinArgs"                     0.02      0.00      0.00     0.00
#"cmpCallArgs"                        0.02      0.00      0.00     0.00
#"cmpCallSymFun"                      0.02      0.00      0.00     0.00
#"cmpPrim1"                           0.02      0.00      0.00     0.00
#"doWithOneRestart"                   0.02      0.00      0.00     0.00
#"make.argContext"                    0.02      0.00      0.00     0.00
#"RdTags"                             0.02      0.00      0.00     0.00
#"sapply"                             0.02      0.00      0.00     0.00
#"simpleWarning"                      0.02      0.00      0.00     0.00
#"simplify2array"                     0.02      0.00      0.00     0.00
#"sort"                               0.02      0.00      0.00     0.00
#"structure"                          0.02      0.00      0.00     0.00
#"union"                              0.02      0.00      0.00     0.00
#"withOneRestart"                     0.02      0.00      0.00     0.00
#"withRestarts"                       0.02      0.00      0.00     0.00
#"writeWrapped"                       0.02      0.00      0.00     0.00
#
#$sample.interval
#[1] 0.02
#
#$sampling.time
#[1] 515.56


Rprof("old")
CV(n_v,y,X,Alpha = A,MonteCarlo = b)
Rprof(NULL)
summaryRprof("old")
#$`by.self`
#self.time self.pct total.time total.pct
#"%*%"                   12.58    60.66      12.58     60.66
#"as.matrix"              2.88    13.89      18.08     87.17
#"La.svd"                 0.84     4.05       1.74      8.39
#"solve.default"          0.60     2.89       1.56      7.52
#"t"                      0.50     2.41       0.74      3.57
#"CV"                     0.44     2.12      20.74    100.00
#"solve"                  0.44     2.12       2.66     12.83
#"svd"                    0.40     1.93      20.06     96.72
#"is.finite"              0.28     1.35       0.28      1.35
#"diag"                   0.26     1.25       0.30      1.45
#"colnames<-"             0.24     1.16       0.32      1.54
#"t.default"              0.24     1.16       0.24      1.16
#"matrix"                 0.22     1.06       0.22      1.06
#"norm"                   0.10     0.48      20.16     97.20
#"double"                 0.10     0.48       0.10      0.48
#"rownames"               0.10     0.48       0.10      0.48
#"as.matrix.default"      0.08     0.39       0.12      0.58
#"ncol"                   0.08     0.39       0.08      0.39
#"is.data.frame"          0.06     0.29       0.06      0.29
#"is.matrix"              0.04     0.19       0.04      0.19
#"nrow"                   0.04     0.19       0.04      0.19
#"seq.default"            0.02     0.10       0.08      0.39
#"pmin"                   0.02     0.10       0.06      0.29
#"FUN"                    0.02     0.10       0.04      0.19
#"vapply"                 0.02     0.10       0.04      0.19
#"%in%"                   0.02     0.10       0.02      0.10
#"any"                    0.02     0.10       0.02      0.10
#"c"                      0.02     0.10       0.02      0.10
#"dim"                    0.02     0.10       0.02      0.10
#"lapply"                 0.02     0.10       0.02      0.10
#"length"                 0.02     0.10       0.02      0.10
#"sample.int"             0.02     0.10       0.02      0.10
#
#$by.total
#total.time total.pct self.time self.pct
#"CV"                        20.74    100.00      0.44     2.12
#"norm"                      20.16     97.20      0.10     0.48
#"svd"                       20.06     96.72      0.40     1.93
#"as.matrix"                 18.08     87.17      2.88    13.89
#"%*%"                       12.58     60.66     12.58    60.66
#"solve"                      2.66     12.83      0.44     2.12
#"La.svd"                     1.74      8.39      0.84     4.05
#"solve.default"              1.56      7.52      0.60     2.89
#"t"                          0.74      3.57      0.50     2.41
#"colnames<-"                 0.32      1.54      0.24     1.16
#"diag"                       0.30      1.45      0.26     1.25
#"is.finite"                  0.28      1.35      0.28     1.35
#"t.default"                  0.24      1.16      0.24     1.16
#"matrix"                     0.22      1.06      0.22     1.06
#"as.matrix.default"          0.12      0.58      0.08     0.39
#"double"                     0.10      0.48      0.10     0.48
#"rownames"                   0.10      0.48      0.10     0.48
#"sample"                     0.10      0.48      0.00     0.00
#"ncol"                       0.08      0.39      0.08     0.39
#"seq.default"                0.08      0.39      0.02     0.10
#"seq"                        0.08      0.39      0.00     0.00
#"is.data.frame"              0.06      0.29      0.06     0.29
#"pmin"                       0.06      0.29      0.02     0.10
#"is.matrix"                  0.04      0.19      0.04     0.19
#"nrow"                       0.04      0.19      0.04     0.19
#"FUN"                        0.04      0.19      0.02     0.10
#"vapply"                     0.04      0.19      0.02     0.10
#"cmpfun"                     0.04      0.19      0.00     0.00
#"compiler:::tryCmpfun"       0.04      0.19      0.00     0.00
#"doTryCatch"                 0.04      0.19      0.00     0.00
#"tryCatch"                   0.04      0.19      0.00     0.00
#"tryCatchList"               0.04      0.19      0.00     0.00
#"tryCatchOne"                0.04      0.19      0.00     0.00
#"%in%"                       0.02      0.10      0.02     0.10
#"any"                        0.02      0.10      0.02     0.10
#"c"                          0.02      0.10      0.02     0.10
#"dim"                        0.02      0.10      0.02     0.10
#"lapply"                     0.02      0.10      0.02     0.10
#"length"                     0.02      0.10      0.02     0.10
#"sample.int"                 0.02      0.10      0.02     0.10
#"cb$putconst"                0.02      0.10      0.00     0.00
#"cmp"                        0.02      0.10      0.00     0.00
#"cmpBuiltinArgs"             0.02      0.10      0.00     0.00
#"cmpCall"                    0.02      0.10      0.00     0.00
#"cmpCallArgs"                0.02      0.10      0.00     0.00
#"cmpCallSymFun"              0.02      0.10      0.00     0.00
#"cmpComplexAssign"           0.02      0.10      0.00     0.00
#"cmpForBody"                 0.02      0.10      0.00     0.00
#"cmpPrim2"                   0.02      0.10      0.00     0.00
#"findLocalsList"             0.02      0.10      0.00     0.00
#"findLocalsList1"            0.02      0.10      0.00     0.00
#"funEnv"                     0.02      0.10      0.00     0.00
#"genCode"                    0.02      0.10      0.00     0.00
#"getInlineInfo"              0.02      0.10      0.00     0.00
#"h"                          0.02      0.10      0.00     0.00
#"make.functionContext"       0.02      0.10      0.00     0.00
#"tryInline"                  0.02      0.10      0.00     0.00
#
#$sample.interval
#[1] 0.02
#
#$sampling.time
#[1] 20.74