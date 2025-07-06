library(glmnet)
library(MASS)

#' This function provides change point estimation by Binary Segmentation Algorithm
#' for the high dimensional generalized linear regression model.
#' @export
#' @param X an n*p data matrix, n=observation, p=dimension
#' @param Y an n dimensional data vector, n=observation
#' @param kmax the maximum change point number.
#' @param c1 a constant specified by users.
#' @param m1 the sample size of the minimal interval.
#' @param delta the minimal interval length specified by users.
#' @return A list including: the estimated change point number, the estimated
#'   change point locations.
#' @importFrom stats coef
#' @importFrom glmnet glmnet
######## change point estimation by bsa ########
bsa_chgpt <- function(X, Y, kmax, c1, m1, delta) {
  print("Entered bsa_chgpt")
  n <- length(Y)
  p <- ncol(X)
  alpha <- c(1, n)
  shat  <- numeric()
  
  for (z in seq_len(kmax)) {
    previous_alpha <- alpha
    for (i in seq_len(length(previous_alpha) - 1)) {
      start_i <- previous_alpha[i]
      end_i   <- previous_alpha[i + 1]
      if ((start_i + m1) > (end_i - m1)) next
      
      losses <- rep(Inf, end_i)
      
      # full-segment fit
      len_full <- end_i - start_i + 1
      λ_full   <- c1 * ( sqrt(2 * log(2 * p) / len_full) + log(2 * p) / len_full )
      fit_full <- glmnet(
        X[start_i:end_i, ], Y[start_i:end_i],
        intercept = FALSE, family = "binomial",
        lambda    = λ_full
      )
      beta_full <- as.vector(coef(fit_full))[-1]
      I_full    <- rep(1, len_full)
      losses[start_i] <- sum(log(I_full + exp(X[start_i:end_i, ] %*% beta_full)) -
                             Y[start_i:end_i] * (X[start_i:end_i, ] %*% beta_full)) / n +
                        delta * c1 * ( sqrt(2 * log(2 * p) / len_full) + log(2 * p) / len_full )
      
      # try every possible split
      for (s in (start_i + m1):(end_i - m1)) {
        # left segment
        len_l   <- s - start_i + 1
        λ_l     <- c1 * ( sqrt(2 * log(2 * p) / len_l) + log(2 * p) / len_l )
        fit_l   <- glmnet(
          X[start_i:s, ], Y[start_i:s],
          intercept = FALSE, family = "binomial",
          lambda    = λ_l
        )
        beta_l  <- as.vector(coef(fit_l))[-1]
        I_l     <- rep(1, len_l)
        loss_l  <- sum(log(I_l + exp(X[start_i:s, ] %*% beta_l)) -
                       Y[start_i:s] * (X[start_i:s, ] %*% beta_l)) / n +
                   delta * c1 * ( sqrt(2 * log(2 * p) / len_l) + log(2 * p) / len_l )
        
        # right segment
        len_r   <- end_i - s
        λ_r     <- c1 * ( sqrt(2 * log(2 * p) / len_r) + log(2 * p) / len_r )
        fit_r   <- glmnet(
          X[(s+1):end_i, ], Y[(s+1):end_i],
          intercept = FALSE, family = "binomial",
          lambda    = λ_r
        )
        beta_r  <- as.vector(coef(fit_r))[-1]
        I_r     <- rep(1, len_r)
        loss_r  <- sum(log(I_r + exp(X[(s+1):end_i, ] %*% beta_r)) -
                       Y[(s+1):end_i] * (X[(s+1):end_i, ] %*% beta_r)) / n +
                   delta * c1 * ( sqrt(2 * log(2 * p) / len_r) + log(2 * p) / len_r )
        
        losses[s] <- loss_l + loss_r
      }
      
      shat[i] <- which.min(losses)
      if (shat[i] != start_i) {
        alpha <- sort(unique(c(alpha, shat[i])))
      }
    }
    # stop early if no new splits added
    if (identical(alpha, previous_alpha)) break
  }
  
  cpnumber.estimator   <- length(alpha) - 2
  cplocation.estimator <- alpha
  return(list(
    cpnumber.estimator   = cpnumber.estimator,
    cplocation.estimator = cplocation.estimator
  ))
}


#' This function provides change point estimation by Binary Segmentation Algorithm
#' for the high dimensional generalized linear regression model.
#' @export
#' @param X an n*p data matrix, n=observation, p=dimension
#' @param Y an n dimensional data vector, n=observation
#' @param kmax the maximum change point number.
#' @param c1 a constant specified by users.
#' @param m1 the sample size of the minimal interval.
#' @param delta the minimal interval length specified by users.
#' @return A list including: the estimated change point number, the estimated
#'   change point locations.
#' @importFrom stats coef
#' @importFrom glmnet glmnet
######## change point estimation by bsa ########
bsa_exactchgpt <- function(X, Y, kmax, c1, m1, delta) {
  print("Entered bsa_exactchgpt")
  n <- length(Y)
  p <- ncol(X)
  alpha <- c(1, n)
  shat  <- numeric()
  
  # enforce exactly kmax splits
  for (z in kmax:kmax) {
    previous_alpha <- alpha
    for (i in seq_len(length(previous_alpha) - 1)) {
      start_i <- previous_alpha[i]
      end_i   <- previous_alpha[i + 1]
      if ((start_i + m1) > (end_i - m1)) next
      
      losses <- rep(Inf, end_i)
      
      # full-segment fit
      len_full <- end_i - start_i + 1
      λ_full   <- c1 * ( sqrt(2 * log(2 * p) / len_full) + log(2 * p) / len_full )
      fit_full <- glmnet(
        X[start_i:end_i, ], Y[start_i:end_i],
        intercept = FALSE, family = "binomial",
        lambda    = λ_full
      )
      beta_full <- as.vector(coef(fit_full))[-1]
      I_full    <- rep(1, len_full)
      losses[start_i] <- sum(log(I_full + exp(X[start_i:end_i, ] %*% beta_full)) -
                             Y[start_i:end_i] * (X[start_i:end_i, ] %*% beta_full)) / n +
                        delta * c1 * ( sqrt(2 * log(2 * p) / len_full) + log(2 * p) / len_full )
      
      # try every possible split
      for (s in (start_i + m1):(end_i - m1)) {
        # left segment
        len_l   <- s - start_i + 1
        λ_l     <- c1 * ( sqrt(2 * log(2 * p) / len_l) + log(2 * p) / len_l )
        fit_l   <- glmnet(
          X[start_i:s, ], Y[start_i:s],
          intercept = FALSE, family = "binomial",
          lambda    = λ_l
        )
        beta_l  <- as.vector(coef(fit_l))[-1]
        I_l     <- rep(1, len_l)
        loss_l  <- sum(log(I_l + exp(X[start_i:s, ] %*% beta_l)) -
                       Y[start_i:s] * (X[start_i:s, ] %*% beta_l)) / n +
                   delta * c1 * ( sqrt(2 * log(2 * p) / len_l) + log(2 * p) / len_l )
        
        # right segment
        len_r   <- end_i - s
        λ_r     <- c1 * ( sqrt(2 * log(2 * p) / len_r) + log(2 * p) / len_r )
        fit_r   <- glmnet(
          X[(s+1):end_i, ], Y[(s+1):end_i],
          intercept = FALSE, family = "binomial",
          lambda    = λ_r
        )
        beta_r  <- as.vector(coef(fit_r))[-1]
        I_r     <- rep(1, len_r)
        loss_r  <- sum(log(I_r + exp(X[(s+1):end_i, ] %*% beta_r)) -
                       Y[(s+1):end_i] * (X[(s+1):end_i, ] %*% beta_r)) / n +
                   delta * c1 * ( sqrt(2 * log(2 * p) / len_r) + log(2 * p) / len_r )
        
        losses[s] <- loss_l + loss_r
      }
      
      shat[i] <- which.min(losses)
      if (shat[i] != start_i) {
        alpha <- sort(unique(c(alpha, shat[i])))
      }
    }
    # break if no change
    if (identical(alpha, previous_alpha)) break
  }
  
  print("Reached end of bsa_exactchgpt")
  cpnumber.estimator   <- length(alpha) - 2
  cplocation.estimator <- alpha
  return(list(
    cpnumber.estimator   = cpnumber.estimator,
    cplocation.estimator = cplocation.estimator
  ))
}

#' This function wraps around bsachgpt.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rnorm runif
#' @param n p=dimension
#' @param p n=observation
#' @param Sigma pxp covariance of rows of X
#' @param kmax the maximum change point number.
#' @param delta the minimal interval length specified by users.
#' @return  A list including: the estimated change point number, the estimated
#'   change point locations.
bsawrapper<-function(n, p, Sigma, kmax, delta){
  #=================== initial settings====================
  # n=200       # sample size 
  # p=10    # data dimension
  # Sigma<-matrix(0,nrow = p,ncol = p)   #covariance matrix 
  # for (i in 1:p){
  #   for (j in 1:p){
  #     Sigma[i,j]<-0.6^(abs(i-j))
  #   }
  # }
  print("Entered bsawrapper")
  tau0<-0.5  # true change point location
  # kmax=6   # upper bound of the number of change point
  #==== tuning parameter(specifed by users)==========
  # delta=0.1  # fractional length of the shortest interval 
  # (i.e. distance between nearest adjacent changepoints as a fraction of n)
  c1<-0.18 # lambda related parameter
  #===============signal jump size==================
  m1=ceiling(n*delta)   # size of the shortest interval
  signaljump<-20*sqrt(log(p)/(delta*n))  # signal jump
  signalsupport<-ceiling(log(p))
  signalsupport.range<-0.3*p

  #===============generating change point model settings============
  beta1<-rep(0,p)   #  regression coefficient 1
  beta2<-rep(0,p)   #  regression coefficient 2
  for (i in sample(1:(signalsupport.range),signalsupport)) {
    beta1[i]<-runif(1,min=0,max=2)
  }
  for (i in sample(1:(signalsupport.range),signalsupport)) {
    beta2[i]<-runif(1,min=0,max=2)+runif(1,min=0.5,max=signaljump)
  }

  X<-mvrnorm(n,rep(0,p),Sigma)  # design matrix
  error<-rnorm(n)      # error term
  Y1<-vector(mode = "numeric",length=0L)  # response variable before change point
  Y2<-vector(mode = "numeric",length=0L)  # response variable after change point
  O1<-vector(mode = "numeric",length=0L)
  O2<-vector(mode = "numeric",length=0L)
  R1=X[1:(tau0*n),]%*%beta1+error[1:(tau0*n)]
  R2=X[(tau0*n+1):n,]%*%beta2+error[(tau0*n+1):n]
  for (i in 1:(n/2)){
    O1[i]<-exp(R1[i])/(1+exp(R1[i]))
    O2[i]<-exp(R2[i])/(1+exp(R2[i]))
    Y1[i]=rbinom(1,1,O1[i])
    Y2[i]=rbinom(1,1,O2[i])
  }

  Y=c(Y1,Y2)   # response variable 


  #=======estimate change points by bsa=======
  reg<-bsa_chgpt(X,Y,kmax,c1,m1,delta)
  cpnumber<-reg$cpnumber.estimator
  cplocation<-reg$cplocation.estimator
  #=========output==========
  print("Reached end of bsawrapper")
  print(cpnumber)
  print(cplocation)
  list<-list(cpnumber=cpnumber,cplocation=cplocation)
  return(list)
}



#' This function wraps around bsa.chgpt.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rnorm runif
#' @param X an n*p data matrix, n=observation, p=dimension
#' @param Y an n dimensional data vector, n=observation
#' @param kmax the maximum change point number.
#' @param delta the minimal interval length specified by users.
#' @return  A list including: the estimated change point number, the estimated
#'   change point locations.
bsa_wrapper<-function(X, Y, kmax, delta){
  #=================== initial settings====================
  # n=200       # sample size 
  # p=10    # data dimension
  # Sigma<-matrix(0,nrow = p,ncol = p)   #covariance matrix 
  # for (i in 1:p){
  #   for (j in 1:p){
  #     Sigma[i,j]<-0.6^(abs(i-j))
  #   }
  # }
  print("Entered bsa_wrapper")
  n<-length(Y)
  # tau0<-0.5  # true change point location
  # kmax=6   # upper bound of the number of change point
  #==== tuning parameter(specifed by users)==========
  # delta=0.1  # fractional length of the shortest interval 
  # (i.e. distance between nearest adjacent changepoints as a fraction of n)
  c1<-0.18 # lambda related parameter
  #===============signal jump size==================
  m1=ceiling(n*delta)   # size of the shortest interval

  #=======estimate change points by bsa=======
  reg<-bsa_chgpt(X,Y,kmax,c1,m1,delta)
  cpnumber<-reg$cpnumber.estimator
  cplocation<-reg$cplocation.estimator
  #=========output==========
  print("Reached end of bsa_wrapper")
  print(cpnumber)
  print(cplocation)
  list<-list(cpnumber=cpnumber,cplocation=cplocation)
  return(list)
}



#' This function wraps around bsa.exactchgpt.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rnorm runif
#' @param X an n*p data matrix, n=observation, p=dimension
#' @param Y an n dimensional data vector, n=observation
#' @param kmax the maximum change point number.
#' @param delta the minimal interval length specified by users.
#' @return  A list including: the estimated change point number, the estimated
#'   change point locations.
bsa_exact_wrapper<-function(X, Y, kmax, delta){
  #=================== initial settings====================
  # n=200       # sample size 
  # p=10    # data dimension
  # Sigma<-matrix(0,nrow = p,ncol = p)   #covariance matrix 
  # for (i in 1:p){
  #   for (j in 1:p){
  #     Sigma[i,j]<-0.6^(abs(i-j))
  #   }
  # }
  print("Entered bsa_exact_wrapper")
  n<-length(Y)
  # tau0<-0.5  # true change point location
  # kmax=6   # upper bound of the number of change point
  #==== tuning parameter(specifed by users)==========
  # delta=0.1  # fractional length of the shortest interval 
  # (i.e. distance between nearest adjacent changepoints as a fraction of n)
  c1<-0.18 # lambda related parameter
  #===============signal jump size==================
  m1=ceiling(n*delta)   # size of the shortest interval
  
  #=======estimate change points by bsa=======
  reg<-bsa_exactchgpt(X,Y,kmax,c1,m1,delta)
  cpnumber<-reg$cpnumber.estimator
  cplocation<-reg$cplocation.estimator
  #=========output==========
  print("Reached end of bsa_exact_wrapper")
  print(cpnumber)
  print(cplocation)
  list<-list(cpnumber=cpnumber,cplocation=cplocation)
  return(list)
}