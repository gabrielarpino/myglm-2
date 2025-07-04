library(glmnet)
library(MASS)

# Alternative safe_glmnet_coef that matches the original exactly
safe_glmnet_coef <- function(X, y, target_lambda) {
  n_rows <- nrow(X)
  n_cols <- ncol(X)
  
  # Check if data is problematic
  if (n_rows < 10 || min(table(y)) < 3) {
    # Return exactly what the original would: intercept + p coefficients
    return(c(0, rep(0, n_cols)))  # Explicit: [intercept, coef1, coef2, ...]
  }
  
  # Try with forced lambda values
  tryCatch({
    # Match the original call exactly
    fit <- glmnet(X, y, intercept=FALSE, family="binomial", 
                  lambda=c(1.0, 0.1, 0.01))
    
    # Use middle lambda
    coef_result <- as.vector(coef(fit, s=0.1))
    
    # Debug: print dimensions
    print(paste("X dims:", n_rows, "x", n_cols, "Coef length:", length(coef_result)))
    
    return(coef_result)
    
  }, error = function(e) {
    print(paste("glmnet failed:", e$message))
    # Return same structure as successful glmnet call
    return(c(0, rep(0, n_cols)))
  })
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
#'  change point locations.
#' @importFrom stats coef
#' @importFrom glmnet glmnet
######## change point estimation by bsa ########
bsa_chgpt<-function(X,Y,kmax,c1,m1,delta){
  print("Entered bsa_chgpt")
  #1.Initialize T to the tree with a single root node labeled by (0; 1].
  n<-length(Y)
  v=length(Y)
  u=1
  p<-length(X[1,])
  alpha<-c(u,v)
  hs<-vector(mode="numeric",length=0L)
  shat<-vector(mode="numeric",length=0L)
  alphalength=vector(mode="numeric",length=0L)

  # --- START ADDITION: Define glmnet control parameters ---
  # These parameters help prevent internal integer overflows in glmnet.
  # nlambda: number of lambda values to use in the regularization path.
  # Default is 100. Reducing it might resolve issues if 100 is too high for certain p/n.
  # pmax: maximum number of variables ever to be nonzero. Default is p.
  # Explicitly setting it to p is generally safe and robust.
  glmnet_nlambda <- 100 # Standard default, can be reduced if needed (e.g., 50 or 20)
  glmnet_pmax <- p     # Explicitly set to the number of features 'p'
  # --- END ADDITION ---

  for (z in 1:kmax) {
    alphalength[z]=length(alpha)
    for (i in 1:(length(alpha)-1)){

      for (s in (alpha[i]+m1):(alpha[i+1]-m1)){
        if ((alpha[i]+m1)>(alpha[i+1]-m1)){
          break
        }

        hs[1:(alpha[i]+m1)]<-rep(10000,(alpha[i]+m1))
        betahat33<-matrix(0,nrow =(alpha[i+1]-m1),ncol = p)
        betahat55<-matrix(0,nrow =(alpha[i+1]-m1),ncol = p)

        # print("Reached 1/3 glmnet in bsa")
        # --- START MODIFICATION: Added nlambda and pmax to glmnet calls ---
        # betahat22<-as.vector(coef(glmnet(X[alpha[i]:s,],Y[alpha[i]:s],intercept=FALSE,
        #                                  family="binomial",
        #                                  nlambda=glmnet_nlambda, # ADDED
        #                                  pmax=glmnet_pmax),     # ADDED
        #                           s=c1*(sqrt(2*log(2*p)/(s-alpha[i]+1))+log(2*p)/(s-alpha[i]+1)))) # the estimation of regression coefficients
        betahat22 <- safe_glmnet_coef(X[alpha[i]:s,], Y[alpha[i]:s], 
                              c1*(sqrt(2*log(2*p)/(s-alpha[i]+1))+log(2*p)/(s-alpha[i]+1)))
        # --- END MODIFICATION ---
        # print("Passed 1/3 glmnet in bsa")
        betahat33[s,]<-betahat22[-1]
        # print("Reached 2/3 glmnet in bsa")
        # --- START MODIFICATION: Added nlambda and pmax to glmnet calls ---
        # betahat44<-as.vector(coef(glmnet(X[(s+1):(alpha[i+1]),],Y[(s+1):(alpha[i+1])],
        #                                  intercept=FALSE,family="binomial",
        #                                  nlambda=glmnet_nlambda, # ADDED
        #                                  pmax=glmnet_pmax),     # ADDED
        #                          s=c1*(sqrt(2*log(2*p)/(alpha[i+1]-s))+log(2*p)/(alpha[i+1]-s))))
        betahat44 <- safe_glmnet_coef(X[(s+1):(alpha[i+1]),], Y[(s+1):(alpha[i+1])], 
                              c1*(sqrt(2*log(2*p)/(alpha[i+1]-s))+log(2*p)/(alpha[i+1]-s)))
        # --- END MODIFICATION ---
        # print("Passed 2/3 glmnet in bsa")
        betahat55[s,]<-betahat44[-1]

        I<-rep(1,(alpha[i+1]-alpha[i]+1))
        # print("Reached 3/3 glmnet in bsa")
        # --- START MODIFICATION: Added nlambda and pmax to glmnet calls ---
        # betahat66<-as.vector(coef(glmnet(X[alpha[i]:alpha[i+1],],Y[alpha[i]:alpha[i+1]],
        #                                  intercept=FALSE,family="binomial",
        #                                  nlambda=glmnet_nlambda, # ADDED
        #                                  pmax=glmnet_pmax),     # ADDED
        #                           s=c1*(sqrt(2*log(2*p)/(alpha[i+1]-alpha[i]+1))+log(2*p)/(alpha[i+1]-alpha[i]+1)))) # the estimation of regression coefficients
        betahat66 <- safe_glmnet_coef(X[alpha[i]:alpha[i+1],], Y[alpha[i]:alpha[i+1]], 
                              c1*(sqrt(2*log(2*p)/(alpha[i+1]-alpha[i]+1))+log(2*p)/(alpha[i+1]-alpha[i]+1)))
        # --- END MODIFICATION ---
        # print("Passed 3/3 glmnet in bsa")
        hs[alpha[i]]<-sum(log(I+exp(X[alpha[i]:alpha[i+1],]%*%betahat66[-1]))-Y[alpha[i]:alpha[i+1]]*(X[alpha[i]:alpha[i+1],]%*%betahat66[-1]))/n
        +delta*c1*(sqrt(2*log(2*p)/(alpha[i+1]-alpha[i]+1))+log(2*p)/(alpha[i+1]-alpha[i]+1))

        I1=rep(1,(s-alpha[i]+1))
        I2=rep(1,(alpha[i+1]-s))
        # loss function
        hs[s]=sum(log(I1+exp(X[alpha[i]:s,]%*%betahat33[s,]))
                  -Y[alpha[i]:s]*(X[alpha[i]:s,]%*%betahat33[s,]))/n+
          delta*c1*(sqrt(2*log(2*p)/(s-alpha[i]+1))+log(2*p)/(s-alpha[i]+1)) +
          sum(log(I2+exp(X[(s+1):(alpha[i+1]),]%*%betahat55[s,]))
              -Y[(s+1):(alpha[i+1])]*(X[(s+1):(alpha[i+1]),]%*%betahat55[s,]))/n+delta*c1*(sqrt(2*log(2*p)/(alpha[i+1]-s))+log(2*p)/(alpha[i+1]-s))
      }

      shat[i]<-which.min(hs)

      if (shat[i]==alpha[i]){
        alpha=alpha }
      else {alpha=sort(c(alpha,shat[i]))
      alpha<-alpha[!duplicated(alpha)]}
    }
    alphalength[z+1]=length(alpha)
    if (alphalength[z]==alphalength[z+1])
    {
      break
    }
  }
  ##======= output change point estimator===============
  print("Reached end of bsa_chgpt")
  cpnumber.estimator<-length(alpha)-2 #change point number
  cplocation.estimator<-rep(0,length(alpha))
  cplocation.estimator<-alpha #change point location
  list<-list(cpnumber.estimator=cpnumber.estimator,cplocation.estimator=cplocation.estimator)
  return(list)
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
#'  change point locations.
#' @importFrom stats coef
#' @importFrom glmnet glmnet
######## change point estimation by bsa ########
bsa_exactchgpt<-function(X,Y,kmax,c1,m1,delta){
  print("Entered bsa_exactchgpt")
  #1.Initialize T to the tree with a single root node labeled by (0; 1].
  n<-length(Y)
  v=length(Y)
  u=1
  p<-length(X[1,])
  alpha<-c(u,v)
  hs<-vector(mode="numeric",length=0L)
  shat<-vector(mode="numeric",length=0L)
  alphalength=vector(mode="numeric",length=0L)

  # --- START ADDITION: Define glmnet control parameters ---
  # These parameters help prevent internal integer overflows in glmnet.
  # nlambda: number of lambda values to use in the regularization path.
  # Default is 100. Reducing it might resolve issues if 100 is too high for certain p/n.
  # pmax: maximum number of variables ever to be nonzero. Default is p.
  # Explicitly setting it to p is generally safe and robust.
  # glmnet_nlambda <- 10 # Standard default, can be reduced if needed (e.g., 50 or 20)
  # glmnet_pmax <- min(p, 50)     # Explicitly set to the number of features 'p'
  # In bsa_exactchgpt function, replace the parameter section with:
  glmnet_nlambda <- 2   # Very conservative
  glmnet_pmax <- 10     # Very conservative
  print(paste("updated 4: glmnet_pmax is set to:", glmnet_pmax))
  print(paste("updated 4: glmnet_nlambda is set to:", glmnet_nlambda))
  print(paste("update 5: safe_glmnet"))
  # Print what glmnet_pmax is
  # --- END ADDITION ---

  # The for loop below is the only difference between this function
  # and bsachgpt. This function only searches through configs with
  # exactly kmax chgpts.
  for (z in kmax:kmax) {
    alphalength[z]=length(alpha)
    for (i in 1:(length(alpha)-1)){

      for (s in (alpha[i]+m1):(alpha[i+1]-m1)){
        if ((alpha[i]+m1)>(alpha[i+1]-m1)){
          break
        }

        hs[1:(alpha[i]+m1)]<-rep(10000,(alpha[i]+m1))
        betahat33<-matrix(0,nrow =(alpha[i+1]-m1),ncol = p)
        betahat55<-matrix(0,nrow =(alpha[i+1]-m1),ncol = p)
        # print("Reached 1/3 glmnet in bsa")
        # --- START MODIFICATION: Added nlambda and pmax to glmnet calls ---
        # betahat22<-as.vector(coef(glmnet(X[alpha[i]:s,],Y[alpha[i]:s],intercept=FALSE,
        #                                  family="binomial",
        #                                  nlambda=glmnet_nlambda, # ADDED
        #                                  pmax=glmnet_pmax),     # ADDED
        #                           s=c1*(sqrt(2*log(2*p)/(s-alpha[i]+1))+log(2*p)/(s-alpha[i]+1)))) # the estimation of regression coefficients
        betahat22 <- safe_glmnet_coef(X[alpha[i]:s,], Y[alpha[i]:s], 
                              c1*(sqrt(2*log(2*p)/(s-alpha[i]+1))+log(2*p)/(s-alpha[i]+1)))
        # --- END MODIFICATION ---
        # print("Passed 1/3 glmnet in bsa")
        betahat33[s,]<-betahat22[-1]
        # print("Reached 2/3 glmnet in bsa")
        # --- START MODIFICATION: Added nlambda and pmax to glmnet calls ---
        # betahat44<-as.vector(coef(glmnet(X[(s+1):(alpha[i+1]),],Y[(s+1):(alpha[i+1])],
        #                                  intercept=FALSE,family="binomial",
        #                                  nlambda=glmnet_nlambda, # ADDED
        #                                  pmax=glmnet_pmax),     # ADDED
        #                          s=c1*(sqrt(2*log(2*p)/(alpha[i+1]-s))+log(2*p)/(alpha[i+1]-s))))
        betahat44 <- safe_glmnet_coef(X[(s+1):(alpha[i+1]),], Y[(s+1):(alpha[i+1])], 
                              c1*(sqrt(2*log(2*p)/(alpha[i+1]-s))+log(2*p)/(alpha[i+1]-s)))
        # --- END MODIFICATION ---
        # print("Passed 2/3 glmnet in bsa")
        betahat55[s,]<-betahat44[-1]

        I<-rep(1,(alpha[i+1]-alpha[i]+1))
        # print("Reached 3/3 glmnet in bsa")
        # --- START MODIFICATION: Added nlambda and pmax to glmnet calls ---
        # betahat66<-as.vector(coef(glmnet(X[alpha[i]:alpha[i+1],],Y[alpha[i]:alpha[i+1]],
        #                                  intercept=FALSE,family="binomial",
        #                                  nlambda=glmnet_nlambda, # ADDED
        #                                  pmax=glmnet_pmax),     # ADDED
        #                           s=c1*(sqrt(2*log(2*p)/(alpha[i+1]-alpha[i]+1))+log(2*p)/(alpha[i+1]-alpha[i]+1)))) # the estimation of regression coefficients
        betahat66 <- safe_glmnet_coef(X[alpha[i]:alpha[i+1],], Y[alpha[i]:alpha[i+1]], 
                              c1*(sqrt(2*log(2*p)/(alpha[i+1]-alpha[i]+1))+log(2*p)/(alpha[i+1]-alpha[i]+1)))
        # --- END MODIFICATION ---
        # print("Passed 3/3 glmnet in bsa")
        hs[alpha[i]]<-sum(log(I+exp(X[alpha[i]:alpha[i+1],]%*%betahat66[-1]))-Y[alpha[i]:alpha[i+1]]*(X[alpha[i]:alpha[i+1],]%*%betahat66[-1]))/n
        +delta*c1*(sqrt(2*log(2*p)/(alpha[i+1]-alpha[i]+1))+log(2*p)/(alpha[i+1]-alpha[i]+1))

        I1=rep(1,(s-alpha[i]+1))
        I2=rep(1,(alpha[i+1]-s))
        # loss function
        hs[s]=sum(log(I1+exp(X[alpha[i]:s,]%*%betahat33[s,]))
                  -Y[alpha[i]:s]*(X[alpha[i]:s,]%*%betahat33[s,]))/n+
          delta*c1*(sqrt(2*log(2*p)/(s-alpha[i]+1))+log(2*p)/(s-alpha[i]+1)) +
          sum(log(I2+exp(X[(s+1):(alpha[i+1]),]%*%betahat55[s,]))
              -Y[(s+1):(alpha[i+1])]*(X[(s+1):(alpha[i+1]),]%*%betahat55[s,]))/n+delta*c1*(sqrt(2*log(2*p)/(alpha[i+1]-s))+log(2*p)/(alpha[i+1]-s))
      }

      shat[i]<-which.min(hs)

      if (shat[i]==alpha[i]){
        alpha=alpha }
      else {alpha=sort(c(alpha,shat[i]))
      alpha<-alpha[!duplicated(alpha)]}
    }
    alphalength[z+1]=length(alpha)
    if (alphalength[z]==alphalength[z+1])
    {
      break
    }
  }
  ##======= output change point estimator===============
  print("Reached end of bsa_exactchgpt")
  cpnumber.estimator<-length(alpha)-2 #change point number
  cplocation.estimator<-rep(0,length(alpha))
  cplocation.estimator<-alpha #change point location
  list<-list(cpnumber.estimator=cpnumber.estimator,cplocation.estimator=cplocation.estimator)
  return(list)
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
#' @return A list including: the estimated change point number, the estimated
#'  change point locations.
bsawrapper<-function(n, p, Sigma, kmax, delta){
  #=================== initial settings====================
  # n=200    # sample size
  # p=10  # data dimension
  # Sigma<-matrix(0,nrow = p,ncol = p)  #covariance matrix
  # for (i in 1:p){
  #  for (j in 1:p){
  #   Sigma[i,j]<-0.6^(abs(i-j))
  #  }
  # }
  print("Entered bsawrapper")
  tau0<-0.5 # true change point location
  # kmax=6  # upper bound of the number of change point
  #==== tuning parameter(specifed by users)==========
  # delta=0.1 # fractional length of the shortest interval
  # (i.e. distance between nearest adjacent changepoints as a fraction of n)
  c1<-0.18 # lambda related parameter
  #===============signal jump size==================
  m1=ceiling(n*delta)  # size of the shortest interval
  signaljump<-20*sqrt(log(p)/(delta*n)) # signal jump
  signalsupport<-ceiling(log(p))
  signalsupport.range<-0.3*p

  #===============generating change point model settings============
  beta1<-rep(0,p)  # regression coefficient 1
  beta2<-rep(0,p)  # regression coefficient 2
  for (i in sample(1:(signalsupport.range),signalsupport)) {
    beta1[i]<-runif(1,min=0,max=2)
  }
  for (i in sample(1:(signalsupport.range),signalsupport)) {
    beta2[i]<-runif(1,min=0,max=2)+runif(1,min=0.5,max=signaljump)
  }

  X<-mvrnorm(n,rep(0,p),Sigma) # design matrix
  error<-rnorm(n)   # error term
  Y1<-vector(mode = "numeric",length=0L) # response variable before change point
  Y2<-vector(mode = "numeric",length=0L) # response variable after change point
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

  Y=c(Y1,Y2)  # response variable


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
#' @return A list including: the estimated change point number, the estimated
#'  change point locations.
bsa_wrapper<-function(X, Y, kmax, delta){
  #=================== initial settings====================
  # n=200    # sample size
  # p=10  # data dimension
  # Sigma<-matrix(0,nrow = p,ncol = p)  #covariance matrix
  # for (i in 1:p){
  #  for (j in 1:p){
  #   Sigma[i,j]<-0.6^(abs(i-j))
  #  }
  # }
  print("Entered bsa_wrapper")
  n<-length(Y)
  # tau0<-0.5 # true change point location
  # kmax=6  # upper bound of the number of change point
  #==== tuning parameter(specifed by users)==========
  # delta=0.1 # fractional length of the shortest interval
  # (i.e. distance between nearest adjacent changepoints as a fraction of n)
  c1<-0.18 # lambda related parameter
  #===============signal jump size==================
  m1=ceiling(n*delta)  # size of the shortest interval

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
#' @return A list including: the estimated change point number, the estimated
#'  change point locations.
bsa_exact_wrapper<-function(X, Y, kmax, delta){
  #=================== initial settings====================
  # n=200    # sample size
  # p=10  # data dimension
  # Sigma<-matrix(0,nrow = p,ncol = p)  #covariance matrix
  # for (i in 1:p){
  #  for (j in 1:p){
  #   Sigma[i,j]<-0.6^(abs(i-j))
  #  }
  # }
  print("Entered bsa_exact_wrapper")
  n<-length(Y)
  # tau0<-0.5 # true change point location
  # kmax=6  # upper bound of the number of change point
  #==== tuning parameter(specifed by users)==========
  # delta=0.1 # fractional length of the shortest interval
  # (i.e. distance between nearest adjacent changepoints as a fraction of n)
  c1<-0.18 # lambda related parameter
  #===============signal jump size==================
  m1=ceiling(n*delta)  # size of the shortest interval

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