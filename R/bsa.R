library(glmnet)
library(MASS)

#' This function provides change point estimation by Binary Segmentation Algorithm
#'  for the high dimensional generalized linear regression model.
#' @export
#' @param X an n*p data matrix, n=observation, p=dimension
#' @param Y an n dimensional data vector, n=observation
#' @param kmax the maximum change point number.
#' @param c1 a constant specified by users.
#' @param m1 the sample size of the minimal interval.
#' @param delta the minimal interval length specified by users.
#' @return  A list including: the estimated change point number, the estimated
#'   change point locations.
#' @importFrom stats coef
#' @importFrom glmnet glmnet
######## change point estimation by bsa ########
bsachgpt<-function(X,Y,kmax,c1,m1,delta){
  print("Entered bsachgpt")
  #1.Initialize T to the tree with a single root node labeled by (0; 1].
  n<-length(Y)
  v=length(Y)
  u=1
  p<-length(X[1,])
  alpha<-c(u,v)
  hs<-vector(mode="numeric",length=0L)
  shat<-vector(mode="numeric",length=0L)
  alphalength=vector(mode="numeric",length=0L)
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
        betahat22<-as.vector(coef(glmnet(X[alpha[i]:s,],Y[alpha[i]:s],intercept=FALSE,
                                         family="binomial"),s=c1*(sqrt(2*log(2*p)/(s-alpha[i]+1))+log(2*p)/(s-alpha[i]+1)))) # the estimation of regression coefficients
        # print("Passed 1/3 glmnet in bsa")
        betahat33[s,]<-betahat22[-1]
        # print("Reached 2/3 glmnet in bsa")
        betahat44<-as.vector(coef(glmnet(X[(s+1):(alpha[i+1]),],Y[(s+1):(alpha[i+1])],
                                         intercept=FALSE,family="binomial"),
                                  s=c1*(sqrt(2*log(2*p)/(alpha[i+1]-s))+log(2*p)/(alpha[i+1]-s))))
        # print("Passed 2/3 glmnet in bsa")
        betahat55[s,]<-betahat44[-1]
        
        I<-rep(1,(alpha[i+1]-alpha[i]+1))
        # print("Reached 3/3 glmnet in bsa")
        betahat66<-as.vector(coef(glmnet(X[alpha[i]:alpha[i+1],],Y[alpha[i]:alpha[i+1]],
                                         intercept=FALSE,family="binomial"),s=c1*(sqrt(2*log(2*p)/(alpha[i+1]-alpha[i]+1))+log(2*p)/(alpha[i+1]-alpha[i]+1))))   # the estimation of regression coefficients
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
        alpha=alpha  }
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
  print("Reached end of bsachgpt")
  cpnumber.estimator<-length(alpha)-2  #change point number
  cplocation.estimator<-rep(0,length(alpha))
  cplocation.estimator<-alpha #change point location
  list<-list(cpnumber.estimator=cpnumber.estimator,cplocation.estimator=cplocation.estimator)
  return(list)
}

#' This function wraps around bsachgpt.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rnorm runif
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
  reg<-bsachgpt(X,Y,kmax,c1,m1,delta)
  cpnumber<-reg$cpnumber.estimator
  cplocation<-reg$cplocation.estimator
  #=========output==========
  print("Reached end of bsawrapper")
  print(cpnumber)
  print(cplocation)
  list<-list(cpnumber=cpnumber,cplocation=cplocation)
  return(list)
}



#' This function wraps around bsachgpt.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rnorm runif
bsawrapper1<-function(X, Y, kmax, delta){
  #=================== initial settings====================
  # n=200       # sample size 
  # p=10    # data dimension
  # Sigma<-matrix(0,nrow = p,ncol = p)   #covariance matrix 
  # for (i in 1:p){
  #   for (j in 1:p){
  #     Sigma[i,j]<-0.6^(abs(i-j))
  #   }
  # }
  print("Entered bsawrapper1")
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
  reg<-bsachgpt(X,Y,kmax,c1,m1,delta)
  cpnumber<-reg$cpnumber.estimator
  cplocation<-reg$cplocation.estimator
  #=========output==========
  print("Reached end of bsawrapper1")
  print(cpnumber)
  print(cplocation)
  list<-list(cpnumber=cpnumber,cplocation=cplocation)
  return(list)
}