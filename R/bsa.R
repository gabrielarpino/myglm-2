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
  x = matrix(rnorm(100 * 20), 100, 20)
  g = sample(c(0,1), 100, replace = TRUE)
  fit = glmnet(x, g, family = "binomial")
  print("Passed first glmnet")
  print(paste(fit, "is the fit solution"))
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
        print("Reached 1/3 glmnet in bsa")
        betahat22<-as.vector(coef(glmnet(X[alpha[i]:s,],Y[alpha[i]:s],intercept=FALSE,
                                         family="binomial"),s=c1*(sqrt(2*log(2*p)/(s-alpha[i]+1))+log(2*p)/(s-alpha[i]+1)))) # the estimation of regression coefficients
        print("Passed 1/3 glmnet in bsa")
        betahat33[s,]<-betahat22[-1]
        print("Reached 2/3 glmnet in bsa")
        betahat44<-as.vector(coef(glmnet(X[(s+1):(alpha[i+1]),],Y[(s+1):(alpha[i+1])],
                                         intercept=FALSE,family="binomial"),
                                  s=c1*(sqrt(2*log(2*p)/(alpha[i+1]-s))+log(2*p)/(alpha[i+1]-s))))
        print("Passed 2/3 glmnet in bsa")
        betahat55[s,]<-betahat44[-1]
        
        I<-rep(1,(alpha[i+1]-alpha[i]+1))
        print("Reached 3/3 glmnet in bsa")
        betahat66<-as.vector(coef(glmnet(X[alpha[i]:alpha[i+1],],Y[alpha[i]:alpha[i+1]],
                                         intercept=FALSE,family="binomial"),s=c1*(sqrt(2*log(2*p)/(alpha[i+1]-alpha[i]+1))+log(2*p)/(alpha[i+1]-alpha[i]+1))))   # the estimation of regression coefficients
        print("Passed 3/3 glmnet in bsa")
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
  cpnumber.estimator<-length(alpha)-2  #change point number
  cplocation.estimator<-rep(0,length(alpha))
  cplocation.estimator<-alpha #change point location
  list<-list(cpnumber.estimator=cpnumber.estimator,cplocation.estimator=cplocation.estimator)
  return(list)
}
