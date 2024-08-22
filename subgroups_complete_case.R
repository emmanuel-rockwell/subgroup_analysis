rm(list=ls())
library(lme4)
library(Matrix)
library(rootSolve)
set.seed(1)
G0=3
G=3
N=100
T=50

p=0.4
sigma2=1
mu_x1=1
mu_c1=mu_x1+sqrt(2*sigma2)*qnorm(1-p)

# score_function <- function(X, Y, Z, sigma2, theta) {
#   beta0 <- theta[1]
#   beta1 <- theta[2]
#   beta2 <- theta[3]
#   
#   # sigma2 <- theta[3]
#   n <- length(Y)
#   resid <- Y - beta0 - beta1*X - beta2*Z
#   
#   score_beta0 <- sum(resid)/sigma2
#   score_beta1 <- sum(resid*X)/sigma2
#   score_beta2 <- sum(resid*Z)/sigma2
#   # score_sigma <- -n/(2*sigma2) + sum(resid^2)/(2*sigma2^2)
#   return(c(score_beta0, score_beta1, score_beta2))
# }

score_function <- function(X, Y, Z, sigma2, theta) {
  a <- rbind(1, X, Z)
  return(a %*% (rowSums(Y - t(a) %*% theta)/sigma2))
}

score_function <- function(X, Y, Z, sigma2, theta) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  
  n <- length(Y)
  resid <- Y - beta0 - beta1*X - beta2*Z
  
  score_beta0 <- sum(resid) / sigma2
  score_beta1 <- sum(resid * X) / sigma2
  score_beta2 <- sum(resid * Z) / sigma2
  
  return(c(score_beta0, score_beta1, score_beta2))
}


num=10 #number of repetition
max_it=10 #max number of iterations
tol=10^{-6} #iteration tolerance
rate=rep(0,num)
for (k in 1: num){
  b=list(c(1,0.4,1.6),c(1,1,1),c(1,1.6,0.4))
  datalist=list()
  true_index=c()
  for (i in 1:N){
    x1=rnorm(T, mean=mu_x1, sd=sqrt(sigma2))
    c1=rnorm(T, mean=mu_c1, sd=sqrt(sigma2))
    w1=pmin(x1, c1)
    d1=ifelse(x1<c1, 1, 0)
    x2=rnorm(T)
    if (i < N*0.3){
      true_index[i]=1
    }else if(i>=N*0.3 & i<N*0.6){
      true_index[i]=2
    }else {
      true_index[i]=3
    }
    beta=b[[true_index[i]]]
    y=beta[1]+beta[2]*x1+beta[3]*x2+rnorm(T)
    df=data.frame(x1=x1,c1=c1,w1=w1,d1=d1,x2=x2,y=y,subject=rep(i,T))
    datalist[[i]]=df
  }
  df=do.call(rbind,datalist)
  
  single_lm=function(df){
    # lmod=lm(y~x1+x2,df)
    # lmod$coefficients
    dfcc <- df[df$d1 == 1, ]
    mle.cc <- multiroot(
      f = function(b){score_function(X = dfcc$w1, Y = dfcc$y, Z = dfcc$x2, sigma2 = sigma2, theta = b)},
      start = c(0,0,0))$root
    mle.cc
  }
  
  lkh=function(df,beta){#likelihood function
    dfcc = df[df$d1 == 1, ]
    temp=dfcc$y
    temp=(dfcc$y-beta[1]-beta[2]*dfcc$x1-beta[3]*dfcc$x2)^2
    sum(-1*temp)
  }
  
  lkh <- function(df, beta) {
    dfcc = df[df$d1 == 1, ]
    sigma2 <- 1
    beta0 <- beta[1]
    beta1 <- beta[2]
    beta2 <- beta[3]
    n <- length(dfcc$y)
    resid <- dfcc$y - beta0 - beta1*dfcc$x1 - beta2*dfcc$x2
    log_likelihood <- sum(resid^2)/(2*sigma2)
    return(-log_likelihood)  # Return log-likelihood for optimization
  }
  
  #  ahat=rnorm(N)
  #  bhat=matrix(rnorm(2*G), nrow=2)
  
  bhat=sapply(datalist,single_lm)
  bhat=kmeans(t(bhat),centers=G)$centers
  bhat=t(bhat)
  it=1
  delta=1000
  while (it<max_it && delta>tol){
    temp=c()
    obhat=bhat
    for (i in 1:G){ #calculate likelihood
      beta=bhat[,i]
      temp=cbind(temp,sapply(datalist,lkh,beta=beta))
    }
    index=apply(temp,1,which.max) #determine which class belonging to
    for (i in 1:G){
      sublist=datalist[which(index==i)] #extract ith class
      df=do.call(rbind,sublist) #combine as a dataframe
      if (length(df)>0){
        # lmod=lm(y~x1+x2,data=df)
        # bhat[,i]=as.matrix(coef(lmod)[2:3])
        bhat[,i]=as.matrix(single_lm(df))
      }
    }
    it=it+1
    delta=mean((obhat-bhat)^2)
  }
  
  b=do.call(cbind,b)
  permu=rep(0,3)
  for (i in 1:G){
    permu[i]=which.min(colSums((b-bhat[,i])^2))
  }
  permu_index=index
  for (i in 1:G){
    permu_index[which(index==i)]=permu[i]
  }
  permu_bhat=bhat[,permu]
  cfrate=mean(permu_index==true_index)
  rate[k]=cfrate
  
  for (i in 1:G){
    permu[i]=which.min(colSums((b[,i]-bhat[])^2))
  }
  permu_bhat=bhat[,permu]
}
rate
