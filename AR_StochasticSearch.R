library(MASS)
set.seed(2920)

PValue <- function(u) {
  model <- lm(Y~X[,(u*(1:p))])
  PValue <- pchisq(as.numeric(2*(-logLik(model)+F_logLik)),df = p-sum(u),lower.tail=FALSE,ncp=0)
}

ModelSampling <- function(w) {
  ModelSampling <- rep(0,p)
  for (j in 1:p){
    ran <- runif(1)
    if (ran<w[j]) {ModelSampling[j] <- 1}
  }
  ModelSampling
}

zeta <- 0.25 #Use p-value[(1-zeta)B] to update next alphe
alpha_star <- 0.05 #Target Level alpha
xi <- 0.2 #Smoothing Parameter
d <- 10 #The algorithm will stop when alpha is stablized for d iterations
p <- 100 #Largest AR order
B <- 300

ts <- arima.sim(n = 500, list(order = c(3,0,0),ar = c(0.8897, -0.4858, 0.3222)))
n <- length(ts)
Y <- ts[(p+1):n]
X <- matrix(rep(0,p*(n-p)),nrow = n-p)
sapply(1:(n-p),function(k){X[k,] <<- ts[k:(k+p-1)]})

F_Model <- lm(Y~X) #Full Model
F_logLik <- logLik(F_Model) #Log Likelihood of Full Model

t <- 1
w <- matrix(rep(0,50*p),nrow=50)
w[1,] <- rep(1/2,p)
stable <- 0
alpha <- rep(0,50)
c <- matrix(rep(0,p*50),nrow =50)
while (stable < 10) {
  t <- t + 1
  M <- NULL
  pval <- rep(0,B)
  for (i in 1:B){
    M <- rbind(M,ModelSampling(w[t-1,])) #Update Models using w^(t-1)
  }
  pval <- sapply(1:B, function(k){PValue(M[k,])})
  alpha[t] <- min(alpha_star, quantile(pval,1-zeta)) #Update alpha
  numer <- matrix(rep(0,p*B),nrow = B)
  denom <- matrix(rep(0,p*B),nrow = B)
  for (i in 1:B){
    for (j in 1:p){
      numer[i,j] <- as.numeric((pval[i]>alpha[t]) & (M[i,j] == 1))
      denom[i,j] <- as.numeric(pval[i]>alpha[t])
    }
  }
  
  c[t,] <- apply(numer,2,sum)/apply(denom,2,sum)
  w[t,] <- xi*c[t,] + (1-xi) * w[t-1,]
  if (alpha[t] == alpha_star) {stable <- stable + 1} else {stable = 0}
}

####Graph Plot#####
plot(1:t,w[1:t,1],ylim=c(0,1),type = "l",col="red",main=paste("AR(3) Model (p=",p,",n=",n,",xi=",xi,")"))
sapply(2:(p-3),function(k){lines(1:t,w[1:t,k],col="red")})
sapply((p-2):p,function(k){lines(1:t,w[1:t,k],col="black")})
