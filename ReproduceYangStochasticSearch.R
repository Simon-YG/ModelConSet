library(MASS)
set.seed(290)

PValue <- function(u) {
  lreg <- glm(Y~X[,(u*(1:p))],family = binomial)
  PValue <- pchisq(as.numeric(2*(-logLik(lreg)+F_logLik)),df = p-sum(u),lower.tail=FALSE,ncp=0)
}

ModelSampling <- function(w) {
  ModelSampling <- rep(0,p)
  for (j in 1:p){
    ran <- runif(1)
    if (ran<w[j]) {ModelSampling[j] <- 1}
  }
  ModelSampling
}
################### Logistic Model Stochastic Searching########################

B <- 300 #Number of Model to be Sampled
n <- 200 #Sample Size
p <- 100 #Number of Covariates
mu_0 <- rep(0,p) #Mean Vector for X
Sigma_0 <- matrix(sapply(0:(p^2-1),function(k) {0.5^(abs(k%%p-k%/%p))}),nrow = p) #Correlation Matrix for X
theta_0 <- c(c(0.5,-0.9,0.8,-2.2,3),rep(0,p-5)) #True paramter
zeta <- 0.25 #Use p-value[(1-zeta)B] to update next alphe
alpha_star <- 0.05 #Target Level alpha
xi <- 0.2 #Smoothing Parameter
d <- 30 #The algorithm will stop when alpha is stablized for d iterations

# Data Generating
X <- mvrnorm(n ,mu = mu_0, Sigma = Sigma_0)
logit <- (X) %*% theta_0
pr <- exp(logit)/(exp(logit)+1)
Y <- sapply(1:n, function(k) {rbinom(1,1,pr[k])})

#Full Model
F_Model <- glm(Y~X, family = binomial) #Full Model
F_logLik <- logLik(F_Model) #Log Likelihood of Full Model

t <- 1
w <- matrix(rep(0,50*p),nrow=50)
w[1,] <- rep(1/2,p)
stable <- 0
alpha <- rep(0,50)
c <- matrix(rep(0,p*50),nrow =50)
while (stable < d) {
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
plot(1:t,w[1:t,1],ylim=c(0,1),type = "l",main=paste("Logistic Regression (p=",p,",n=",n,"xi=",xi,")"))
sapply(2:5,function(k){lines(1:t,w[1:t,k])})
sapply(6:p,function(k){lines(1:t,w[1:t,k],col="red")})

