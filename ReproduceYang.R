library(MASS)
################### Logistic Model ########################
Card = 0
Cover = 0
N = 1000
for (rept in 1:N) {

# Data Generating
n <- 25  0 # number of sample
p <- 12 #number of covariates
mu_0 <- rep(0,p)
Sigma_0 <- diag(1,p)
x <- matrix(rep(NA,p*n),nrow = p)
x <- mvrnorm(n ,mu = mu_0, Sigma = Sigma_0)
theta_0 <- c(rep(2,p/2)/(1:(p/2)),rep(0,p/2))
logit <- - x %*% theta_0
pr <- exp(logit)/(exp(logit)+1)
Y <- rep(NA,n)
sapply(1:n, function(k) {Y[k] <<- rbinom(1,1,pr[k])})
index <- rep(0,p)
F_Model <- glm(Y~x, family = binomial)
F_loglik <- logLik(F_Model)
MSCS <- list(rep(1,p))

for (i in 1:(2^p-2)){
  num <- i
  k <- 1
  index <- rep(0,p)
  while (num > 0){
    index[k] <- num %% 2
    num <- num %/% 2
    k = k + 1
  }
  lreg <- glm(Y~x[,(index*(1:p))],family=binomial)
  if (2*(-logLik(lreg)+F_loglik)<= qchisq(p = 0.90,df = p-sum(index))) {
    MSCS <- c(MSCS,list(index))
  }
}
if (list(c(rep(1,p/2),rep(0,p/2))) %in% MSCS) Cover = Cover + 1 
Card = Card + length(MSCS)/N
}
