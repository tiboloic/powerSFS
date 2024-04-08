# RTMB multinomial test

n = 100
p = 20

# simulate data
bs = 1 + rbeta(n, 2, 2)
ps <- sapply(bs, function(b) {1/(1:p)^b})
obs <- matrix(
  sapply(1:n, function(i) rmultinom(1, rpois(1, 1000), prob=ps[,i])),
  p,
  n)

par <- list(
  mu = 2,
  sig = 1,
  bs = rep(2, n)
)

nLL <- function(params) {
  
  dm <- function(x,p) {
    lgamma(sum(x) + 1) - sum(lgamma(x+1)) + sum(x*log(p))
  }
  getAll(params)
  
  obs <- OBS(obs)
  
  nLL <- 0.0
  nLL <- nLL - sum(dnorm(bs, mu, sig, TRUE))
  
  for (i in 1:ncol(obs)) {
    beta <- bs[i]
    preds <- 1/(1:p)^beta
    preds <- preds/sum(preds)
    
    nLL <- nLL - dmultinom(obs[,i], prob=preds, log = TRUE)
    #nLL <- nLL - dm(obs[,i], preds)
  }
  nLL
}

obj <- MakeADFun(nLL, par, random = c("bs"))

fit <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
ranef <- summary(rep, "random")
plot(ranef[,1] ~ bs)
abline(a=0, b=1, col=2, lwd=2)
