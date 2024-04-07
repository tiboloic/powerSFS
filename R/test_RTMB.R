# RTMB multinomial test

n = 10
p = 20

# simulate data
obs <- t(rmultinom(n, 100, prob=p:1))

par <- list(
  mu = 2,
  sig = 1,
  bs = rep(2, n)
)

nLL <- function(params) {
  
  dm <- function(x,p) {
    lgamma(sum(x) + 1) - sum(lgamma(x+1)) + sum(x*log(p))
  }
  getAll(par)
  
  nLL <- 0.0
  nLL <- nLL - sum(dnorm(bs, mu, sig, TRUE))
  
  for (i in 1:nrow(obs)) {
    beta <- bs[i]
    preds <- 1/(1:p)^beta
    preds <- preds/sum(preds)
    
    #nLL <- nLL - dmultinom(obs[i,], prob=preds, log = TRUE)
    nLL <- nLL - dm(obs[i,], preds)
  }
  nLL
}

obj <- MakeADFun(nLL, par, random = c("bs"))

fit <- nlminb(obj$par, obj$fn, obj$gr)
