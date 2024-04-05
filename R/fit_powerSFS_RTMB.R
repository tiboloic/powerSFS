# LT 3/04/2024
# 
# first implementation of powerSFS with RTMB
# 
#

library(RTMB)
library(tidyverse)

sfs <- read_delim("gnomADv4_protein_altering_sfs_prototype.tsv.gz")
n_genes <- sfs %>% count(gene) %>% tally() %>% as.numeric()

obs <- sfs %>% select(gene, octave, count)  %>% spread(key=octave, value=count, fill=0)

parameters <- list(
  mub = -1,
  sdb = 1,
  bs = rep(0, n_genes)
)

# negative log-likelihood
nLL <- function(params, eps = 1e-8) {
  nLL <- 0.0
  
  # random deviations in SFS slopes
  nLL <- nLL - sum(dnorm(bs, mean = mub, sd - sdb, log = TRUE))
  
  # Bernouilli polunomials evaluated at zero
  B2 <- c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510)
  n_B <- length(B2)
  
  # proportion of missense predicted by bin
  # initialise with small probability to avoid limit case where p = 0
  preds <- rep(eps, ncol(obs) - 1)
  
  # calculate first 15 terms of the series
  for (i in 1:15) {
    preds[floor(log2(i))] = preds[floor(log2(i))] + i^beta
  }
  
  # use asymptotic approximation for the remaining 1 000 000+ terms
  #  build lower and upper limits of bins
  
  last_bin <- ceiling(log2(an[i]))
  lower <- 2 ^ (1:last_bin - 1)
  upper <- 2 ^ (1:last_bin) -1
  zeta <- 
  
  }