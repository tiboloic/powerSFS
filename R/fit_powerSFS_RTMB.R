# LT 3/04/2024
# 
# first implementation of powerSFS with RTMB
# 
#

library(RTMB)
library(tidyverse)

sfs <- read_delim("gnomADv4_protein_altering_sfs_prototype.tsv.gz")

# remove genes with low number of observations
sfs <- sfs %>% group_by(gene) %>% filter(max(octave) > 3)

obs <- sfs %>% 
  select(gene, octave, count)  %>% 
  spread(key=octave, value=count, fill=0) %>%
  arrange(gene)

# matrix of counts per gene
counts <- as.matrix(obs %>% ungroup %>% select(!gene))

dat <- data.frame(
  genes = obs %>% pull(gene),
  obs = I(counts),
  an = sfs %>% group_by(gene) %>%
    summarise(an = max(AN_max)) %>%
    arrange(gene) %>%
    pull(an)
)

parameters <- list(
  mub = 2.0,
  sdb = 1.0,
  bs = rep(2, nrow(obs))
)

# negative log-likelihood
nLL <- function(params) {

  eps <- 1e-8
  
  getAll(params, dat)
  
  nLL <- 0.0
  
  # random deviations in SFS slopes
  nLL <- nLL - sum(dnorm(bs, mean = mub, sd = sdb, log = TRUE))
  
  # Bernoulli polynomials evaluated at zero
  Bs <- c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510)
  #n_B <- length(Bs)
  
  for (i in 1:nrow(obs)) {
    
    # slope of SFS for gene i
    beta <- bs[i]

    #  build upper limits of bins
    last_bin <- floor(log2(an[i]))
    lower <- 2 ^ (0:last_bin)

    # upper limit of last bin is the maximum number of alleles observed
    # for that particular gene
    lower[last_bin + 2] <- an[i]

    # proportion of missense predicted by bin
    # initialise with small probability to avoid limit case where p = 0
    preds <- rep(eps, last_bin + 1)

    # calculate first 15 terms of the series
    for (j in 1:15) {
      preds[floor(log2(j)) + 1] <- preds[floor(log2(j)) + 1] + j ^ -beta
    }

    ## use asymptotic approximation for the remaining 1 000 000+ terms

    # zeta[i] = \sum_n^\infty{1/i^\beta}
    zetas <- rep(0, last_bin - 2)

    for (j in 1:(last_bin - 2)) {

      n <- lower[j + 4]

      ## calculate first 9 terms of the Euler-Maclaurin formula
      # term 1
      zeta <- n ^ (1 - beta) / (beta - 1) + 0.5 / n ^ beta

      # term 2
      term <- beta / 2 / n ^ (beta + 1)
      zeta <-  zeta + term * Bs[1];

      # terms 3 to 9
      for (k in 2:8) {
        term <- term * (beta + 2 * k - 3) * (beta + 2 * k - 2) / ((2 * k -1) * 2 * k * n^2);
        zeta <- zeta + term * Bs[k];
      }

      zetas[j] <- zeta
    }
    preds[5:(last_bin + 1)] <- -diff(zetas)

    # proportion of variants in each bin
    preds <- preds / sum(preds)

    # multinomial likelihood of observations given expected
    counts <- as.vector(obs[i, 1:(last_bin + 1)])
    nLL <- nLL - dmultinom(x = counts, prob = preds, log = TRUE)
    }
  nLL
}

obj <- MakeADFun(nLL, parameters, random = c("bs"))

fit <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
ranef <- summary(rep, "random")
fixef <- summary(rep, "fixed")
res <- tibble(
  genes = dat$genes,
  beta = ranef[, "Estimate"],
  sd_beta = ranef[, "Std. Error"],
  mub = fixef["mub", "Estimate"],
  sd_mub = fixef["mub", "Std. Error"],
  sdb = fixef["sdb", "Estimate"])

res <- res %>% mutate(
  z = (beta - mub) / sd_beta,
  p = 2 * pnorm(abs(z), lower.tail = FALSE)
)

write_tsv(res, file = "gnomADv4_protein_altering_sfs_prototype_fit.tsv.gz")
