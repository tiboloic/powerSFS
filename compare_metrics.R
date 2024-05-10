# LT 8/05/2024
# 
# Compare powerSFS to other intolerance metric
# 

library(tidyverse)

# read alpha missense data
# 

alpha <- read_tsv("data/AlphaMissense_gene_hg38.tsv.gz", comment='#') 
constraint <- read_tsv("data/gnomad.v4.1.constraint_metrics.tsv")

# transcript Ensembl ID in alpha missense includes version of transcript (eg .xx)
mean(str_extract(alpha$transcript_id, "[^\\.]+") %in% constraint$transcript)
# 90.4 % of transcripts match gnomADv4.1 constraint

alpha <- alpha %>% mutate(transcript = str_extract(alpha$transcript_id, "[^\\.]+"))
dat <- alpha %>% inner_join(constraint)

psfs <- read_tsv("gnomADv4_protein_altering_sfs_prototype_fit.tsv.gz")
names(psfs)[1] <- "gene"
psfs <- psfs %>% mutate(beta = beta - mub)
dat <- dat %>% inner_join(psfs)

plot(mean_am_pathogenicity ~ mis.oe, data = dat)
plot(mean_am_pathogenicity ~ mis.oe_ci.upper, data = dat)
plot(beta ~ mis.oe, data = dat)
plot(beta ~ lof.oe, data = dat)
plot(-beta ~ lof.oe_ci.upper, data = dat, xlab = "LOEUF 2", ylab = "powerSFS")
m.1 <- lm(-beta ~ lof.oe_ci.upper, data = dat)
abline(m.1, lwd = 2, col = 2)
summary(m.1)
# R2 = 17%
summary(lm(beta ~ lof.oe, data = dat))
# R2 = 11%
summary(lm(beta ~ mis.oe, data = dat))
# 11.5 %
summary(lm(beta ~ mis.oe_ci.upper, data = dat))
# 11.2 %

qqplot(dat$beta, dat$lof.oe_ci.upper)

plot(beta ~ mean_am_pathogenicity, data = dat, xlab = "alpha missense", ylab = "powerSFS")
m.2 <- lm(beta ~ mean_am_pathogenicity, data = dat)
summary(m.2)
abline(m.2, lwd = 2, col = 2)
# R2 = 31%
summary(lm(z ~ mean_am_pathogenicity, data = dat))
# qqplot(dat$beta, dat$mean_am_pathogenicity)

m.am <- lm(beta ~ mean_am_pathogenicity, data = dat)
plot(residuals(m.am) ~ dat$mean_am_pathogenicity)
plot(residuals(m.am) ~ log(dat$cds_length))
summary(lm(residuals(m.am) ~ log(dat$cds_length)))

# is beta length biased compared to am ?
summary(lm(beta~ log(cds_length), data = dat))
summary(lm(mean_am_pathogenicity~ log(cds_length), data = dat))

# essential genes analysis
ess_genes <- read_tsv("data/AchillesCommonEssentialControls.csv") %>%
  transmute(gene = str_extract(Gene, "[^ ]+"), essential = 1)
control_genes <- read_tsv("data/AchillesNonessentialControls.csv") %>%
  transmute(gene = str_extract(Gene, "[^ ]+"), essential = 0)
ess_dat <- rbind(ess_genes, control_genes) %>% inner_join(dat)

# try an auRoc curve

par(cex = 1.2)
ess_dat_psfs <- ess_dat %>% arrange(desc(beta)) %>% transmute(beta.x = 1:nrow(ess_dat)/nrow(ess_dat), ess.y = cumsum(essential)/sum(ess_dat$essential))  
plot(ess.y ~ beta.x, data = ess_dat_psfs, type="l", col = 2, lwd = 2,
     xlab  = "Quantile of intolerance metric", ylab = "proportion of essential genes")
abline(a=0, b=1, col=1)
ess_dat_loeuf <- ess_dat %>% arrange(lof.oe_ci.upper) %>% transmute(loeuf.x = 1:nrow(ess_dat)/nrow(ess_dat), ess.y = cumsum(essential)/sum(ess_dat$essential))
lines(ess.y ~ loeuf.x, data = ess_dat_loeuf, col=3, lwd = 2)
ess_dat_am <- ess_dat %>% arrange(desc(mean_am_pathogenicity)) %>% transmute(am.x = 1:nrow(ess_dat)/nrow(ess_dat), ess.y = cumsum(essential)/sum(ess_dat$essential))
lines(ess.y ~ am.x, data = ess_dat_am, col=4, lwd = 2)
legend("topleft", c("LOEUF 2", "alpha missense", "powerSFS"), lwd = 2, col = c(3, 4, 2))

# restricting to underpowered genes
ess_dat_under <- filter(ess_dat, lof.exp < 10)
ess_dat_psfs <- ess_dat_under %>% arrange(desc(beta)) %>% transmute(beta.x = 1:nrow(ess_dat_under)/nrow(ess_dat_under), ess.y = cumsum(essential)/sum(ess_dat_under$essential))  
plot(ess.y ~ beta.x, data = ess_dat_psfs, type="l", col = 2, lwd = 2,
     xlab  = "Quantile of intolerance metric", ylab = "proportion of essential genes")
abline(a=0, b=1, col=1)
ess_dat_loeuf <- ess_dat_under %>% arrange(lof.oe_ci.upper) %>% transmute(loeuf.x = 1:nrow(ess_dat_under)/nrow(ess_dat_under), ess.y = cumsum(essential)/sum(ess_dat_under$essential))
lines(ess.y ~ loeuf.x, data = ess_dat_loeuf, col=3, lwd = 2)
ess_dat_am <- ess_dat_under %>% arrange(desc(mean_am_pathogenicity)) %>% transmute(am.x = 1:nrow(ess_dat_under)/nrow(ess_dat_under), ess.y = cumsum(essential)/sum(ess_dat_under$essential))
lines(ess.y ~ am.x, data = ess_dat_am, col=4, lwd = 2)
legend("topleft", c("alpha missense", "powerSFS", "LOEUF 2"), lwd = 2, col = c(4, 2, 3))

# Clinvar pathogenic de novo  
# 
denovo_patho <- read_tsv("data/clinvar_denovo_pathogenic.tsv")
dn <- denovo_patho %>% separate_longer_delim(`Gene(s)`, delim="|") %>%
  transmute(gene = `Gene(s)`) %>% group_by(gene) %>% summarise(n_patho=n()) %>%
  inner_join(dat)

plot(y ~ x, data = dn %>% arrange(desc(beta)) %>%
       transmute(x = 1:nrow(dn)/nrow(dn),
                 y = cumsum(n_patho)/sum(n_patho)),
     type="l", col = 2, lwd = 2,
     xlab  = "Quantile of intolerance metric", ylab = "proportion of denovo pathogenic ClinVar variants")
abline(a=0, b=1, col=1)
lines(y ~ x, data = dn %>% arrange(lof.oe_ci.upper) %>%
        transmute(x = 1:nrow(dn)/nrow(dn),
                  y = cumsum(n_patho)/sum(n_patho)), col=3, lwd = 2)
lines(y ~ x, data = dn %>% arrange(desc(mean_am_pathogenicity)) %>%
        transmute(x = 1:nrow(dn)/nrow(dn),
                  y = cumsum(n_patho)/sum(n_patho)), col=4, lwd = 2)
legend("topleft", c("LOEUF 2", "powerSFS", "alpha missense"), lwd = 2, col = c(3, 2, 4), bty = "n")

#
# de novo pathogenic and likely pathogenic
#

denovo_patho <- read_tsv("data/clinvar_denovo_pathogenic_likelypatho.tsv")
dn <- denovo_patho %>% separate_longer_delim(`Gene(s)`, delim="|") %>%
  transmute(gene = `Gene(s)`) %>% group_by(gene) %>% summarise(n_patho=n()) %>%
  inner_join(dat)

plot(y ~ x, data = dn %>% arrange(desc(beta)) %>%
       transmute(x = 1:nrow(dn)/nrow(dn),
                 y = cumsum(n_patho)/sum(n_patho)),
     type="l", col = 2, lwd = 2,
     xlab  = "Quantile of intolerance metric", ylab = "proportion of denovo pathogenic ClinVar variants")
abline(a=0, b=1, col=1)
lines(y ~ x, data = dn %>% arrange(lof.oe_ci.upper) %>%
        transmute(x = 1:nrow(dn)/nrow(dn),
                  y = cumsum(n_patho)/sum(n_patho)), col=3, lwd = 2)
lines(y ~ x, data = dn %>% arrange(desc(mean_am_pathogenicity)) %>%
        transmute(x = 1:nrow(dn)/nrow(dn),
                  y = cumsum(n_patho)/sum(n_patho)), col=4, lwd = 2)
legend("topleft", c("LOEUF 2", "powerSFS", "alpha missense"), lwd = 2, col = c(3, 2, 4), bty = "n")

#
# all patho and likely patho
#

denovo_patho <- read_tsv("data/clinvar_pathogenic_likelypatho.tsv")
dn <- denovo_patho %>% separate_longer_delim(`Gene(s)`, delim="|") %>%
  transmute(gene = `Gene(s)`) %>% group_by(gene) %>% summarise(n_patho=n()) %>%
  inner_join(dat)

plot(y ~ x, data = dn %>% arrange(desc(beta)) %>%
       transmute(x = 1:nrow(dn)/nrow(dn),
                 y = cumsum(n_patho)/sum(n_patho)),
     type="l", col = 2, lwd = 2,
     xlab  = "Quantile of intolerance metric", ylab = "proportion of denovo pathogenic ClinVar variants")
abline(a=0, b=1, col=1)
lines(y ~ x, data = dn %>% arrange(lof.oe_ci.upper) %>%
        transmute(x = 1:nrow(dn)/nrow(dn),
                  y = cumsum(n_patho)/sum(n_patho)), col=3, lwd = 2)
lines(y ~ x, data = dn %>% arrange(desc(mean_am_pathogenicity)) %>%
        transmute(x = 1:nrow(dn)/nrow(dn),
                  y = cumsum(n_patho)/sum(n_patho)), col=4, lwd = 2)
legend("topleft", c("LOEUF 2", "powerSFS", "alpha missense"), lwd = 2, col = c(3, 2, 4), bty = "n")

#
# all missense patho and likely patho
#

denovo_patho <- read_tsv("data/clinvar_missense_pathogenic_likelypatho.tsv")
dn <- denovo_patho %>% separate_longer_delim(`Gene(s)`, delim="|") %>%
  transmute(gene = `Gene(s)`) %>% group_by(gene) %>% summarise(n_patho=n()) %>%
  inner_join(dat)

plot(y ~ x, data = dn %>% arrange(desc(beta)) %>%
       transmute(x = 1:nrow(dn)/nrow(dn),
                 y = cumsum(n_patho)/sum(n_patho)),
     type="l", col = 2, lwd = 2,
     xlab  = "Quantile of intolerance metric", ylab = "proportion of denovo pathogenic ClinVar variants")
abline(a=0, b=1, col=1)
lines(y ~ x, data = dn %>% arrange(lof.oe_ci.upper) %>%
        transmute(x = 1:nrow(dn)/nrow(dn),
                  y = cumsum(n_patho)/sum(n_patho)), col=3, lwd = 2)
lines(y ~ x, data = dn %>% arrange(desc(mean_am_pathogenicity)) %>%
        transmute(x = 1:nrow(dn)/nrow(dn),
                  y = cumsum(n_patho)/sum(n_patho)), col=4, lwd = 2)
legend("topleft", c("LOEUF 2", "powerSFS", "alpha missense"), lwd = 2, col = c(3, 2, 4), bty = "n")

