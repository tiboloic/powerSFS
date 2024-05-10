# LT 9/05/2024
# 
# Preston plots
# 
# 

psfs <- read_tsv("gnomADv4_protein_altering_sfs_prototype_fit.tsv.gz")
names(psfs)[1] <- "gene"

sfs <- read_tsv("gnomADv4_protein_altering_sfs_prototype.tsv.gz")
obs <- sfs %>% 
  select(gene, octave, count)  %>% 
  spread(key=octave, value=count, fill=0) %>%
  arrange(gene)

# calculate average gene ps:
beta <- psfs$mub[1]
pred = rep(0,21)
for (i in 0:20) {
  pred[i+1] = sum((2^i:(2^(i + 1) - 1)) ^ -beta)
}
ps <- pred / sum(pred)

preston_plot <- function(gene_name) {
  one_gene_sfs <- unlist(filter(obs, gene == gene_name)[-1])
  this_gene <- filter(psfs, gene == gene_name)
  if (nrow(this_gene) == 0) stop("gene not found")
  beta <- this_gene$beta
  pred = rep(0,21)
  for (i in 0:20) {
    pred[i+1] = sum((2^i:(2^(i + 1) - 1)) ^ -beta)
  }
  pred <- pred / sum(pred) * sum(one_gene_sfs)
  plot(one_gene_sfs, log="y", main = gene_name, xlab = "octave", ylab = "missense count", xaxt = "n")
  lines(pred, lwd = 2, col = 2)
  lines( ps * sum(one_gene_sfs))
  axis(1, at = 1:21, las = 1, labels = 0:20)
  pval <- ifelse(this_gene$p < 0.000001, "< 0.000001", round(this_gene$p, 6))
  text(20, one_gene_sfs[2], substitute(paste(italic("p = "), pval)))
}

preston_plot("MIB1")
preston_plot("ACTA1")
preston_plot("DNMT3A")
preston_plot("SF3B2")
preston_plot("OR10G4")
preston_plot("A1BG")
