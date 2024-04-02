# LT 19/08/2020
# use variable allele number

# LT 3/06/2020

# refit model using variant list from gene wich HGNC id
# use refitted pade approximation

# LT 13/12/2019

#
# TMB multinomial power law with random slope

# 18/12/2019
# padé approximation of power porbabilities

library(TMB)
library(dplyr)
library(tidyr)


# for debugging
fixwinpath <- function() {
  PATH <- Sys.getenv("PATH")
  PATH <- paste0(R.home(), "/bin/x64;", PATH)
  PATH <- paste0("C:/RBuildTools/3.5/mingw_64/bin;", PATH)
  Sys.setenv(PATH=PATH)
}

gnomad = read.delim(file='gnomADSFS.4.tsv')

# remove an
SFS = gnomad[,-4] %>% spread(key=cat, value=obs, fill=0)

# get AN: group by gene, get sum of an and total number of variants
AN = gnomad %>% group_by(gene) %>% summarize(variant_count=sum(obs), an=sum(an))
AN = gnomad %>% group_by(gene) %>% summarize(an=sum(an)/sum(obs))

any(AN$an<16)
which(AN$an<16)
SFS = SFS[AN$an>16,]
AN=AN[AN$an>16,]

# remove gene information: keep only counts
#dat = as.matrix(SFS[,-1, drop=FALSE])
dat = as.matrix(SFS[,-1, drop=FALSE])

#compile("powerSFSan.cpp", "-O1 -g",DLLFLAGS="")
compile("powerSFSan.cpp")

dyn.load(dynlib("powerSFSan"))

f = MakeADFun(
  data=list(obs=dat, an=as.vector(AN$an, "integer")),
  parameters=list(b=1.8, lsigb=-2, bdevs=rep(0,nrow(dat))),
  random=c("bdevs"),
  DLL="powerSFSan", silent=FALSE)


#fit = optim(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))
fit = nlminb(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=5000))
mysd = sdreport(f)

SFS$beta = f$env$parList()$bdevs
SFS$an = AN$an
save(f, SFS, mysd, file='fit.AN.RData')
