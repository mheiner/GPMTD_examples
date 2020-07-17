
id = commandArgs(trailingOnly = TRUE)

library("coda")

### user input
simregexp = paste0(".*", id, ".*\\.rda")
dctory = "./postsim"
allfiles = list.files(dctory)
(selfiles = allfiles[grep(simregexp, allfiles)])
(nfiles = length(selfiles))

### end user input

i = 1
for ( i in 1:nfiles ) {
  load( paste0(dctory, "/", selfiles[i]) )

  nuse = 2000
  whichiter = floor(seq(1, nsim, length=nuse))

  pdf(paste0("plots/trace/", mesg, ".pdf"), height=4, width=9)

  traceplot(as.mcmc(sims_llik[whichiter]), main="log-likelihood")

  for (j in 1:(L+1)) {
    traceplot(as.mcmc(sims_lam[whichiter,j]), ylim=c(0,1),
              main=paste0("lambda ", j-1))
  }

  traceplot(as.mcmc(sims_int_mu[whichiter]), main="intercept mu")
  traceplot(as.mcmc(sims_int_sig2[whichiter]), main="intercept sig2")

  dev.off()

  cat(i, " of ", nfiles, "\r")
}

rm(list=ls())
