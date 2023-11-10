library(bayefdr)
postMCMC <- function(postParm){
  
  da.prob <- 1-postParm$posteriorProb.hj
  efdr <- efdr_search(da.prob, target_efdr = 0.05)
  bayes.fdr <- efdr$threshold[optimal(efdr)]
  # What is the bayesian FDR cutoff?
  plot(da.prob)
  DA.status.BayesFDR <- ifelse(da.prob>=bayes.fdr,2,1)
  DA.status.BayesFDR[1] <- 1 # make cahnegs to the root of this problem as the first taxa is always non-DA
  table(DA.status.BayesFDR)
  
  taxa.list.da.ZIBNP <- data.frame(taxa.name =postParm$data.sim$taxonomy ,
                                   da.status.ZIBNP = DA.status.BayesFDR[-1])
  return(taxa.list.da.ZIBNP)
}