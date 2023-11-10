rm(list=ls())

tt <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Import Data
#load("~/Desktop/Desktop september 2021/Everything on Desktop/OralExam_March2021/Final Documents/Future work/Simulation2/Sim2_ZIBNP/Data3/sparseD_simdata3.Rdata")

dat.name <- paste("simulate_data/sim2Data_",tt,".Rdata",sep="")

load(dat.name)

data <- data_sim2
 

options(warn=0)
library(mvtnorm)
library(DirichletReg)
#library(collapse)
library(gplots)
library("RColorBrewer")
library(tidyverse)
#library(RcppArmadillo)
library(invgamma)
library(LaplacesDemon) #-
library(mniw)  
library(ars)
library(gridExtra)
#library(grid)
library(ggplot2)
library(lattice)
library(ggfortify)
library(nlme)
library("factoextra")
library(gtools)
options(warn=2)
options(error=recover)




source("init_func_ZGP.R") 
source("real.dataInit.R")
source("mcmc0_ZGP.R")
source("clustAlloc_mcmc_new.R") # is this function even required?
source("update_c.R")
source("update_eta_newv3.R")
source("posteriorProb_hj.R")
source("postMCMC_bayesFDR.R")

data.type <- "real.data"
ref.clust <- FALSE #TRUE  # FALSE is to keep the ref taxa as 1 only
use.init.true.cj  <- FALSE
truncated.sig.e <- TRUE
range.R2 <-c(0.999,0.9999)
c.j.perc <- 1
data$c.j.perc <- c.j.perc
###### IMP!!!
# Adding the reference column of all 1's in the data


#data$Zmat <- cbind(taxa_ref=rep(1,data$n), data$Zmat)
data$Zmat.w0 <- data$Zmat
data$lib.size <- apply(data$Zmat,1,sum)
# since data$lib.size  will change based on S_i.tilde, store original lib.size in  lib.size.w0
data$lib.size.w0 <- data$lib.size
data$p <- ncol(data$Zmat)
dat <-data.frame(grp=data$grp,data$X1)
mmat <- model.matrix(grp~.,dat)
data$X <- mmat[,-1]
head(data$X)

table(data$grp)
data$grp <- as.numeric(as.factor(data$grp))
table(data$grp)
class((data$grp))
data$grp <- as.factor(data$grp)

# make a true parm object as well
true_parm <-  NULL
true_parm$K <- 2
true_parm$M <- 5 #14
true_parm$C <- 20
init.M <- true_parm$M

# define hyper parameters
hyper.param <- NULL
prior.mu <-NULL
tau = 1 #### why??? Does the following values make sense? ----
prior.mu$mean <- rep(0,ncol(data$X1)+1) # change ----
prior.mu$tau <- tau
x.intercept <- cbind(1,data$X) # change ----
prior.mu$xtx <- t(x.intercept) %*% x.intercept  # change ----
hyper.param$prior.mu <- prior.mu
hyper.param$a.e <- hyper.param$a.tau <- 0.001
hyper.param$b.e <- hyper.param$b.tau <-  0.001 
hyper.param$a.clus <- 1  # This is for c.j
hyper.param$a0 <- 1 # This is for beta.pi FMM
hyper.param$a0.int <- 0.001
hyper.param$b0.int <- 0.001
hyper.param$nu0.int <- 30
hyper.param$m0.int <- 0
data$hyper.param <- hyper.param
data$range.R2 <-range.R2 


# zero % 

mean(data$Zmat==0) * 100

niter= 150
nburnin=50

## ---------------  ----------##
# CASE 1: M and G are fixed , rest of them will be updated 
## ---------------  ----------##
############## ---------------
############## INITIALIZATION  #################
############## ---------------
#if use.init.true.cj=TRUE then cj are also fixed 
# Initialize 


data.init <- initialize.new.case.w0(data, true_parm,init.M)
data <- data.init

############## ---------------
# COPIED FROM fn_mcmc_ZGP.R , function : fn.mcmc.w0
init.parm <- data$init
# Collect all the parameters and the data (eta.hat) in a new list parm 
data$parm<- NULL
data$parm$beta.parm <- data$init$beta.parm
data$parm$h.j <- data$init$h.j
data$parm$sig2e <- data$init$sig2e
data$parm$tau <- data$init$tau
data$parm$c.j <- data$init$c.j
data$parm$m.u <- data$init$m.u
data$parm$eta.hat <- data$init$eta.hat
data$parm$qstar.hat.mat <- data$init$qstar.hat.mat 
data$parm$log.l <- c(rep(0,data$G))
#data$parm$beta0iu.est <- data$init$beta0iu.est # change ----
data$parm$Ji <- data$init$Ji
data$parm$no.nontech.0 <- data$init$no.nontech.0


# continue running fn_mcmc_ZGP.R
# MCMC iterations
data11 <- NULL
pprob.mat<- array(,c(niter,data$p))
data11[[1]] <- mcmc.update.case2_randintercept.w0(data,hyper.param)

# check post ----------------------
posteriorProb.hj.sum <- func.posteriorProb.hj(data1=data11[[1]])
plot(1-posteriorProb.hj.sum)
pprob.mat[1,] <- posteriorProb.hj.sum
#true.hj <- c(1,true_parm$h.j) #truth$h.j
#pred.hj <- ifelse((1-posteriorProb.hj.sum)>0.5,2,1)
#table(pred.hj,true.hj)
# check post ----------------------

posteriorProb.hj.sum <- 0 #func.posteriorProb.hj(b=1)
START.TIME <- Sys.time()
for(nn in 2:niter){
  print(paste("************ ------------ ************ -------------- ************   ", nn))
  data11[[nn]] <- mcmc.update.case2_randintercept.w0(data11[[nn-1]],hyper.param)
  # posteriorProb.hj.sum <- posteriorProb.hj.sum + func.posteriorProb.hj(b=nn)
  pprob.mat[nn,] <- func.posteriorProb.hj(data1=data11[[nn]])
  if(nn>nburnin){
    posteriorProb.hj.sum <- posteriorProb.hj.sum + pprob.mat[nn,]
  } 
  
}
END.TIME <- Sys.time()
total.time.mcmc <- END.TIME-START.TIME
total.time.mcmc
#out.write <- paste("Total time to run the MCMC for C=", data$C, " and M=",init.M," is : ",total.time.mcmc," mins")
#write.table(out.write,file="tot_time.txt")

posteriorProb.hj <- posteriorProb.hj.sum/(niter-nburnin)




posterior.sample <- NULL
posterior.sample$sigma2 <- posterior.sample$tau2 <- posterior.sample$tau0.int <- posterior.sample$m0.int <- vector()
posterior.sample$da.status.clust <- array(,c(niter,(data$G)))
posterior.sample$vku <- array(,c(niter,(true_parm$K*data$C)))
posterior.sample$mu.m.beta <- array(,c(true_parm$M,ncol(data$X)+1,niter)) # change ----
#posterior.sample$beta0iu <- array(,c(data$n,data$C,niter)) # change ----
posterior.sample$h.j <- array(,c(niter,(data$p)))
posterior.sample$m.u <- array(,c(niter,(data$G)))
posterior.sample$eta.mt <- posterior.sample$qstar <- array(,c(data$n,data$G,niter))
posterior.sample$c.j <-array(,c(niter,(data$p)))
posterior.sample$pi.mix <-array(,c(niter,(data$parm$beta.parm$M)))
posterior.sample$lib.size <- array(,c(niter,data$n))
posterior.sample$no.0 <- array(,c(niter))
posterior.sample$Si.tilde <- array(,c(niter,data$n))

for(n in 1:niter){
  posterior.sample$sigma2[n] <- data11[[n]]$parm$sig2e
  posterior.sample$tau2[n] <- data11[[n]]$parm$tau
  posterior.sample$tau0.int <- data11[[n]]$parm$tau0
  posterior.sample$m0.int <- data11[[n]]$parm$mu0
  posterior.sample$da.status.clust[n,] <- data11[[n]]$parm$beta.parm$da.status
  posterior.sample$vku[n,] <- data11[[n]]$parm$beta.parm$v.ku
  posterior.sample$h.j[n,] <- data11[[n]]$parm$h.j
  posterior.sample$mu.m.beta[,,n] <- data11[[n]]$parm$beta.parm$mu.m
  posterior.sample$m.u[n,] <-data11[[n]]$parm$m.u # this is actually fixed since we do not change c.j 
  posterior.sample$eta.mt[,,n] <-data11[[n]]$parm$eta.hat
  posterior.sample$qstar[,,n] <-data11[[n]]$parm$qstar.hat.mat
  posterior.sample$c.j[n,] <-data11[[n]]$parm$c.j
  # posterior.sample$beta0iu[,,n] <- data11[[n]]$parm$beta0iu.est # change ----
  posterior.sample$pi.mix[n,] <-data11[[n]]$parm$beta.parm$pi
  posterior.sample$lib.size[n,] <- data11[[n]]$lib.size
  posterior.sample$no.0[n] <- data11[[n]]$parm$no.0
  posterior.sample$Si.tilde[n,] <- data11[[n]]$parm$Si.tilde
}


posterior.sample$pprob.hj <- pprob.mat

postParm <- NULL
postParm$posterior.sample <- posterior.sample
postParm$posteriorProb.hj <- posteriorProb.hj
postParm$true_parm <- true_parm
postParm$data.sim <- data
postParm$ref.z.ind <- data$init$ref.z.ind
postParm$init.parm <- init.parm

taxa.list.da.ZIBNP <- postMCMC(postParm)

filename <- paste("postSample_sim2_self_",tt,".Rdata",sep="")
save(data.init, postParm,true_parm,file=filename)

save(taxa.list.da.ZIBNP, file=paste("ZIBNPresults_sim2_self_dat",tt,".Rdata",sep="")) 
