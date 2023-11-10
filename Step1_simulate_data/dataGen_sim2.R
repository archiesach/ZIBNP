rm(list=ls())

tt <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# IMPORTANT PARAMETERS FOR SIMULATION 2 - SELF

R2.aij <- 0.97
pi0 <- 0.3
psi0 <- 6
######## ---- % of zeros
lambda0 <- -0.1
# 0.02 # 58% zero
  #-0.1 # 15% zero
######## ---- Fold change
fold.change0 <- 2



options(warn=0)
library(mvtnorm)
library(DirichletReg)
library(gplots)
library("RColorBrewer")
library(tidyverse)
library(invgamma)
library(LaplacesDemon) #-
library(mniw)  
library(ars)
library(gridExtra)
library(ggplot2)
library(lattice)
library(ggfortify)
library(nlme)
library("factoextra")
library(gtools)
library(rpartitions)
options(warn=2)
options(error=recover)


source("func_sim_v2.R")
source("init_func_v1.R")
source("ZGP.R")
source("sim2_self.R")


# setwd("~/Desktop/Desktop september 2021/Everything on Desktop/OralExam_March2021/Final Documents/Future work/model_simulation_ZeroGenP/ZIBNP_allcodes/Data_simulation")

###### -------- datagenRun_ZGP.R ---- start
######  ----------- ##########
#  USER INPUT
######  ----------- ##########

p.cov = 1000
K = 2
M <- 5 #7 # 8  # 7 for seed = 1125
n <- 100

tau = 1
range.R2 <-c(0.999,0.9999)
R2.eta <- 0.999
alpha.p <- 1
alpha.cj <- 2 #10
perc <- 1 # % of c.j's to update
status.da.prob <- c(0.003,0.474,0.523) # c(ref,non-DA,DA) # c(0.003,0.474,0.523) #0.003 0.374 0.623
A <-100000 #100000

# Other indicators

# these two are initialization flags 
clust.method = "nonDP"  #"DP", nonDP
data.type <- "sim.data"
use.init.true.cj  <- FALSE
truncated.sig.e <- TRUE

need.pca<- FALSE
plot.eta <- FALSE
with.intercept <- FALSE
do.subjects.same= TRUE
type.mu=  "zellner" #"zellner_diffmean" 


# Data
pheno.all <- read.csv("data/covariates_1.csv")
lib.size <- read.csv("data/libSize.csv")
all.equal(lib.size$sample.id,pheno.all$X)

li.factor <- 1
libsize= "rsamp" #"poisson2"  # "rsamp", "poisson1", "poisson2" (poisson 2 was how we were generating the lib.size earlier)
pois.lambda <- 10000  # for poisson1
max.nonDP  <- 12#12
min.nonDP <- 8#8
######  ----------- ##########
# Generate covariate matrix

set.seed(tt+100)

source("genXcode_v3.R")
hist(data$lib.size)
range(lib.size$L_i)
range(data$lib.size)

# To check if the covariates have been standardized
apply(data$X,2,mean)
apply(data$X,2,sd) 

if(with.intercept ==TRUE){
  data$X <- cbind(1,data$X)
}        


data$c.j.perc <- perc
data$p <- p.cov +1
data$n <- nrow(data$X)
data$tau <- tau

# Generate parameters, c.j, mu.m, pi, beta

true_parm <- NULL
true_parm$p <- p.cov
true_parm$K <-K
true_parm$R2.eta <- R2.eta
true_parm$clust <- NULL
true_parm$clust$a1 <- alpha.cj  # PDP alpha # move it to hiper prior ----

if(clust.method =="DP"){
  true_parm <- clust.alloc(true_parm)
}else if(clust.method =="nonDP"){
  true_parm <- clust.alloc.nonDP(true_parm)
}



###### -------- datagenRun_ZGP.R ---- end

XX <- data$X
dim(XX)

cj <- true_parm$clust$c.v
C <- true_parm$clust$C

true_parm$da.status <- da.status.gen(pi0,C)


#check

table(true_parm$da.status )

# da status for all taxa
true_parm <- da.status.true.parm(true_parm)

psi.vec <- psi(psi0 = psi0,fold.change =fold.change0 )
bet <- gen.beta.norm(C,tau,XX)
beta.ku <- beta.ku(psi.vec,bet,K, true_parm$da.status,C)
p= data$p
cj<- cj[-1] -1
#lam0=30000


data.sim2 <- counts.ij(XX,beta.ku,grp,p,cj)

data.sim2$p <- data$p
data.sim2$n <- data$n
data.sim2$X1 <- data$X1
data.sim2$X <- data$X
data.sim2$grp <- data$grp
data.sim2$true_parm <- true_parm
data.sim2$lib.size <- apply(data.sim2$Zij,1,sum)
data.sim2$Zmat <- cbind(1,data.sim2$Zij)
data.sim2$taxonomy <- paste("taxa_",1:(ncol(data.sim2$Zmat)-1),sep="")
#View(data.sim2)

# Zero-inflation



data.w0 <- fn.zero.gen(data=data.sim2,lambda0)
hist(data.w0$lib.size.w0)
hist(data.w0$lib.size.wout0)
# Percentage of zeros


fname <- paste("sim2Data_",tt,".Rdata",sep="")

data_sim2 <- NULL
data_sim2$Zmat <- data.w0$Zmat
data_sim2$X1 <- data.w0$X1
data_sim2$grp <- data.w0$grp
data_sim2$n <- data.w0$n
data_sim2$taxonomy <- data.sim2$taxonomy

true.DA <- true_parm$h.j[-1]
truth <- true_parm
save(fold.change0,lambda0,data_sim2,true.DA,truth,file=fname)
#save(data_sim2,true.DA,truth,file=fname)
