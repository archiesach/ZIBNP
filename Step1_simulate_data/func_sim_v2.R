

which.m <- function(x){
  which(x==1)
}


mixture.sim.zellner <- function(M,alpha.p,X,type.mu,tau){
  
  T1 <- ncol(X)
  ll <- floor((M/2))
  Aa <- 0.6*seq(-ll,ll,by=1)
  
  xtx <- t(X) %*% X
  ss <- tau * solve(xtx)
  mu <- matrix(0,nrow = M,ncol = T1)
  for( i in 1:M){
    if(type.mu=="zellner"){
      mu[i,] <- rmvnorm(n=1,mean = rep(0,T1),sigma = ss) 
    }
    if(type.mu=="zellner_diffmean"){
      mu[i,] <- rmvnorm(n=1,mean = rep(Aa[i],T1),sigma = ss) 
    }
    
  }
  prior.mu <-NULL
  prior.mu$mean <- rep(0,T1)
  prior.mu$tau <- tau
  prior.mu$xtx <- xtx
  a.vec <- rep(alpha.p/M,M)
  #ed(511115)
  pi <- rep(1/M,M)  # rdirichlet(n=1, a.vec)  #
  return(list(mu=mu,pi=pi,prior.mu=prior.mu))
  
}


beta_gen <- function(K,M,T1,mu.beta,pi.beta,alpha.p){
  mu <- mu.beta   # this is a M x T1 matrix
  pi <- pi.beta
  m <- matrix(0,nrow=M, ncol=K)
  for(i in 1:K){
    m[,i] <- rmultinom(n=1,size=1,prob = pi)
  }
  comp.m <- apply(m , 2, which.m)
  return(list(comp.m=comp.m,beta=mu[comp.m,]))
}

clust.alloc.nonDP <- function(true_parm){
  library(rpartitions)
  #divide_and_conquer(c(5, 4), 5, 4, hash(), 2, TRUE, FALSE)
  p <- true_parm$p
  rand.part <- round(runif(n=1,min = min.nonDP,max = max.nonDP),0)
  c.m.vec <- rand_partitions(p, rand.part, 1)
  rand.ind <- sample(c(1:p),size = p,replace = F)
  stack <- rep(c(1:rand.part),times=c.m.vec)
  
  rand.cj <- stack[rand.ind]
  print(all.equal(as.vector(c(table(rand.cj))),c(as.vector(c(table(rand.cj))))))
  
  true_parm$clust$c.v <- c(1, (rand.cj+1) )
  true_parm$clust$C.m.vec <- c(1,as.vector(c.m.vec))
  true_parm$clust$G <- rand.part+1
  true_parm$clust$C <- rand.part
  true_parm$hyperparam$a0 <- 1 # hyper prior paramater for pi in FMM
  true_parm
  
}

clust.alloc <- function(true_parm){
  #set.seed(1251)
  #set.seed(1125)
  #set.seed(751)
  #set.seed(51)
  p <- true_parm$p
  true_parm$clust$c.v <- array(,p)
  true_parm$clust$c.v[1] <- 1   # c1, c2, ... , cp
  true_parm$clust$C.m.vec <- 1  # n1, n2, ... , nq
  true_parm$clust$G <- 1 # q
  
  for (xx in 2:p)
  {
    prob.v <- true_parm$clust$C.m.vec 
    prob.v <- c(prob.v, true_parm$clust$a1 )
    new.c <- sample(1:(true_parm$clust$G+1), size=1, prob=prob.v)
    
    true_parm$clust$c.v[xx] <- new.c
    new.flag <- (new.c > true_parm$clust$G)
    
    if (new.flag)
    {true_parm$clust$G <- true_parm$clust$G + 1
    true_parm$clust$C.m.vec <- c(true_parm$clust$C.m.vec, 1)
    }
    if (!new.flag)
    {true_parm$clust$C.m.vec[new.c] <- true_parm$clust$C.m.vec[new.c] + 1
    }
  }
  true_parm$clust$c.v <- c(1, (true_parm$clust$c.v+1) )
  true_parm$clust$C.m.vec <- c(1,true_parm$clust$C.m.vec)
  true_parm$clust$G <- true_parm$clust$G+1
  true_parm$clust$C <- true_parm$clust$G-1
  true_parm$hyperparam$a0 <- 1 # hyper prior paramater for pi in FMM 
  true_parm
}


r.sq.lm <-function(x,data){
  tmp = lm( x ~ -1 + data$X)
  tmp2 = summary(tmp)
  R.sq = tmp2$r.squared
  return(R.sq)
}



gen.eta.new <- function(true_parm,data){
  eta <- array(,c(true_parm$n, true_parm$clust$C))
  for(i in 1:true_parm$n){
    for(c in 1:true_parm$clust$C){
      if(data$grp[i]==1){
        eta[i,c] = true_parm$beta[1,,c] %*% data$X[i,]
      }else if(data$grp[i]==0){
        eta[i,c] = true_parm$beta[2,,c] %*% data$X[i,]
      }
    }
  }
  
  sigm.e <- var(as.vector(eta)) *((1/true_parm$R2.eta) -1)
  data$sigm.e <- sigm.e
  eta1 <- eta +  matrix(rnorm(n=data$n*(true_parm$clust$C),sd=sqrt(data$sigm.e)), nrow=data$n)   #matrix(rnorm(n=data$n*(true_parm$clust$C),sd=sqrt(data$sigm.e)), nrow=data$n)  
  #R.sq <- apply(eta1,2, r.sq.lm,data)
  #R.sq <- apply(eta,2,function(x) var(x)/(var(x)+true_parm$sigm.e))
  data$eta.true <- eta
  data$eta <- cbind(0,eta1)
  #data$R.sq <- R.sq
  data$R.sq.overall <- var(as.vector(eta))/var(as.vector(eta1))
  return(data)
}


check.simplex.sim <- function(true_parm,data){
  # qstar
  m.u <- true_parm$clust$C.m.vec
  check <- apply(apply(data$qstar ,1, function(x) m.u*x),2,sum)
  #print(plot(check,main="sim"))
}


check.simplex.mcmc <- function(data){
  # qstar
  m.u <- data$parm$m.u
  check <- apply(apply(data$parm$qstar.hat.mat ,1, function(x) m.u*x),2,sum)
  #print(plot(check,main="mcmc"))
}


q.star.new <- function(data,true_parm){
  which_nD = setdiff(which(true_parm$da.status==1),1)   #setdiff(which(true_parm$h.j==1), 1)
  which_D = which(true_parm$da.status==2)
  
  # find alpha
  prob_1.v = true_parm$rho[,1]
  prob_nD.v = true_parm$rho[,2]
  prob_D.v = true_parm$rho[,3]
  alpha.mt = array(,c(n,2))
  m.v <- true_parm$clust$C.m.vec
  alpha.mt[,1] = prob_nD.v/as.vector((exp(data$eta[,which_nD,drop=FALSE])%*% m.v[which_nD]))
  alpha.mt[,2] = prob_D.v/as.vector((exp(data$eta[,which_D,drop=FALSE])%*%m.v[which_D]))
  
  qstar <- array(,c(true_parm$n,(true_parm$clust$G)))
  qstar[,1] <- prob_1.v
  qstar[,which_nD] = alpha.mt[,1]*exp(data$eta[,which_nD])
  qstar[,which_D] = alpha.mt[,2]*exp(data$eta[,which_D])  
  data$alpha.mt <- alpha.mt
  data$qstar <- qstar
  return(data)
}




gen.rho <- function(status.da.prob,A,n){
  rho.i <- rdirichlet(n = n,c(status.da.prob*A))
  return(rho.i)
}

# get number of clusters from c, DP_cluster alloc function

trueParm <- function(){
  true_parm <- NULL
  true_parm$p <- p.cov
  true_parm$K <-K
  true_parm$R2.eta <- R2.eta
  true_parm$clust <- NULL
  true_parm$clust$a1 <- alpha.cj  # PDP alpha # move it to hiper prior ----
  true_parm <- clust.alloc(true_parm)
  true_parm$mix.comp <- mix.comp
  # generate cluster allocation vector
  
  T1 <- ncol(data$X)
  all.beta <- replicate(true_parm$clust$C,beta_gen(K=K,M,T1,mu.beta,pi.beta,alpha.p))
  true_parm$beta <- array(,c( true_parm$K,T1,true_parm$clust$C))
  true_parm$beta <-array(as.numeric(unlist(all.beta[2*(1:true_parm$clust$C)])), dim=c(true_parm$K,T1,true_parm$clust$C))
  true_parm$m_beta <- all.beta[(2*(0:(true_parm$clust$C-1)))+1]
  true_parm$da.status <- c(1,as.numeric(ifelse(lapply(true_parm$m_beta,function(x) length(unique(x)))==1,1,2)))
  true_parm$n <- nrow(data$X)
  true_parm$M <-M
  true_parm$v.ku.mat <- matrix(unlist(true_parm$m_beta),nrow=2)
  
  return(true_parm)
}

beta_gen.0 <- function(M,K,T1,C,mu){
  mu <- mu.beta   # this is a M x T1 matrix
  combn <- permutations(n=M,r=K,repeats.allowed = T)
  sss1 <- which(combn[,1]==combn[,2])
  sss2 <- setdiff(c(1:nrow(combn)),sss1)
  c1 <- round(0.2*C,0)
  c2 <- C-c1
  index.sam1 <- sample(sss1,size =c1,replace = FALSE,prob = NULL)
  index.sam2 <- sample(sss2,size =c2,replace = FALSE,prob = NULL)
  index.sam0 <- c(index.sam1,index.sam2)
  index.sam <- sample(index.sam0,size=length(index.sam0),replace = FALSE,prob = NULL) # reshuffle again
  comp.m <- t(combn[index.sam,])
  
  beta <- array(,c( K,T1,C))
  for(CC in 1:C){
    beta[,,CC] <- mu[comp.m[,CC],]
  }
  return(list(comp.m=comp.m,beta=beta))
}

trueParm.0 <- function(){
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
  
  
  
  true_parm$mix.comp <- mix.comp
  # generate cluster allocation vector
  true_parm$M <-M
  T1 <- ncol(data$X)
  
  all.beta <- beta_gen.0(M=true_parm$M,K=true_parm$K,T1=T1,C=true_parm$clust$C, mu=mu.beta) #replicate(true_parm$clust$C,beta_gen(K=K,M,T1,mu.beta,pi.beta,alpha.p))
  true_parm$beta <- all.beta$beta
  true_parm$m_beta <- all.beta$comp.m
  true_parm$da.status <- c(1,apply(true_parm$m_beta,2, function(x) length(unique(x))  ))
  true_parm$n <- nrow(data$X)
  true_parm$v.ku.mat <- matrix(unlist(true_parm$m_beta),nrow=2)
  
  return(true_parm)
}


qfromqstar <- function(true_parm,data){
  q.mat <- matrix(0,nrow = nrow(data$qstar),ncol = (data$p))
  q.mat[,1] <- data$qstar[,1]
  for(u in 2: true_parm$clust$G ){
    indx <- which(true_parm$clust$c.v==(u))
    q.mat[,indx] <- data$qstar[,u]
  }
  data$qmat <- q.mat
  return(data)
}



Z.abundance <- function(true_parm,data){
  Zmat <- matrix(0,nrow = (true_parm$n),ncol = (data$p))
  for(i in 1:true_parm$n){
    Zmat[i,] = rmultinom(n=1,size = data$lib.size[i],prob = data$qmat[i,])
  }
  colnames(Zmat) <-c(1:(data$p))
  data$Zmat <- Zmat
  return(data)
}

