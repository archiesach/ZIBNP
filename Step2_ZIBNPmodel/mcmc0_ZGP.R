
ssquare <- function(y,x,bet){  # ****This function can be written in Rcpp ----
  ssq <- 0
  x <- cbind(1,x) # change ----
  for(i in 1:ncol(y)){
    ssq <- ssq + sum(( y[,i] - x %*% bet[i,] )^2)
  }
  return(ssq)
}

ssquare.new <- function(y,x,bet){  # ****This function can be written in Rcpp ----
  ssq <- 0
  for(i in 1:ncol(y)){
    ssq <- ssq + sum(( y[,i] - x %*% bet )^2)
  }
  return(ssq)
}


range_sigmasq <- function(ll,ul,y){
  #Rsq <- 1- s^2/ var(y)
  range <- NULL
  range$sig.ll <- ((1-ul)*var(y))
  range$sig.ul <- ((1-ll)*var(y))
  return(range)
}


truncated_InvGamma3 <- function(astar,bstar,range.sigma,flag.unif){
  ll <- pinvgamma(range.sigma$sig.ll,shape = astar,scale=bstar)
  ul <- pinvgamma(range.sigma$sig.ul,shape = astar,scale=bstar)
  if(round(ul,10)==0){
    s <- runif(1,range.sigma$sig.ll,range.sigma$sig.ul)
    flag.unif <- flag.unif +1
  } else if(round(ll,10)==1){
    s <- runif(1,range.sigma$sig.ll,range.sigma$sig.ul)
    flag.unif <- flag.unif +1
  } else {
    u <- runif(1,ll,ul)
    s <- qinvgamma(u,shape = astar,scale = bstar)
  }
  
  return(list(ssq=s,ll=ll,ul=ul,flag.unif=flag.unif))
}


update.sigma.e0 <- function(data,hyper.param){
  
  ###### sigma.e 
  len <- length(data$parm$beta.parm$v.ku)
  eta.hat.C <- data$parm$eta.hat[,-1]
  vv <- matrix(data$parm$beta.parm$v.ku,nrow=true_parm$K)
  # grp k beta 
  # eval(parse(text= paste("yy",1:true_parm$K,"<- (eta.hat.C[data$grp==",1:true_parm$K,",]-data$parm$beta0iu.est[data$grp==",1:true_parm$K,",])",sep="" )  ))
  eval(parse(text= paste("yy",1:true_parm$K,"<- (eta.hat.C[data$grp==",1:true_parm$K,",])",sep="" )  )) # change ----
  
  bet <- NULL
  bet <- lapply(1:true_parm$K, function(x) data$parm$beta.parm$mu.m[(vv[x,]),] )
  
  eval(parse(text= paste("s",1:true_parm$K,"=ssquare(y=yy",1:true_parm$K,", x=data$X[data$grp==",1:true_parm$K,
                         ",] , bet=bet[[",1:true_parm$K,"]])",sep="")    ))
  
  sum.sq.sum = eval(parse(text=  paste(paste("s",1:true_parm$K,sep=""),collapse="+")   ))
  
  a.post = hyper.param$a.e + (data$C*data$n)/2
  b.post = hyper.param$b.e  + sum.sq.sum/2     # * generate sigma.e # beta$parm$sig2e ----
  
  # sample sigma.e sq directly from InvGamma - CASE 1 
  if(truncated.sig.e){
    # sample sigma.e sq from truncated InvGamma - CASE 2 
    # NEW truncated dist (start) ----
    y.1 <- as.vector(eta.hat.C)
    range_sig.sq <- range_sigmasq(ll=data$range.R2[1],ul=data$range.R2[2],y=y.1)
    s.out <- truncated_InvGamma3(astar=a.post,bstar=b.post,range.sigma=range_sig.sq,flag.unif=0)
    data$parm$sig2e <- s.out$ssq
    # NEW truncated dist (end)----
  }else{
    s.out <- rinvgamma(n=1,shape = a.post,scale = b.post)
    data$parm$sig2e <- s.out
  }
  return(data)
}


update.v_ku <- function(pi.post){
  #print(pi.post)
  pi.post1 <- t(apply(pi.post,1, function(x)  exp(x-max(x))))
  pi.sum <- apply(pi.post1,1,sum)
  pi.post1 <- pi.post1/pi.sum
  # replace zeros by "small"
  small2 <- 1e-5
  pi.post1[pi.post1 < small2] <- small2
  # again normalize
  pi.sum1 <- apply(pi.post1,1,sum)
  pi.post <- pi.post1/pi.sum1   
  v.ku <- vector()
  for(i in 1:nrow(pi.post)){
    v__ku <- rmultinom(n= 1, size=1,prob = pi.post[i,] )  
    v.ku[i]  <- apply(v__ku,2, function(x) which(x==1) )
  }
  return(list(v.ku=v.ku,pi.post=pi.post))
}

update.beta.FMM.new0 <- function(data){
  # This function must update elements of data$parm$beta.parm  except mu.m
  len <- length(data$parm$beta.parm$v.ku)
  v.ku.mat <- matrix(data$parm$beta.parm$v.ku,nrow=true_parm$K)
  #grp <- ifelse(data$grp==1,1,2)
  xx <- NULL
  
  # xx[[1]] <- data$X[grp==1,]
  # xx[[2]] <- data$X[grp==2,]
  eval(parse(text= paste("xx[[",1:true_parm$K,"]] <- data$X[data$grp==",1:true_parm$K,",]",sep="")   ))
  log.lik.bet <- NULL
  pi.post <-pi.post.prob <- NULL
  #log.lik.bet[[1]] <- log.lik.bet[[2]] <- array(,c(data$C,data$parm$beta.parm$M))
  eval(parse(text = paste( "log.lik.bet[[", 1:true_parm$K,"]] <- array(,c(data$C,data$parm$beta.parm$M))",sep=""   )    ))
  
  # for(k in 1:2){
  #   for(u in 2:data$G){
  #     y_kiu0 = matrix(data$parm$eta.hat[grp==k,u],ncol=1)
  #     y_kiu = y_kiu0 - data$parm$beta0iu.est[grp==k,(u-1)]
  #     
  #     for(m in 1:data$parm$beta.parm$M){
  #       bet <- matrix(data$parm$beta.parm$mu.m[m,],nrow=ncol(data$X))
  #       log.lik.bet[[k]][u-1,m] = (-1/(2*data$parm$sig2e))*ssquare.new(y=y_kiu,x=xx[[k]],bet=bet)
  #     }
  #   }
  #   pi.post[[k]] <-  log(data$parm$beta.parm$pi)  + log.lik.bet[[k]]
  # }
  grp <- data$grp
  
  for(k in 1:true_parm$K){
    for(u in 2:data$G){
      y_kiu0 = matrix(data$parm$eta.hat[grp==k,u],ncol=1)
      y_kiu = y_kiu0 # y_kiu0   - data$parm$beta0iu.est[grp==k,(u-1)] # change ----
      
      for(m in 1:data$parm$beta.parm$M){
        bet <- matrix(data$parm$beta.parm$mu.m[m,],nrow=(ncol(data$X)+1)) # change ----
        log.lik.bet[[k]][u-1,m] = (-1/(2*data$parm$sig2e))*ssquare.new(y=y_kiu,x=cbind(1,xx[[k]]),bet=bet) # change  cbind(1,x) ----
      }
    }
    pi.post[[k]] <-  log(data$parm$beta.parm$pi)  + log.lik.bet[[k]]
  }
  
  vv <- array(,c(true_parm$K,data$C))
  # sample v_ku
  for(k in 1:true_parm$K){
    updatevku.out <- update.v_ku(pi.post=pi.post[[k]] )
    vv[k,] <- updatevku.out$v.ku  
    pi.post.prob[[k]] <- updatevku.out$pi.post
  }
  
  data$parm$beta.parm$v.ku <- c(vv)
  data$parm$beta.parm$d.m <- tabulate(data$parm$beta.parm$v.ku,nbins = data$parm$beta.parm$M)  
  post.alpha <- ((data$hyper.param$a0)/data$parm$beta.parm$M)+ data$parm$beta.parm$d.m
  data$parm$beta.parm$pi <- as.vector(rdirichlet(n=1,alpha = post.alpha))
  data$parm$beta.parm$beta.mu <- data$parm$beta.parm$mu.m[data$parm$beta.parm$v.ku,] 
  data$parm$beta.parm$da.status <- c(1,ifelse(apply(matrix(data$parm$beta.parm$v.ku,nrow=true_parm$K),2,function(x)  sum(length(unique(x))==1)  )==1,1,2))
    #c(1, apply(matrix(data$parm$beta.parm$v.ku,nrow=true_parm$K),2, function(x) length(unique(x)))  )
  data$parm$check$pi.post.vku <- pi.post.prob
  return(data)
}


# posterior.mu <- function(y,x,data){ # ****This function can be written in Rcpp ----
#   xtx <- t(x) %*% x
#   xtx.inv <- solve(xtx)
#   mu.hat.m <- xtx.inv  %*% t(x) %*% y
#   W <- solve( (xtx) + (solve(true_parm$mix.comp$prior.mu$xtx)* (data$parm$sig2e/ data$tau )) ) %*% xtx
#   A.star <- W %*% mu.hat.m 
#   B.star <- data$parm$sig2e * (W %*% xtx.inv)
#   
#   mu.m <- rmvnorm(n=1,mean = A.star,sigma = B.star)
#   return(mu.m)
# }

redesign.x <- function(x){
  x <- data.frame(x)
  degen <- which(apply(x,2,function(i)  length(unique(i)))==1)
  ind <- names(degen)
  x.new <- x %>% dplyr::select(-all_of(ind))
  x.new <- as.matrix(x.new)
  redesign.flag <- (length(ind)>0)
  x.des <-NULL
  x.des$x <- x.new
  x.des$redesign.flag <- redesign.flag
  x.des$degen <- degen
  return(x.des)
}


posterior.mu.new <- function(y,x,data){ # ****This function can be written in Rcpp ----
  # Before this remove the degenerate variables from the X matrix to avoid singularity issues in the inverse
  x.orig <- data.frame(x)
  x.redes <- redesign.x(x)
  x <- x.redes[["x"]]
  x <- cbind(1,x) # change ----
  xtx <- t(x) %*% x
  xtx.inv <- solve(xtx)
  mu.hat.m <- xtx.inv  %*% t(x) %*% y
  W <-  1/(1+ (data$parm$sig2e/data$parm$tau))  #solve( (xtx) + ( solve(data$Sigma.mu.m ) * data$parm$sig2e ) ) %*% xtx
  
  # If the dim of new X is same as mu_m then we do the usual update
  # But when no. of cov in x are LESS (i.e. after removing teh degenrate ones's),
  # Then the mu_mt of those covariates can be considered as zero, i.e., for degereate x, the beta coefficient will be simply zero
  A.star <- W * mu.hat.m  #+ (diag(nrow(W)) - W) %*% matrix(data$hyper.param$prior.mu$mean,nrow =nrow(xtx) )
  B.star <- data$parm$sig2e * (W * xtx.inv)
  
  mu.m <- rmvnorm(n=1,mean = A.star,sigma = B.star)
  if(!x.redes[["redesign.flag"]]){
    mu.m = mu.m
  }else{
    reduced<- data.frame(cov.name = colnames(x),mu.m=t(mu.m))
    temp.0 <- data.frame(cov.name=names(x.redes[["degen"]]),mu.m= 0)
    temp.mu <- rbind(reduced,temp.0)
    temp.mu_1 <- temp.mu[match(names(x.orig), temp.mu$cov.name),]
    mu.m <- temp.mu_1$mu.m
  }
  return(mu.m)
}

update.tau.mu <- function(data){
  a.tau <- data$hyper.param$a.tau  +  (data$parm$beta.parm$M * ncol(data$X))/2
  xtx <- data$hyper.param$prior.mu$xtx
  mu.m <- data$parm$beta.parm$mu.m
  sum.sq <- sum(sapply(1:nrow(mu.m) , function(x)  t(mu.m[x,]) %*% xtx %*% mu.m[x,] ))
  b.tau <- data$hyper.param$b.tau + sum.sq/2
  data$parm$tau <- rinvgamma(n=1,shape = a.tau,scale = (b.tau) )
  return(data)
}

# update.mu.m <- function(data){
#   v.ku.mat <- matrix(data$parm$beta.parm$v.ku,nrow=true_parm$K)
#   index.m <- which(data$parm$beta.parm$d.m!=0)
#   post.mu <- array(,c(data$parm$beta.parm$M,ncol(data$X)))
#   
#   # x1=data$X[data$grp==1,]
#   # x0=data$X[data$grp==0,]
#   eval(parse(text= paste( "x",1:true_parm$K,"=data$X[data$grp==",1:true_parm$K,",]",sep="" )  ))
#   
#   eta.hat.C <- data$parm$eta.hat[,-1]
#   for(m in index.m){
#     
#     # ind.m1 <- which(v.ku.mat[1,]==m)  # this=0 means grp 2, = 1 means grp 1
#     # ind.m2 <- which(v.ku.mat[2,]==m)
#     
#     eval(parse(text=  paste("ind.m",1:true_parm$K,"<- which(v.ku.mat[",1:true_parm$K , ",]==m)",sep="") ))
#     
#     if(length(ind.m1)!=0 ){
#       y1=matrix(eta.hat.C[data$grp==1,ind.m1],ncol=length(ind.m1))
#       y11 <- matrix(c(y1),ncol=1)
#     } 
#     if(length(ind.m2)!=0 ){
#       y2=matrix(eta.hat.C[data$grp==0,ind.m2],ncol=length(ind.m2))
#       y22 <- matrix(c(y2),ncol=1)
#     }
#     
#     if(length(ind.m1)!=0 && length(ind.m2)!=0){
#       y <- rbind(y11,y22)
#       x <- eval(parse(text= paste("rbind(",paste(paste(rep("x1",length(ind.m1)),collapse=","),
#                                                  paste(rep("x0",length(ind.m2)),collapse=","),sep=","),")",sep="")  ))
#     }
#     if(length(ind.m1)==0 && length(ind.m2)!=0){
#       y <- (y22)
#       x <- eval(parse(text= paste("rbind(",paste(paste(rep("x0",length(ind.m2)),collapse=","),sep=","),")",collapse="")  ))
#     }
#     if(length(ind.m1)!=0 && length(ind.m2)==0){
#       y <- (y11)
#       x <- eval(parse(text= paste("rbind(",paste(paste(rep("x1",length(ind.m1)),collapse=","),sep=","),")",collapse="")  ))
#     }
#     
#     post.mu[m,] <- posterior.mu(y,x,data)
#   }
#   
#   data$parm$beta.parm$mu.m <- post.mu
#   return(data)
# }
# 
# update.mu.from.prior <- function(data){
#   ss <- data$parm$tau * solve(data$hyper.param$prior.mu$xtx)
#   mea <- rep(0,ncol(data$X))
#   mu <- rmvnorm(n=1,mean=mea, sigma = ss)
#   return(mu)
# }



# update.mu.m.new0 <- function(data){
#   v.ku.mat <- matrix(data$parm$beta.parm$v.ku,nrow=true_parm$K)
#   index.m <- which(data$parm$beta.parm$d.m!=0)
#   post.mu <- array(,c(data$parm$beta.parm$M,ncol(data$X)))
#   not.in.vku <- setdiff(c(1:data$parm$beta.parm$M),index.m)
#   x1=data$X[data$grp==1,]
#   x0=data$X[data$grp==0,]
#   eta.hat.C <- data$parm$eta.hat[,-1]
#   for(m in index.m){
#     ind.m1 <- which(v.ku.mat[1,]==m)  # this=0 means grp 2, = 1 means grp 1
#     ind.m2 <- which(v.ku.mat[2,]==m)
#     if(length(ind.m1)!=0 ){
#       y1 = matrix(eta.hat.C[data$grp==1,ind.m1],ncol=length(ind.m1))
#       b01 = matrix(data$parm$beta0iu.est[data$grp==1,ind.m1],ncol=length(ind.m1))
#       y11 <- matrix(c(y1)-c(b01),ncol=1)
#     } 
#     if(length(ind.m2)!=0 ){
#       y2=matrix(eta.hat.C[data$grp==0,ind.m2],ncol=length(ind.m2))
#       b02 = matrix(data$parm$beta0iu.est[data$grp==0,ind.m2],ncol=length(ind.m2))
#       y22 <- matrix(c(y2)-c(b02),ncol=1)
#     }
#     
#     if(length(ind.m1)!=0 && length(ind.m2)!=0){
#       y <- rbind(y11,y22)
#       x <- eval(parse(text= paste("rbind(",paste(paste(rep("x1",length(ind.m1)),collapse=","),
#                                                  paste(rep("x0",length(ind.m2)),collapse=","),sep=","),")",sep="")  ))
#     }
#     if(length(ind.m1)==0 && length(ind.m2)!=0){
#       y <- (y22)
#       x <- eval(parse(text= paste("rbind(",paste(paste(rep("x0",length(ind.m2)),collapse=","),sep=","),")",collapse="")  ))
#     }
#     if(length(ind.m1)!=0 && length(ind.m2)==0){
#       y <- (y11)
#       x <- eval(parse(text= paste("rbind(",paste(paste(rep("x1",length(ind.m1)),collapse=","),sep=","),")",collapse="")  ))
#     }
#     
#     post.mu[m,] <- posterior.mu.new(y,x,data)
#   }
#   
#   for(mm in not.in.vku){
#     post.mu[mm,] <- update.mu.from.prior(data)
#   }
#   
#   data$parm$beta.parm$mu.m <- post.mu
#   return(data)
# }

update.mu.m.new0 <- function(data){
  v.ku.mat <- matrix(data$parm$beta.parm$v.ku,nrow=true_parm$K)
  index.m <- which(data$parm$beta.parm$d.m!=0)
  post.mu <- array(,c(data$parm$beta.parm$M,(ncol(data$X)+1))) # change ----
  not.in.vku <- setdiff(c(1:data$parm$beta.parm$M),index.m)
  # x1=data$X[data$grp==1,]
  # x0=data$X[data$grp==0,]
  
  x.k <- lapply(1:true_parm$K, function(x) data$X[data$grp==x,])
  
  eta.hat.C <- data$parm$eta.hat[,-1]
  for(m in index.m){

    ind.mk <- lapply(1:true_parm$K, function(x) which(v.ku.mat[x,]==m))
    
    ykk <- NULL
    for(k in 1:true_parm$K){
      if(length(ind.mk[[k]])!=0 ){
        yk=matrix(eta.hat.C[data$grp==k,ind.mk[[k]]],ncol=length(ind.mk[[k]]))
        #b0k = matrix(data$parm$beta0iu.est[data$grp==k,ind.mk[[k]]],ncol=length(ind.mk[[k]])) # change ----
        ykk[[k]] <- matrix(c(yk),ncol=1) #matrix(c(yk)-c(b0k),ncol=1)    # change ----
      }
      
    }

    #length(ykk)
    # check this  !!!!! ----
    y <- 0
    x <- array(,c(1,ncol(data$X)))
    for(ll in 1: length(ykk)){
      if(!is.null(ykk[[ll]])){
        y <- rbind(y,ykk[[ll]])
        x=rbind(x,eval(parse(text=paste("rbind(",paste(rep(paste("x.k[[",ll,"]]",sep=""),length(ind.mk[[ll]])  ),collapse=","),")")  )))
      }
    }
    y <- y[-1]
    x <- x[-1,]
    print(paste("-----------  ",m))
    post.mu[m,] <- posterior.mu.new(y,x,data)
  }
  
  for(mm in not.in.vku){
    post.mu[mm,] <- update.mu.from.prior(data)
  }
  
  data$parm$beta.parm$mu.m <- post.mu
  return(data)
}


update.intercept.hyperParms <- function(data){
  NN <- (data$n*data$C)
  beta0.mean <- mean(c(data$parm$beta0iu.est))
  a1 <- data$hyper.param$a0.int + NN/2
  nu1 <- data$hyper.param$nu0.int  + NN
  m1 <- ((data$hyper.param$m0.int * data$hyper.param$nu0.int) + (NN * beta0.mean))/ nu1
  term2 <- ((data$hyper.param$nu0.int * NN)/nu1) *( beta0.mean- data$hyper.param$m0.int )^2
  b1 <- data$hyper.param$b0.int + ( 0.5*((NN-1)*var(c(data$parm$beta0iu.est)) ))  + (0.5*term2) 
  data$parm$tau0 <- rinvgamma(n=1, shape = a1,scale = b1)
  tau02 <- data$parm$tau0/nu1
  data$parm$mu0 <- rnorm(n=1,mean = m1,sd =sqrt(tau02) )
  return(data)
}

#######
update.intercept <- function(data){
  
  xbeta <- beta0iu <- Riu <- array(,c(data$n,data$C))
  xx <- data$X
  grp <- data$grp #ifelse(data$grp==1,1,2)
  v.ku.mat <- matrix(data$parm$beta.parm$v.ku,nrow=true_parm$K)
  for(i in 1:data$n){
    for(uu in 1:data$C){
      xbeta[i,uu] <-   xx[i,]  %*% data$parm$beta.parm$mu.m[ (v.ku.mat[grp[i],uu]),  ]
    }
  }
  
  Riu <- data$parm$eta.hat[,-1] - xbeta
  
  vv.inv <- ((1/data$parm$tau0) + (1/data$parm$sig2e))
  for(i in 1:data$n){
    for(uu in 1:data$C){
      riu <- ((Riu[i,uu]/data$parm$sig2e)  + (data$parm$mu0/data$parm$tau0)) /vv.inv
      vv <- 1/vv.inv
      beta0iu[i,uu] <- rnorm(n = 1,mean = riu,sd = sqrt(vv))
    }
  }
  data$parm$beta0iu.est <- beta0iu
  return(data)
}



####### Update allocation vector c.j
# code in clustAlloc_mcmc.R

# This function updates the parameters in data$parm


mcmc.update.case2_randintercept.w0 <- function(data,hyper.param){
  # update sigm.e^2
  data <- update.sigma.e0(data,hyper.param)
  data <- update.tau.mu(data)   # AFTER UPDATING TAU2
  # update tau^2 and mu.m , m=1,...,M
  data <- update.mu.m.new0(data)  # updates tau and mu.m
  #data <- update.intercept.hyperParms(data) # change ----
  #data <- update.intercept(data)
  
  # update beta_ku, pi, v_ku, da.status, d.m
  data <- update.beta.FMM.new0(data)
  da.info <- da.status.new(data)
  data$parm$h.j <-da.info # update h.j
  
  # c.j's and the related parameters remain the same in this case
  # update c.j's 
  # data <- update.cj.new(data)
  # data$parm$c.j, data$parm$m.u doesnt change
  # AFTER UPDATING C'J's update QSTAR AS WELL (IMP!!!)
  data.check.c.j.full <- update.cj.fullLik.corrected(data)
  data <- data.check.c.j.full$data
  data$parm$check$full.log.lik <- data.check.c.j.full$full.log.lik
  
  data <- update.qstar.new(data)
  #check simplex
  #plot(apply(data$parm$qstar.hat.mat, 1, function(x) sum(data$parm$m.u*x)),main = "check simplex qstar.hat.mat",ylab="")
  
  data <- update.eta.iu0(data)
  #check.simplex.mcmc(data)  
  data <- update.qstar.new(data)
  
  data <- impute.zeroes.mcmc(data)
  
  return(data)
}




