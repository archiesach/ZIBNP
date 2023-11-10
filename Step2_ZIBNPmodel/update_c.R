update.likl.c.j <- function(data,u,j.ind){
  data.temp <- data
  #print(u)
  data.temp$parm$c.j[j.ind] <- u
  data.temp$parm$m.u <- tabulate(data.temp$parm$c.j,nbins = data$G)
  xi.mat.temp <- summary.log.f.eta(data.temp)
  qstar.temp <- array(,c(data$n,data$G))
  qstar.temp[,1] <- 1/(apply(exp(data.temp$parm$eta.hat),1, function(x)  sum(x * data.temp$parm$m.u)))
  
  for(uu in 2:data$G){
    qstar.temp[,uu] <- exp(data.temp$parm$eta.hat[,uu]) * qstar.temp[,1]
  }
  #print(j.ind)
  #print(sum(apply(qstar.temp, 1 , function(x)  sum(x * data.temp$parm$m.u))))
  #print(sum(apply(xi.mat.temp,1,sum) == data.temp$lib.size))
  data.temp$parm$qstar.hat.mat <- qstar.temp
  
  LL <- NULL
  LL$likl <- sum(log(data.temp$parm$qstar.hat.mat) * xi.mat.temp)
  #print(u)
  #print(data.temp$parm$m.u)
  LL$m.u <- data.temp$parm$m.u[u]
  # changing the data parameters Qstar, m.u, and c.j
  return(LL)
}



update.cj.fullLik.corrected <- function(data){
  
  if(data$c.j.perc<1){
    c.update <- sample(x=c(1:(data$p)),size =(data$c.j.perc*data$p),replace = FALSE)
    #c.update <- setdiff(c.update,1)
  }
  if(data$c.j.perc==1){
    c.update <- c(1:(data$p))
    c.update <- setdiff(c.update,1)
  }
  
  log.lik.mat <- array(,c(length(c.update),data$C))
  full.log.lik <- array(,c(length(c.update),data$C))
  mu.prior <- array(,c(length(c.update),data$C))
  change =0
  c.j.temp <-vector()
  c.j.temp[1] <- 1
  ct <- 0
  #prior.mu <- data$parm$m.u[-1] + (hyper.param$a.clus/data$C)
  
  for(ll in c.update){
    #print(ll)
    ct<- ct+1
    j.ind <- ll
    u.ind <- data$parm$c.j[j.ind]
    
    for( u in 2:data$G){
      LL <- update.likl.c.j(data,u,j.ind)
      full.log.lik[ct,(u-1)] <- LL$likl
      mu.prior[ct,(u-1)] <- LL$m.u
    }
    
    log.lik.u <- full.log.lik[ct,]  # log likl for a particular j
    prior.mu <- (mu.prior[ct,])
    
    log.P.u <- log.lik.u + log(prior.mu)
    log.lik.mat[ct,] <- log.P.u
    lp <- exp(log.P.u -max(log.P.u))
    post.pi.u <- lp/sum(lp)
    #round(post.pi.u,2)
    
    c.j.post <- which(rmultinom(n=1,size=1,post.pi.u)==1)
    c.j.temp[j.ind] <- c.j.post + 1
    
    #update.current.state
    data$parm$c.j[j.ind] <- c.j.temp[j.ind]
    data$parm$m.u <- tabulate(data$parm$c.j,nbins = (data$G))
    data <- update.qstar.new(data)
    if(ll %% 10 ==0){
      print(paste("   ****** updating c.j ********   ",ll,sep=""))
    }
    
  }
  
  print(table(c.j.temp))
  # update m.u and qstar
  #data$parm$c.j <- c.j.temp
  print(table(data$parm$c.j))
  #update m.u and qstar
  #data$parm$m.u <- tabulate(data$parm$c.j,nbins = (data$G))
  #data <- update.qstar.new(data)
  
  data.out <- NULL
  data.out$data <- data
  data.out$full.log.lik <- full.log.lik

  return(data.out)
  
  #return(data)
}
