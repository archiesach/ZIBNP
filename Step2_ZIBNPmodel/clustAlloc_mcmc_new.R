

# Think about c.j's whether to include the indexing for the reference taxa, it is getting messy here. 
# Anyway the dimension of the below likelihood is C and c.j has C+1 indexes, so decide how to take care of this


# This function is required only at the beginning of an mcmc iteration 
# Likhood.new <- function(data){
#   log.l <- vector()  #data$parm$log.l
#   # this indexing doesnt change with j
#   ind.j <- lapply(1:(data$G),function(x) which(data$parm$c.j==x))   # indexing ----
#   #for(j in c.update){
#   for(uu in 1:(data$G)){
#     #u_c.j <- data$parm$c.j[j]
#     # since ind.j[[1]] = 1 (reference cluster), we start with uu=2
#     if(length(ind.j[[uu]])>1){
#       xi.u <- apply(data$Zmat[,ind.j[[uu]] ],1,sum)
#     }
#     if(length(ind.j[[uu]])==1){
#       xi.u <- data$Zmat[,ind.j[[uu]]]
#     }
#     if(length(ind.j[[uu]])==0){
#       xi.u <- 0
#     }
#     
#     log.qu <- log(data$parm$qstar.hat.mat[,uu])
#     log.l[(uu)] <- sum(xi.u * log.qu) 
#   }
#   
#   #}
#   data$parm$log.l <- log.l
#   return(data)
# }




# update.cj <- function(data){
#   if(data$c.j.perc<1){
#     c.update <- sample(x=c(1:(data$p+1)),size =(data$c.j.perc*data$p),replace = FALSE)
#     #c.update <- setdiff(c.update,1)
#   }
#   if(data$c.j.perc==1){
#     c.update <- c(1:(data$p+1))
#     #c.update <- setdiff(c.update,1)
#   }
# 
#   data$parm$c.j.update <- c.update
# 
#   data <- Likhood.new(data) # Likelihood updates based on the previous c.j allocations
#   change =0
#   for(ll in c.update){
#     #print(ll)
#     j.ind <- ll
#     u.ind <- data$parm$c.j[j.ind]
#     #print(j.ind)
#     Z.j <- data$Zmat[,j.ind]
# 
#     log.star.j.ind <- log(data$parm$qstar.hat.mat[,u.ind])  # dim data$parm$qstar.hat.mat is n x (C+1)
#     log.u <- log(data$parm$qstar.hat.mat)-log.star.j.ind
#     Q.u <- apply(log.u,2,function(x) sum(x*Z.j))
# 
#     log.lik.u <- data$parm$log.l + Q.u  # dim data$parm$log.l is  C x 1
# 
# 
#     prior.mu <- data$parm$m.u+ (hyper.param$a.clus/data$C)
#     log.P.u <- log.lik.u + log(prior.mu)
# 
#     lp <- exp(log.P.u -max(log.P.u))
#     post.pi.u <- lp/sum(lp)
# 
#     c.j.post <- which(rmultinom(n=1,size=1,post.pi.u)==1)
# 
#     data$parm$c.j[j.ind] <- c.j.post #+1  # updated allocation for j.ind, that's because ref cluster =1 and c.j=2 actually means cluster 1 out of C clusters
#     # if c.j.post == u.ind then no change in allocation
#     change <- change + (sum((c.j.post+1) != u.ind))
#   }
# 
#   print(paste("***** No. of changes in c.j allocation vector : ", change, " out of ", (data$c.j.perc*(data$p+1))))
#   # update m.u
#   data$parm$m.u <- tabulate(data$parm$c.j,nbins = (data$G))
# 
#   #check.simplex.mcmc(data)
#   data <- update.qstar.new(data)
#   return(data)
# }



update.cj.new <- function(data){
  if(data$c.j.perc<1){
    c.update <- sample(x=c(1:(data$p)),size =(data$c.j.perc*data$p),replace = FALSE)
    #c.update <- setdiff(c.update,1)
  }
  if(data$c.j.perc==1){
    c.update <- c(1:(data$p))
    #c.update <- setdiff(c.update,1)
  }
  
  data$parm$c.j.update <- c.update
  
  
  data <- Likhood.new(data) # Likelihood updates based on the previous c.j allocations
  change =0
  for(ll in c.update){
    #print(ll)
    j.ind <- ll
    u.ind <- data$parm$c.j[j.ind]
    #print(j.ind)
    Z.j <- data$Zmat[,j.ind]
    
    log.star.j.ind <- log(data$parm$qstar.hat.mat[,u.ind])  # dim data$parm$qstar.hat.mat is n x (C+1)
    log.u <- log(data$parm$qstar.hat.mat)-log.star.j.ind
    Q.u <- apply(log.u,2,function(x) sum(x*Z.j)) 
    
    log.lik.u <- data$parm$log.l + Q.u  # dim data$parm$log.l is  C x 1
    
    
    prior.mu <- data$parm$m.u+ (hyper.param$a.clus/data$C)
    log.P.u <- log.lik.u + log(prior.mu)
    
    lp <- exp(log.P.u -max(log.P.u))
    post.pi.u <- lp/sum(lp)
    
    c.j.post <- which(rmultinom(n=1,size=1,post.pi.u)==1)
    
    data$parm$c.j[j.ind] <- c.j.post #+1  # updated allocation for j.ind, that's because ref cluster =1 and c.j=2 actually means cluster 1 out of C clusters
    # if c.j.post == u.ind then no change in allocation
    data <- Likhood.new(data)
    data$parm$m.u <- tabulate(data$parm$c.j,nbins = (data$G))
    change <- change + (sum((c.j.post+1) != u.ind))
  }
  print(paste("****** updating c.j ********   "))
  print(paste("***** No. of changes in c.j allocation vector : ", change, " out of ", (data$c.j.perc*(data$p))))
  # update m.u and qstar
  #data$parm$m.u <- tabulate(data$parm$c.j,nbins = (data$G))
  data <- update.qstar.new(data)
  
  # applt the swap function
  swap.fun <- swap.index.refclust(c.j= data$parm$c.j, qstar=data$parm$qstar.hat.mat)
  data$parm$qstar.hat.mat <- swap.fun$qstar.swap
  data$parm$c.j <- swap.fun$c.j.swap
  data$parm$m.u <- tabulate(data$parm$c.j,nbins =data$G )
  
  
  return(data)
}



update.cj.reduced <- function(data){
  if(data$c.j.perc<1){
    c.update <- sample(x=c(1:(data$p)),size =(data$c.j.perc*data$p),replace = FALSE)
    #c.update <- setdiff(c.update,1)
  }
  if(data$c.j.perc==1){
    c.update <- c(1:(data$p))
    c.update <- setdiff(c.update,1)
  }
  
  # c.update <- which((data$parm$c.j-1)==4)
  #true.c.j <- data$parm$c.j[c.update]-1
  log.lik.mat <- Q.u.mat <- array(,c(length(c.update),data$C))
  
  likld.mat.init <- data$parm$log.l
  c.j.temp <-vector()
  ct <- 0
  qstar.C <- data$parm$qstar.hat.mat[,-1]
  prior.mu <- data$parm$m.u[-1]+ (hyper.param$a.clus/data$C)
  #c.update <- which(data$parm$c.j==ll)
  #data.temp <- data
  for(ll in c.update){
    #print(ll)
    ct<- ct+1
    j.ind <- ll
    #u.ind <- data$parm$c.j[j.ind]
    #print(j.ind)
    Z.j <- data$Zmat[,j.ind]
    
    Q.u <-apply(log(qstar.C),2,function(x) sum(x*Z.j)) 
    Q.u.mat[ct,] <- Q.u
    log.lik.u <- Q.u #data$parm$log.l + Q.u  # dim data$parm$log.l is  C x 1
    
    log.P.u <- log.lik.u + log(prior.mu)
    log.lik.mat[ct,] <- log.P.u
    lp <- exp(log.P.u -max(log.P.u))
    post.pi.u <- lp/sum(lp)
    
    c.j.post <- which(rmultinom(n=1,size=1,post.pi.u)==1)
    c.j.temp[j.ind] <- c.j.post
    
    #change <- change + (sum((c.j.post+1) != u.ind))
  }
  # c.j.temp[c.update]
  print(paste("****** updating c.j ********   "))
  #print(paste("***** No. of changes in c.j allocation vector : ", change, " out of ", (data$c.j.perc*(data$p+1))))
  #update m.u and qstar
  data$parm$c.j <- c.j.temp
  m.u.local <- tabulate(data$parm$c.j,nbins = (data$C))
  data$parm$m.u <- c(1, m.u.local)
  #data <- update.qstar.new(data)
  
  # applt the swap function
  # swap.fun <- swap.index.refclust(c.j= data$parm$c.j, qstar=data$parm$qstar.hat.mat)
  # data$parm$qstar.hat.mat <- swap.fun$qstar.swap
  # data$parm$c.j <- swap.fun$c.j.swap
  # data$parm$m.u <- tabulate(data$parm$c.j,nbins =data$G )
  
  
  return(data)
}

