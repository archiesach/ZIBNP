


impute.zeroes <- function(data){
  
  cj <-  data$init$c.j
  
  Ji <- data$init$Ji
  
  #data$Zmat[x,Ji[[x]]]
  # compute qi.tilde
  qij.Ji <- lapply(1:data$n, function(x) data$init$qstar.hat.mat[x,cj[Ji[[x]]]])
  qi.tilde <- sapply(1:data$n, function(x) sum( qij.Ji[[x]] ) )  
  length(qi.tilde)
  plot(qi.tilde)
  # generate Si.tilde
  Si.tilde <- sapply(1:data$n, function(x) rnbinom(n=1,size = data$lib.size.w0[x] ,prob = (1-qi.tilde[x])) )
  
  plot(Si.tilde,data$lib.size,xlim=c(min(Si.tilde,data$lib.size), max(Si.tilde,data$lib.size)), ylim=c(min(Si.tilde,data$lib.size), max(Si.tilde,data$lib.size)) )
  # generate Zij.tilde
  
  Zij.tilde <- lapply(1:data$n, function(x)  rmultinom(n=1,size = Si.tilde[x] ,prob = qij.Ji[[x]] ) )
  

  # impute Zij
  mean(data$Zmat==0)
  impute.z <- data$Zmat
  print(mean(impute.z==0))
  sum(impute.z==0)
    for(x in 1:data$n){
      impute.z[ x,Ji[[x]]] = c(Zij.tilde[[x]])
    }
  
  mean(data$Zmat==0)
  print(mean(impute.z==0))
  no.0 <- sum(impute.z==0)
  # if there is a zero after imputing for Technical zeroes, it could be a sampling zero
  # therefore 
  # if(sum(impute.z==0)>0){
  #   which(impute.z==0,arr.ind = T)
  #   impute.z[impute.z==0] =1
  # }
  
  data$Zmat <- impute.z
  data$lib.size <- apply(data$Zmat,1,sum)
  data$init$Ji <- Ji
  data$init$Si.tilde <- Si.tilde
  data$init$no.0 <- no.0
  data$init$qi.tilde <- qi.tilde
  return(data)
}


impute.zeroes.mcmc <- function(data){
  cj <-  data$parm$c.j
  
  # Ji <- lapply(1:data$n, function(x) which(data$Zmat[x,]==0))  
  Ji <- data$parm$Ji
  
  # compute qi.tilde
  
  qi.tilde <- sapply(1:data$n, function(x) sum( data$parm$qstar.hat.mat[x,cj[Ji[[x]]]] ) )  
  length(qi.tilde)
  plot(qi.tilde)
  # generate Si.tilde
  Si.tilde <- sapply(1:data$n, function(x) rnbinom(n=1,size = data$lib.size.w0[x] ,prob = (1-qi.tilde[x])) )
  
  #plot(Si.tilde,data$lib.size,xlim=c(min(Si.tilde,data$lib.size), max(Si.tilde,data$lib.size)), ylim=c(min(Si.tilde,data$lib.size), max(Si.tilde,data$lib.size)) )
  # generate Zij.tilde
  
  Zij.tilde <- lapply(1:data$n, function(x)  rmultinom(n=1,size = Si.tilde[x] ,prob = data$parm$qstar.hat.mat[x,cj[Ji[[x]]]]) )
  
  
  # impute Zij
  #mean(data$Zmat==0)
  impute.z <- data$Zmat.w0 #cbind(rep(1,data$n),data$Zmat.w0)
  print(mean(impute.z==0))
  sum(impute.z==0)
  for(x in 1:data$n){
    impute.z[ x,Ji[[x]]] = c(Zij.tilde[[x]])
  }
  
  #mean(data$Zmat==0)
  print(mean(impute.z==0))
  no.0 <- sum(impute.z==0)
  # if there is a zero after imputing for Technical zeroes, it could be a sampling zero
  # therefore 
  # if(sum(impute.z==0)>0){
  #   which(impute.z==0,arr.ind = T)
  #   impute.z[impute.z==0] =1
  # }
  
  data$Zmat <- impute.z
  data$lib.size <- apply(data$Zmat,1,sum)
  data$parm$Si.tilde <- Si.tilde
  data$parm$no.0 <- no.0
  return(data)
}

init.till.qstar <- function(data, true_parm,init.M){
  if(data.type=="sim.data"){
    data$G <- true_parm$clust$G + 1 # true_parm$updated.clust$C
    data$C <- true_parm$clust$C + 1
    data$init$init.M <- init.M
  } else if (data.type=="real.data"){
    data$C <- true_parm$C  #round(log(ncol(data$Zmat)))
    data$G <- data$C + 1
    data$init$init.M <- init.M
  }
  
  
  data <- rel.abund.of.Z(data)   # This gives data$init$qhat.mat
  plot(apply(data$init$qhat.mat,1,sum))
  
  
  # now find the reference cluster by the logic of minimum variance 
  # to tag that as the reference column and consider the cluster it belongs to -to be reference cluster
  
    if(ref.clust){
      
      # z.mat.var <- apply(data$Zmat,2,var)
      # ref.z.ind <- which(z.mat.var==min(z.mat.var))
      ind1 <- data$init$qhat.mat[1,]
      which(ind1==0)
      mat2 <- t(apply(data$init$qhat.mat[-1,], 1, function(x) x/ind1))
      dev <- (mat2-1)^2
      qty2 <- apply(dev, 2, mean )
      min.ind.ratio <- which(qty2==min(qty2,na.rm = TRUE))
      ref.z.ind <- min.ind.ratio
      print(paste("  *********   index of reference cluster  : ",ref.z.ind))
      if(length(ref.z.ind)==1){
        z1 <- data$Zmat[,1]
        data$Zmat[,1]  <- data$Zmat[,ref.z.ind]
        data$Zmat[,ref.z.ind]<- z1
      }else{
        print(paste("Z mat have more than one columns with minimum variance. No. of such taxa col.= ",length(ref.z.ind)))
        ref.z.ind <- ref.z.ind[1]
        z1 <- data$Zmat[,1]
        data$Zmat[,1]  <- data$Zmat[,ref.z.ind]
        data$Zmat[,ref.z.ind]<- z1
        print("picking one of them")
      }
      data$init$ref.z.ind <- ref.z.ind
    }else{
      data$init$ref.z.ind <- 1
    }
  

  
  
  # INSTEAD OF init.clust.new we are using true values  if use.init.true.cj is TRUE
  if(use.init.true.cj == TRUE){
    data$init$c.j <- true_parm$clust$c.v
    data$init$G <- true_parm$clust$G
    data$init$m.u <- true_parm$clust$C.m.vec
    data$init$qstar.hat.mat <- data$qstar  # or can take avergae of q_ij's as per c.j'c and using matrix qhat.mat
  } else{
    data <- init.clust.new0(data)       # This gives data$init$ ,"qstar.hat.mat","c.j","Ghat","m.u"
  }
  return(data)
}


initialize.new.case.w0 <- function(data, true_parm,init.M){
  data$init$Ji <- lapply(1:data$n, function(x) which(data$Zmat[x,]==0))  
  orig.zmat <- data$Zmat # this will be used at the end for imputing Zmat
  
  # FIRST clustering with zeros included 
  # is useful for estimating lambda's (n x C)
  data <- init.till.qstar(data, true_parm,init.M)
  
  # estimate lambda_ku (FIGURE THIS OUT) # new ----
  # data <- fn.init.lambda_ku(data) # first element corresponding to ref taxa of lambda will be NA 
  # Just to check lambda estimates
  # lamm <- c(data$init$lambda)
  # plot(lamm[!is.na(lamm)])
  # abline(h=lambda0,col="red")
  
  # SECOND clustering by adding 1's 
  # is useful for estimating all other parameters
  # replace all zero's with 1 first to get qstars 
  data$Zmat <- data$Zmat + 1
  data$lib.size <- apply(data$Zmat,1,sum)

  # after adding 1, redo the clustering to get q.hat.mat, c.j's,  qstar
  data <- init.till.qstar(data, true_parm,init.M)
  
  data <- init.eta.star(data)    # This gives data$init$eta.hat
  
  data <- init.sig2(data)        # this initializes sigma_e^2
  
  
  data <- init.get.beta(data)    # this initializes beta : data$init$beta.init
  
  data <- init.da.new(data)          # This provides the parameters of beta FMM : data$init$beta.parm using K-MEANS clustering
  
  data <- init.tau2(data) 
  #data$init$beta.parm$beta.mu <- data$init$beta.parm$mu.m[data$init$beta.parm$v.ku,]
  
  # Initialize beta0iu
  # data <- init.beta.intercept(data) # change ----
  
  # Impute the zeroes
  data$Zmat <- orig.zmat
  data$lib.size <- apply(data$Zmat,1,sum)
  data <- impute.zeroes(data)
  
  
  return(data)
}