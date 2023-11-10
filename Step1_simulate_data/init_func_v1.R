# Initialization 

rel.abund.of.Z <- function(data){
  data$init$qhat.mat <- t(apply(data$Zmat,1, function(x) x/sum(x)))
  return(data)
}



my_round <- function(x)  {
  x1 <- round(x, 2)
  if(round(x, 2) == 0) {
    n1 <- stringr::str_locate(x, "[^-0.]")[1] -  str_locate(x, fixed("."))[1]
    print(n1)
    x1 <- round(x, n1)
  }
  return(x1)
}


san.qstar.hat <- function(qq,cj){
  cc <- max(cj)
  q.star.hat.san <- array(,c(nrow(qq),cc))
  for(i in 1:cc){
    q.star.hat.san[,i] <- apply(qq[,which(cj==1)],1,mean)
  }
  return(q.star.hat.san)
}




init.clust.new0 <-function(data){
  #data <- psuedo.count(data,dec=10) # removed the psuedo count concept

    # cluster all columns except the refernce column
    qq <- data$init$qhat.mat[,-1]   
    #set.seed(7240)
    km <- kmeans(t(qq),centers = data$C)
    qstar.hat.mat <- t(km$centers) 
    data$init$qstar.hat.mat <- cbind(data$init$qhat.mat[,1], qstar.hat.mat)
    data$init$c.j <- c(1,km$cluster+1)
    data$init$G <- max(data$init$c.j)  
    data$init$m.u <- tabulate(data$init$c.j,nbins =data$init$G )
  

  return(data)
}



init.eta.star <- function(data){
  # note that first column of data$init$qstar.hat.mat is reference cluster 
  
  qq <- data$init$qstar.hat.mat[,-1]
  qq1 <- data$init$qstar.hat.mat[,1]
  eta.hat  <-  log(apply(qq,2,function(x)  x/qq1))
  data$init$eta.hat <- cbind(0, eta.hat)
  return(data)
}


lm.beta.eta <- function(dta){
  coef <- lm(eta~ .- 1,data=as.data.frame(dta ))$coef
  coef[is.na(coef)] <- 0
  return(coef)
}

init.get.beta <- function(data){

  beta.init <- array(,c(true_parm$K,ncol(data$X),(data$C)))
  eta.h <- data$init$eta.hat[,-1]
  grp <- unique(data$grp)
  for( i in 1:data$C){
    eta=eta.h[,i]
    for(k in 1:true_parm$K){
      dta <-cbind(eta=eta[data$grp==grp[k]],data$X[data$grp==grp[k],])
      beta.init[k,,i] <-  lm.beta.eta(dta)
    }
  }
  data$init$beta.init <- beta.init
  return(data)
}

beta.kmeans.init <- function(mat,init.M){
  beta.parm <- NULL
  #set.seed(17240)
  # CHANGE IT ----
  if(data.type=="sim.data"){
    km <- kmeans((mat),centers = init.M)  # same as data generation M (just for initial stage ) ---- 
  } 
  if(data.type=="real.data"){
    km <- kmeans((mat),centers = round(log(true_parm$K*data$C)))
  }
  
  beta.parm$mu.m <- (km$centers)  
  beta.parm$v.ku <- km$cluster
  beta.parm$M <- max(beta.parm$v.ku)  # excluding the reference cluster
  beta.parm$d.m <- tabulate(beta.parm$v.ku,nbins = beta.parm$M )
  #prob.pi <- (true_parm$hyperparam$a0/beta.parm$M) +beta.parm$d.m
  beta.parm$pi <- beta.parm$d.m/sum(beta.parm$d.m)  #rdirichlet(n=1,alpha = prob.pi)
  return(beta.parm)
}



da.status.new <- function(data){
  parm=data$parm
  h.c.status <- parm$beta.parm$da.status
  c.j <- parm$c.j    # that's because it also marks the reference cluster
  h.j <- vector()
  h.j[1] <- 1
  for(c in 2:data$G){
    h.j[c.j==c] =h.c.status[c]
  }
  return(h.j)
}

da.status.true.parm <- function(true_parm){
  h.c.status <- true_parm[["da.status"]]
  c.j <- true_parm$clust$c.v    # that's because it also marks the reference cluster
  h.j <- vector()
  h.j[1] <- 1
  for(c in 2:true_parm$clust$G){
    h.j[c.j==c] =h.c.status[c]
  }
  true_parm$h.j <- h.j  # this extra 1 is for reference cluster which is always nDA
  return(true_parm)
}

which.cat <- function(x,cat){
  cat.sum <- vector()
  for(i in 1:cat){
    cat.sum[i] <- sum(x==i)
  }
  return(cat.sum)
}


S.summ <-function(mu.m_mubar){
  mu.m_mubar1 <- NULL
  for(i in 1:ncol(mu.m_mubar)){
    x <- mu.m_mubar[i,]
    x <- matrix(x,ncol=1)
    mu.m_mubar1[[i]] <- x %*% t(x)
  }
  mu.m_mubar1
}



init.da.new <- function(data){
  beta.ku <- apply(data$init$beta.init,2,rbind) #  b11,b21,...,b1C,b2C
  
  data$init$beta.parm <- beta.kmeans.init(beta.ku,init.M=data$init$init.M)
  data$init$beta.parm$beta.mu <- data$init$beta.parm$mu.m[data$init$beta.parm$v.ku,]
  
  # extra 1 is for reference cluster
  data$init$beta.parm$da.status <- c(1,apply(matrix(data$init$beta.parm$v.ku,2),2, function(x)  length(unique(x))))
  data.temp <- data
  data.temp$parm<- data.temp$init
  temp <-da.status.new(data.temp)
  rm(data.temp)
  data$init$h.j <-temp  #t(apply(b.ku,1,da.status))
  #data$init$beta.parm$da.status <- temp$da.stat
  return(data)
}

init.sig2 <- function(data){
  eta.hat <- data$init$eta.hat[,-1]
  y.1 <- as.vector(eta.hat)
  x.1 <-  eval(parse(text=paste("rbind(",paste(rep("data$X",ncol(eta.hat)),collapse=","),")",sep="")))
  lm.init <- lm(y.1~x.1)
  data$init$sig2e <- (summary(lm.init)$sigma)^2
  return(data)
}


fn.Dist <- function(c.v,G)
{
  tmp.mat <- array(0,c(G,G))
  for (jj in 1:G)
  {indx.jj <- which(c.v==jj)
  tmp.mat[indx.jj,indx.jj] <- 1
  }
  
  Dist.mt = tmp.mat
  
  Dist.mt
}


clust.alloc.func <- function(cv,meandist.mat){
  l <- nrow(meandist.mat)
  tmp.mat <- array(0,c(l,l))
  g<-max(cv)
  for (jj in 1:g)
  {
    indx.jj <- which(cv==jj)
    tmp.mat[indx.jj,indx.jj] <- 1
  }
  
  f.norm <- sum((tmp.mat-meandist.mat)^2)
  return(f.norm)
}





init.tau2 <- function(data){
  mu.m <- data$init$beta.parm$mu.m
  x <- data$X
  xtx <- t(x) %*% x
  ssq <- sum(sapply(1:nrow(mu.m) , function(x)  t(mu.m[x,]) %*% xtx %*% mu.m[x,] ))
  tau2 <- ssq / (data$init$beta.parm$M * ncol(x))
  data$init$tau <- tau2
  return(data)
}



init.beta.intercept <- function(data){
  et.hat <- data$init$eta.hat[,-1]
  h.u <- data$init$beta.parm$da.status[-1]
  xx <- data.frame(id =as.factor(c(1:data$n)),grp=data$grp,data$X)
  colnames(xx) <- c("id", "grp","X1","X2","X3" ,"X4")
  bet0iu.init <- array(,c(data$n,data$C))
  
  for(uu in 1: length(h.u)){
    d <- cbind(eta=et.hat[,uu],xx )
    # glm.fit <- lmer(eta ~ X1+X2+X3 +X4+ (1|id), data = d)
    if(h.u[uu]==1){
      glm.fit1 <- lme(eta ~ X1+X2+X3 +X4, random=  ~ 1|id, data=d)
      inf.model.intercept <- coef(glm.fit1)[,"(Intercept)"]
    }
    if(h.u[uu]==2){
      glm.fit11 <- lme(eta ~ X1+X2+X3 +X4 , random=  ~ 1|id, data=d[d$grp==1,])
      glm.fit12 <- lme(eta ~ X1+X2+X3 +X4 , random=  ~ 1|id, data=d[d$grp==0,])
      inf.model.intercept <- c(coef(glm.fit11)[,"(Intercept)"],coef(glm.fit12)[,"(Intercept)"])
    }
    bet0iu.init[,uu] <- inf.model.intercept
  
  }
  
  data$init$beta0iu.est <- bet0iu.init
  return(data)
  
}



initialize.new.case1 <- function(data, true_parm,init.M){
  if(data.type=="sim.data"){
    data$G <- true_parm$clust$G # true_parm$updated.clust$C
    data$C <- true_parm$clust$C
    data$init$init.M <- init.M
  } else if (data.type=="real.data"){
    data$C <-  round(log(ncol(data$Zmat)))
  }
  
  
  data <- rel.abund.of.Z(data)   # This gives data$init$qhat.mat
  print(apply(data$init$qhat.mat,1,sum))

  
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

  
  data <- init.eta.star(data)    # This gives data$init$eta.hat
 
  data <- init.sig2(data)        # this initializes sigma_e^2
  
  
  data <- init.get.beta(data)    # this initializes beta : data$init$beta.init
  
  data <- init.da.new(data)          # This provides the parameters of beta FMM : data$init$beta.parm using K-MEANS clustering
  
  data <- init.tau2(data) 
  #data$init$beta.parm$beta.mu <- data$init$beta.parm$mu.m[data$init$beta.parm$v.ku,]

  # Initialize beta0iu
  data <- init.beta.intercept(data)
  
  return(data)
}

