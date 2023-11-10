

xi.summary <- function(ind.j,Zmat,G){
  xi.sumi <- lapply(1:G, function(x) sum(Zmat[ind.j[[x]] ])  )
  return(unlist(xi.sumi))
}

# Run this function for each u, 1,...,C and it can be vectorized for i's
summary.log.f.eta <- function(data){
  # calculate xi
  ind.j <- lapply(1:(data$G),function(x) which(data$parm$c.j==x))   # indexing ----
  xi <- t(sapply(1:data$n ,function(x) xi.summary(ind.j,data$Zmat[x,],G=data$G) ))
  return(xi)
}



log.f.eta <- function(x.eta, u ,xi,eta,Li,m.u,sig2.e, X ,betai,b0){
  #print(u)
  xi.iu <- xi[u]
  xi.iw <- xi[-u]
  eta.iw <- eta[-u]
  x_beta <- matrix(X,nrow=1) %*% matrix(betai,ncol=1)
  sum1 <- sum(eta.iw * xi.iw)
  sum2 <- sum(m.u[-u]*exp(eta.iw))
  log.lik <- (x.eta * xi.iu) + sum1 - Li * (log( m.u[u]*exp(x.eta) + sum2))
  log.prior <- -(1/(2*sig2.e)) * (c(x.eta) - b0  - c(x_beta ))^2
  log.post <- log.lik +log.prior
  return(log.post)
}





d_log.f.eta <- function(x.eta,u,xi,eta,Li,m.u,sig2.e, X ,betai,b0){
  xi.iu <- xi[u]
  m.u.u <- m.u[u]
  m.u.w <- m.u[-u]
  eta.iw <- eta[-u]
  
  x_beta <- matrix(X,nrow=1) %*% matrix(betai,ncol=1)
  sum2 <- sum(m.u.w*exp(eta.iw))
  q.iu.star <- exp(x.eta)/ (m.u.u *exp(x.eta) + (sum2))
  # print("x.eta")
  # print(x.eta)
  # print("x_beta")
  # print(x_beta)
  d1.log.post <- xi.iu - (Li*m.u.u * q.iu.star) - ((1/sig2.e) * (c(x.eta)-b0 - c(x_beta )) )
  #print(u)
  return(d1.log.post)
}


update.eta.iu0 <- function(data){
  v.ku.mat <- matrix(data$parm$beta.parm$v.ku,nrow=true_parm$K)
  #label.i <- ifelse(data$grp==1,1,2)
  label.i <- data$grp
  xi.mat <- summary.log.f.eta(data) # this is a n x C matrix
  
  #all.equal(apply(xi.mat,1,sum),data$lib.size)
  eta.samp <- array(,c(data$n,ncol(v.ku.mat)))
  eta.log.concave <- array(,c(data$n,ncol(v.ku.mat)))
  
  for(i in 1:data$n){
    #beta0iu <- data$parm$beta0iu.est[i,] # change ----
    
    xi=xi.mat[i,]
    eta=data$parm$eta.hat[i,] #- c(0,beta0iu)  # This is because eta_i1 = 0
    Li=data$lib.size[i] # change check ---- this is because 5001 is with refernce zi1
    m.u=data$parm$m.u
    X=data$X[i,]
    sig2.e=data$parm$sig2e
    
    if(i%%10==0){
      print(paste("****** updating eta ********   n = ",i))
    }
    
    
    for(u in 1:(ncol(v.ku.mat)) ){ # data$C c(1,2,4)
      inde.mu <- v.ku.mat[label.i[i],u]
      betai <- data$parm$beta.parm$mu.m[inde.mu,]
      
      u1<- u + 1
      b0 <- 0 #beta0iu[u] # change ----
      #print(paste("**** for u = ", u))
      mysample <- ars(1, f=log.f.eta, fprima=d_log.f.eta,u=u1 ,lb=FALSE,ub=FALSE , xi=xi, eta=eta, Li=Li,
                      m.u=m.u,sig2.e=sig2.e ,X=c(1,X) , betai=betai,b0=b0) # change c(1,x) ----
      
      #samp <- ars_iferror(u1, xi=xi, eta=eta, Li=Li, m.u=m.u,sig2.e=sig2.e ,X=X , betai=betai)
      tt_test<- paste("tests_",tt,".txt",sep="")
      capture.output(ars(1, f=log.f.eta, fprima=d_log.f.eta,u=u1 ,lb=FALSE,ub=FALSE , xi=xi, eta=eta, Li=Li,
                         m.u=m.u,sig2.e=sig2.e ,X=c(1,X) , betai=betai,b0=b0), file = tt_test, append = FALSE) # change c(1,x) ----
      
      read.error <- read.csv(tt_test)
      # use this when the function is concave down increasing
      #if( !is.na(read.error[1,]) && (read.error[3,]=="ifault= 3 " | read.error[1,]=="ifault= 3 ") ){
      if( !is.na(read.error[1,]) && (read.error[3,]%in% c("ifault= 3 ","ifault= 4 ") | read.error[1,]%in% c("ifault= 3 ","ifault= 4 "))){
        mysample <- ars(1, f=log.f.eta, fprima=d_log.f.eta,u=u1 , ub=TRUE,lb=TRUE, xlb=-100,xub=100,xi=xi, eta=eta, Li=Li,
                        m.u=m.u,sig2.e=sig2.e ,X=c(1,X) , betai=betai,b0=b0) # change c(1,x) ----
        eta.log.concave[i,u] <- 1
      }
      
      eta.samp[i,u] <- mysample
      #print(mysample)
    }
  }
  
  #data$eta.log.concave <- eta.log.concave
  data$parm$eta.hat <- cbind(0,eta.samp)
  data$parm$check$xi.mat <- xi.mat
  return(data)
}



update.qstar.new <- function(data){
  # after updating eta, update qstar
  qstar <- array(,c(data$n,(data$G)))
  hm <- data$parm$m.u
  eta.hat <- data$parm$eta.hat
  mu_qstar <- apply(apply(exp(eta.hat) ,1, function(x) hm*x),2,sum)
  q.ref = exp(eta.hat[,1]) / (mu_qstar)
  
  qstar[,1] <- q.ref
  for(o in 2:(data$G)){
    qstar[,o] <- exp(eta.hat[,o]) *  q.ref
  }
  
  #check.simplex.mcmc(data)
  simplex.check  <- apply(apply(qstar, 1, function(x)  (x*(data$parm$m.u))),2,sum)
  #print(plot(simplex.check))
  
  data$parm$qstar.hat.mat <- qstar
  return(data)
}
