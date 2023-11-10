



gen.beta.norm <- function(C,tau,X){
  m <- rep(0,ncol(X))
  v.mat <- tau^2 * solve(t(X) %*% X)
  bet <- rmvnorm(n = C,mean = m,sigma = v.mat)
  return(bet)
}

psi <- function(psi0,K,fold.change){
  psi <- vector()
  psi[1] <- psi0 
  psi[2] <- psi0 + log(fold.change) 
  return(psi)
}


da.status.gen <- function(pi0,C){
  dtemp <- rbinom(n = C,size = 1,prob = pi0)
  da.status <- c(1,dtemp+1)
  return(da.status)
}

beta.ku <- function(psi.vec,bet,K,da.status,C){
  bet.ku <- NULL
  da.status <- da.status[-1]
  bet.ku[[1]] <- cbind(psi.vec[1], bet)
  bet.ku[[2]] <- array(,c(nrow(bet.ku[[1]]),ncol(bet.ku[[1]])))
  for(u in 1:C){
    if(da.status[u]==1){
      bet.ku[[2]][u,] <- c(psi.vec[1], bet[u,])
    } else if(da.status[u]==2){
      bet.ku[[2]][u,] <- c(psi.vec[2], bet[u,])
    }
  }
  
  return(bet.ku)
}


sigma.forZ <- function(cj,mu.iu,p){
  aij <- mu.iu[,cj]
  sigm.e <- var(as.vector(log(aij))) *((1/R2.aij) -1)
  return(sigm.e)
}


gen.Zij <- function(lam0,rel.A.ij ){
  #Li <- rpois(n = nrow(rel.A.ij),lambda = lam0 )
  Li <- lam0
  Zij <- t(sapply(1:nrow(rel.A.ij), function(x)  rmultinom(n = 1,size = Li[x],prob = rel.A.ij[x,]) ))
  
  return(Zij)
  
}
  
  

counts.ij <- function(X,beta.ku,grp,p,cj){
  grp <- grp+1
  X1 <- cbind(1,X)
  mu.iu <- array(,c(n,C))
  for(i in 1:n){
    for(u in 1:C){
      b.ki.u <-  beta.ku[[grp[i]]][u,]
      mu.iu[i,u] <- (X1[i,]) %*% (b.ki.u)
    }
  }
  
  p <- p-1

  sig = sigma.forZ(cj,mu.iu,p)
  
  A.ij <- array(,c(n,p))
  for(i in 1:n){
    for(j in 1:(p)){
      A.ij[i,j] <- rlnorm(n=1, meanlog =mu.iu[i,cj[j]] , sdlog = sig)
    }
  }
  
  rel.A.ij <- t(sapply(1:nrow(A.ij), function(x) A.ij[x,]/sum(A.ij[x,])  ))
  lam0 <- round(apply(A.ij,1,sum),0)
  Zij <- gen.Zij(lam0,rel.A.ij )
  
  
  return(list(mu.iu=mu.iu,Zij=Zij,A.ij=A.ij))
  
}



