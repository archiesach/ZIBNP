fn.zero.gen <- function(data,lambda){
  lam.ku <- matrix(lambda * c(rep(1,ncol(data$X1)+2)),ncol=ncol(data$X1)+2)
  dd <- data.frame(cbind(data$X,log.Li =log(data$lib.size)))
  form <- as.formula(paste("~",paste(names(dd),collapse = "+")))
  mm <- model.matrix(form,dd)
  
  X.L <- lam.ku %*% t(mm)
  e.X.L <- exp(X.L)
  rij <- e.X.L/(e.X.L+1) 
  range(rij)
  
  # generate a binary matrix using rij
  ppp <- data$p-1
  delta.zero.mat <- array(,c(nrow(data$X),ncol=(ppp)))
  
  delta.zero.mat <- t(sapply(1:data$n,function(x) rbinom(n = ppp,size = 1,prob = (1-rij[x]) )  ))
  #mean(delta.zero.mat==0)
  data$Zmat.wout0 <- temp.mat <- data$Zmat
  temp.mat <- temp.mat[,-1]
  temp.mat[delta.zero.mat==0] = 0
  data$Zmat.w0 <- cbind(data$Zmat.wout0[,1],temp.mat )
  
  perct.zero <- mean(data$Zmat.w0==0)
  range(apply(data$Zmat.w0,2, function(x) mean(x==0)))
  print(paste("~~~~~~~~~~~~~~ % of zero in the data = ", round(perct.zero*100,2)))
  ifelse(mean(delta.zero.mat==0) == mean(data$Zmat.w0[,-1]==0), print("There are NO sampling zeroes"), print("There may be sampling zeroes")) 
  data$perct.zero <- perct.zero
  data$delta.zero.mat.true <- delta.zero.mat
  data$lib.size.wout0  <- data$lib.size
  data$lib.size <- apply(data$Zmat.w0,1,sum)
  data$lib.size.w0 <- data$lib.size
  data$Zmat <- data$Zmat.w0
  return(data)
}
