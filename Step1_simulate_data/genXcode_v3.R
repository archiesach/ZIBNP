# chnaged the sampling of library size as compared to the genXcode.R script
# changed the scaling of Xmat, not it is not scaled

X.mat1 <- pheno.all 
X.mat1 <- dplyr::select(X.mat1,antibiotics_last_6_months,bmi,brush_teeth_daily,floss_each_day,periodontal_disease,sex,african_american,asian,white  )   # X.mat1[, -which(colnames(X.mat1)=="percent_body_fat") ] # highly corrlated with BMI


grp <- pheno.all %>%  dplyr::select(adult)
data.type <- "sim.data" #"real.data"


#want.T <- 4# 40

#if(want.T<=5){
X.mat <-  dplyr::select(X.mat1, antibiotics_last_6_months,bmi,brush_teeth_daily,sex)
#X.mat <- cbind(1,X.mat)
#}

# CHANGE : using same covariate matrix for both adult and child ----
if(do.subjects.same==TRUE){
  n1 <- sum(grp==1)
  n0 <- sum(grp==0)
  #set.seed(1983)
  ind1 <- sample(1:n1,size=n/2)
  ind0 <- sample(1:n0,size=n/2)
  X.mat<- rbind(X.mat[grp==1,][ind1,], X.mat[grp==1,][ind1,])
  grp <- c(grp[grp==1,][ind1], grp[grp==0,][ind0])
  #L_i <- round((1/li.factor) * (sample(x = lib.size$L_i,size = length(grp),replace = T)),0)  #rpois(n=n, lambda = 10000)*rpois(n=n, lambda = 100) #rep(5000,length(c(ind1,ind0))) #lib.size$L_i[c(ind1,ind0)]
  if(libsize=="rsamp"){
    L_i <- round((1/li.factor) * (sample(x = lib.size$L_i,size = length(grp),replace = T)),0)  #rpois(n=n, lambda = 10000)*rpois(n=n, lambda = 100) #rep(5000,length(c(ind1,ind0))) #lib.size$L_i[c(ind1,ind0)]  
  }
  if(libsize=="poisson1"){
    L_i <- rpois(n=n, lambda = pois.lambda)
  }
  if(libsize=="poisson2"){
    L_i <- rpois(n=n, lambda = 10000)*rpois(n=n, lambda = 100)
  }
  check<- apply(X.mat,2, function(x) length(unique(x)))
  print(check)
  if(length(which(check==1))>=1){
    X.mat <- X.mat[,-which(check==1)]
  }
}

if(do.subjects.same==FALSE){
  n1 <- sum(grp==1)
  n0 <- sum(grp==0)
  #set.seed(1983)
  ind1 <- sample(1:n1,size=n/2)
  ind0 <- sample(1:n0,size=n/2)
  X.mat<- rbind(X.mat[grp==1,][ind1,], X.mat[grp==0,][ind0,])
  grp <- c(grp[grp==1,][ind1], grp[grp==0,][ind0])
  #L_i <- round((1/li.factor) * (sample(x = lib.size$L_i,size = length(grp),replace = T)),0) #rpois(n=n, lambda = 10000)*rpois(n=n, lambda = 100) #rep(5000,length(c(ind1,ind0))) #lib.size$L_i[c(ind1,ind0)]
  if(libsize=="rsamp"){
    L_i <- round((1/li.factor) * (sample(x = lib.size$L_i,size = length(grp),replace = T)),0)  #rpois(n=n, lambda = 10000)*rpois(n=n, lambda = 100) #rep(5000,length(c(ind1,ind0))) #lib.size$L_i[c(ind1,ind0)]  
  }
  if(libsize=="poisson"){
    L_i <- rpois(n=n, lambda = pois.lambda)
  }
  check<- apply(X.mat,2, function(x) length(unique(x)))
  print(check)
  if(length(which(check==1))>=1){
    X.mat <- X.mat[,-which(check==1)]
  } 
  
  # if(ncol(X.mat)<want.T){
  #   gen <- want.T-ncol(X.mat)
  #   x.plus <- rmvnorm(n=nrow(X.mat)/2,mean=c(rep(0,gen)), sigma = diag(gen) )
  #   X.mat <- cbind(X.mat,rbind(x.plus,x.plus))
  # }
}


normalize.minmax <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

dim(X.mat)
data<- NULL
data$X1 <- X.mat
data$X1$antibiotics_last_6_months <- as.factor(data$X1$antibiotics_last_6_months)
data$X1$brush_teeth_daily <- as.factor(data$X1$brush_teeth_daily)
data$X1$sex <- as.factor(data$X1$sex)
#apply(X.mat,2,mean)
#apply(X.mat,2,sd)

X.mat <- data.frame(X.mat) #sapply(1:ncol(X.mat),function(x)  scale(X.mat[,x]))
X.mat$bmi <- X.mat$bmi/10
# 1) X.mat$bmi/10
# 2) min-max scaling = normalize.minmax(X.mat$bmi)
data$grp <- as.vector(unlist(grp))
data$X <-  as.matrix(X.mat) #as.matrix(cbind(int0=1,X.mat))
 
data$lib.size <- L_i
