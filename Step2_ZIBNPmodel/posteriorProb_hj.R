# changed mcmc0.R
# posterior probabilities are stored in the following component
# The prob h.u =1 can be calculated as follows


# for mcmc iteration b=1,
# This gives a vector of length C x 1
# func.h.u.post.prob <- function(data1){
#   post.pi <- data1[["parm"]][["check"]][["pi.post.vku"]]
#   h.u.post.prob <- vector()
#   
#   for(u in 1:data1$C){
#     h.u.post.prob[u] <- sum(post.pi[[1]][u,] * post.pi[[2]][u,])
#     # post.pi[[1]][u,1] * post.pi[[2]][u,1] # pu11 * pu12
#   }
#   return(h.u.post.prob)
# }


func.h.u.post.prob <- function(data1){
  post.pi <- data1[["parm"]][["check"]][["pi.post.vku"]]
  h.u.post.prob <- vector()
  
  for(u in 1:data1$C){
    prod.prob <- eval(parse(text=  paste(paste("post.pi[[",1:true_parm$K,"]][u,]",sep=""),collapse="*")  ))
    # post.pi[[1]][u,] * post.pi[[2]][u,]
    h.u.post.prob[u] <- sum(prod.prob)
    # post.pi[[1]][u,1] * post.pi[[2]][u,1] # pu11 * pu12
  }
  return(h.u.post.prob)
}



func.post.prob.c.j <- function(data1){
  #log.lik.c.j <- data.check.c.j.full[["full.log.lik"]]
  log.lik.c.j <- data1[["parm"]][["check"]][["full.log.lik"]]
  
  log.lik.c.j1 <- t(apply(log.lik.c.j,1, function(x)  exp(x-max(x))))
  pi.sum <- apply(log.lik.c.j1,1,sum)
  log.lik.c.j1 <- log.lik.c.j1/pi.sum
  # replace zeros by "small"
  small2 <- 1e-5
  log.lik.c.j1[log.lik.c.j1 < small2] <- small2
  # again normalize
  pi.sum1 <- apply(log.lik.c.j1,1,sum)
  log.lik.c.j <- log.lik.c.j1/pi.sum1   
  log.lik.c.j0 <- rbind(NA,log.lik.c.j)
  return(log.lik.c.j0)
}



func.posteriorProb.hj <- function(data1){
  h.u.post.prob <- func.h.u.post.prob(data1)
  log.lik.c.j0 <- func.post.prob.c.j(data1)
  h.j.prob <- t(h.u.post.prob %*% t(log.lik.c.j0))
  h.j.prob[1] <- 0 
  return(h.j.prob)
}



# table(ifelse(round(1-h.u.post.prob,4)>0.5,1,0),true_parm$da.status[-1]-1)
# 
# round(h.u.post.prob,4)
# round(h.j.prob,4)
# 
# 
# (true_parm$h.j[1:10])
# (round(h.j.prob,4)[1:10])

# pred.h.j <- ifelse(h.j.prob>=th0,1,0) 
