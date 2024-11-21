#testing how to write function when passing in a vector
#for the endec_adjusted_mortlity() function in biting_processes.R

cov <- c(0.0, 0.1, 0.1)
mu <- c(0.3, 0.3, 0.3)
time <- 10
out <- function(mu, cov){
  if (cov > 0){
    mu_endec <- (mu + 0.1)*cov 
    return(mu_endec)
  }
  
  else {return(mu)}
    
}
out(mu, cov, time)


#DO IT LIKE THIS
out <- function(mu, cov, time){
  mu_me <- mu
  kill_me <- (mu+0.1)*cov
  if(time == 10){
    mu_endec <- ifelse(cov > 0, kill_me, mu_me)
    return(mu_endec)
  }
  #mu_endec <- ifelse(time == 10 & cov > 0, kill_me, mu_me)
  
  
}
