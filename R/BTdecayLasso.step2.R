BTdecayLasso.step2 <- function(dataframe, ability, lambda, weight, decay.rate = 0, fixed = 1, thersh = 1e-5, iter = 100) {
  u <- decay.rate
  n <- nrow(ability) - 1
  theta <- matrix(0, nrow = n, ncol = n) 
  Lagrangian <- matrix(0, nrow = n, ncol = n) 
  ability[,1] <- 0
  
  stop <- 0
  j <- 1
  v <- 10
  while(stop==0){
    ability <- BTdecayLasso.step1(dataframe, ability, weight, Lagrangian, theta, v, lambda, 
                                  decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
    theta <- BTtheta(ability, weight, Lagrangian, v, lambda)
    Lagrangian0 <- BTLagrangian(Lagrangian, ability, theta, v)
    k <- sum(abs(Lagrangian0 - Lagrangian))
    if(k < thersh/10){
      stop <- 1
    } else{
      Lagrangian <- Lagrangian0
      v <- max(Lagrangian^2)
    }
    s <- BTLikelihood.all(dataframe, ability, theta, v, weight, Lagrangian, lambda, decay.rate = decay.rate)
  }
  cat(s,'\n')
  return(ability)
}
