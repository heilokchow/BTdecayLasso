BTdecayLasso.step1 <- function(dataframe, ability, weight, Lagrangian, theta, penalty.Qua, lambda, 
                          decay.rate = 0, fixed = 1, thersh = 1e-5, iter = 100) {
  ##stop <- 0
  ##s1 <- 1000
  ##while(stop==0){
  ##  ability <- BTdecay.Qua(dataframe, ability, theta, penalty.Qua, Lagrangian, decay.rate = decay.rate,
  ##                         fixed = fixed, iter = iter)
  ##  theta <- BTtheta(ability, weight, Lagrangian, penalty.Qua, lambda)
  ##  s <- BTLikelihood.all(dataframe, ability, theta, penalty.Qua, weight, Lagrangian, lambda, decay.rate = decay.rate)
  ##  if(abs(s-s1) < thersh){
  ##    stop <- 1
  ##  } else{
  ##    s1 <- s
  ##  }
  ##}
  ability <- BTdecay.Qua(dataframe, ability, theta, penalty.Qua, Lagrangian, decay.rate = decay.rate,
                                                fixed = fixed, iter = iter)
  return(ability)
}