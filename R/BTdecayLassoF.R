#' @export
BTdecayLassoF <- function(dataframe, ability, penalty, decay.rate = 0, fixed = 1, thersh = 1e-5, max = 100, iter = 100) {
  df <- dataframe
  n <- nrow(ability) - 1
  df[, 5] <- df[, 5] - df[1, 5]
  p <- 1 - penalty
  
  if(p == 1){
    BT <- BTdecay(df, ability, decay.rate = decay.rate, fixed = fixed, iter = iter)
    ability <- BT$ability
    s <- rbind(ability, matrix(c(0), ncol = 1))
    return(s)
  } else{
    weight <- BTLasso.weight(df, ability, decay.rate = decay.rate, fixed = fixed, thersh = thersh, max = max, iter = iter)
    ability <- BTdecayLasso.step2(df, ability, 0, weight, decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
    k0 <- penaltyAmount(ability, weight)
    ability <- BTdecayLasso.step2(df, ability, 0.1, weight, decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
    k1 <- penaltyAmount(ability, weight)
    
    a <- -log(k1/k0)/0.1
    x0 <- 0.1
    m0 <- -log(k1/k0)
    x1 <- -log(p)/a
    mo <- -log(p)
    
    stop <- 0
    while(stop == 0){
      ability <- BTdecayLasso.step2(df, ability, x1, weight, decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
      k2 <- penaltyAmount(ability, weight)
      m1 <- -log(k2/k0)
      
      x2 <- (x0 - x1)/(m0 - m1)*(mo - m0) + x0
      if(abs(k2/k0 - p) < thersh/10){
        stop <- 1
      } else if(((k2/k0) < thersh/10) && (x2 > max(x1, x0))){
        x0 <- max(x1, x0)
        x1 <- x2
      } else if(x2 < 0){
        x1 <- 0.9 * x1
      } else{
        x <- x1
        x1 <- x2
        x0 <- x
        m0 <- m1
      }
    }
    
    s <- list(ability = round(ability, -log10(thersh)), penalty = penalty)
    s
  }
}


