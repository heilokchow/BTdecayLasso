BTLagrangian <- function(Lagrangian, ability, theta, penalty.Qua) {
  n <- nrow(ability)-1
  v <- penalty.Qua
  
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      Lagrangian[i, j] <- Lagrangian[i, j] + v * (theta[i, j] - ability[i, 1] + ability[j, 1])
    }
  }
  return(Lagrangian)
}

