#' @export
BTLasso.weight <- function(dataframe, ability, decay.rate = 0, fixed = 1, thersh = 1e-5, max = 100, iter = 100) {
  BT <- BTdecay(dataframe, ability, decay.rate = decay.rate, fixed = fixed, iter = iter)
  ability <- BT$ability
  n <- nrow(ability) - 1
  weight <- matrix(0, nrow = n, ncol = n)
  
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      if(abs(ability[i, 1] - ability[j, 1]) < 1/max){
        weight[i, j] <- -1
      } else{
        weight[i, j] <- 1/abs(ability[i, 1] - ability[j, 1])
      }
    }
  }
  
  k <- max(weight)
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      if(weight[i, j] == -1){
        weight[i, j] <- k
      }
    }
  }
  weight
}