penaltyAmount <- function(ability, weight){
  n <- nrow(ability) - 1
  s <- 0
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      s <- s + abs(ability[i, 1] - ability[j, 1])*weight[i, j]
    }
  }
  return(s)
}