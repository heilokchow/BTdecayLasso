#' Bradley-Terry Model with Exponential Decayed weighted likelihood and Adaptive Lasso with a given penalty rate
#'
#' @param dataframe Matrix with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' It can be generated using function BTdataframe.
#' @param ability A column vector of teams ability, the last row is the home parameter.
#' The row number is consistent with the team's index shown in dataframe.
#' #' It can be generated using function BTdataframe.
#' @param penalty The amount of Lasso penalty induced (1-s/max(s)) where is the sum of Lasso penalty part.
#' @param decay.rate The exponential decay rate. Usually ranging from (0, 0.1), A larger decay rate weights more
#' importance to most recent matches and the estimated parameters reflect more on recent behaviour.
#' @param fixed A teams index whose ability will be fixed as 0 (usually the team loss most which can be
#' generated using function BTdataframe).
#' @param thersh Thershold for convergency
#' @param max Maximum weight for w_{ij} (weight used for Adaptive Lasso)
#' @param iter Number of iterations used in L-BFGS-B algorithm.
#' @details
#' The objective likelihood function to be optimized is,
#' \deqn{\sum_{k=1}^{n}\sum_{i<j}\exp(-\alpha t_{k})\cdot(y_{ij}(\tau h_{ij}^{t_{k}}+\mu_{i}-\mu_{j})-\log(1+\exp(\tau h_{ij}^{t_{k}}+\mu_{i}-\mu_{j})))}
#' With the Lasso constraint,
#' \deqn{\sum{i<j}w_{ij}\left|\mu_{i}-\mu_{j}\right|\leq s}
#' where n is the number of matches, \eqn{\alpha} is the exponential decay rate, \eqn{\tau} is the home parameter and 
#' \eqn{y_{ij}} takes 0 if i is defeated by j, 1 otherwise. \eqn{\mu_{i}} is the team i's ability score and penalty is 1-s/max(s).
#' This likelihood function is optimized using L-BFGS-B method with package \bold{optimr}.
#' @return The list contains estimated abilities and penalty.
#' \item{ability}{Estimated ability scores}
#' \item{df}{Degree of freedom (number of distinct \eqn{\mu})}
#' \item{penalty}{Amount of Lasso Penalty}
#' @examples
#' ##Initializing Dataframe
#' x <- BTdataframe(NFL2010)
#' 
#' ##BTdecayLasso run with exponential decay rate 0.005 and Lasso penaty 0.5
#' y <- BTdecayLassoF(x$df, x$ability, 0.5, decay.rate = 0.005, fixed = x$worstTeam)
#' summary(y)
#' @export

BTdecayLassoF <- function(dataframe, ability, penalty, decay.rate = 0, fixed = 1, thersh = 1e-5, max = 100, iter = 100) {
  
  
  df <- dataframe
  n <- nrow(ability) - 1
  df[, 5] <- df[, 5] - df[1, 5]
  p <- 1 - penalty
  
  if (p == 1) {
    BT <- BTdecay(df, ability, decay.rate = decay.rate, fixed = fixed, iter = iter)
    ability <- BT$ability
    s <- rbind(ability, matrix(c(0), ncol = 1))
    return(s)
  } else {
    weight <- BTLasso.weight(df, ability, decay.rate = decay.rate, fixed = fixed, thersh = thersh, max = max, iter = iter)
    BT1 <- BTdecayLasso.step2(df, ability, 0, weight, decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
    ability <- BT1$ability
    k0 <- penaltyAmount(ability, weight)
    BT1 <- BTdecayLasso.step2(df, ability, 0.1, weight, decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
    ability <- BT1$ability
    k1 <- penaltyAmount(ability, weight)
    
    a <- -log(k1/k0)/0.1
    x0 <- 0.1
    m0 <- -log(k1/k0)
    x1 <- -log(p)/a
    mo <- -log(p)
    
    stop <- 0
    while (stop == 0) {
      BT1 <- BTdecayLasso.step2(df, ability, x1, weight, decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
      ability <- BT1$ability
      k2 <- penaltyAmount(ability, weight)
      m1 <- -log(k2/k0)
      
      x2 <- (x0 - x1)/(m0 - m1)*(mo - m0) + x0
      if(abs(k2/k0 - p) < thersh/10){
        stop <- 1
      } else if (((k2/k0) < thersh/10) && (x2 > max(x1, x0))) {
        x0 <- max(x1, x0)
        x1 <- x2
      } else if (x2 < 0) {
        x1 <- 0.9 * x1
      } else {
        x <- x1
        x1 <- x2
        x0 <- x
        m0 <- m1
      }
    }
    
    s <- list(ability = round(ability, -log10(thersh)), df = BT1$df, penalty = penalty, decay.rate = decay.rate, lambda = x1)
    class(s) <- "BTF"
    s
  }
}


