#' Bradley-Terry Model with Exponential Decayed weighted likelihood and Adaptive Lasso with a given penalty rate
#' 
#' This function provides a method to computed the estimated abilities and lambda given an intuitive fixed Lasso penalty rate.
#' Since in Lasso method, the selection of lambda varies a lot with respect to different datasets. We can keep the consistency of
#' amount of Lasso penalty induced in different datasets from different period by setting a fixed Lasso penalty rate "p".
#' 
#'
#' @param dataframe Generated using \code{\link{BTdataframe}} given raw data.
#' @param ability A column vector of teams ability, the last row is the home parameter. It can be generated using \code{\link{BTdataframe}} given raw data.
#' The row number is consistent with the team's index shown in dataframe. It can be generated using \code{\link{BTdataframe}} given raw data.
#' @param penalty The amount of Lasso penalty induced (1-s/max(s)) where is the sum of Lasso penalty part.
#' @param decay.rate The exponential decay rate. Usually ranging from (0, 0.01), A larger decay rate weights more
#' importance to most recent matches and the estimated parameters reflect more on recent behaviour.
#' @param fixed A teams index whose ability will be fixed as 0. The worstTeam's index
#' can be generated using \code{\link{BTdataframe}} given raw data.
#' @param thersh Threshold for convergency
#' @param max Maximum weight for w_{ij} (weight used for Adaptive Lasso)
#' @param iter Number of iterations used in L-BFGS-B algorithm.
#' @details
#' The estimated ability given fixed penalty \eqn{p = 1- s/\max(s)} where s is the sum of Lasso penalty part.
#' When p = 0, this model is reduced to a standard Bradley-Terry Model.
#' When p = 1, all ability scores are shrinking to 0.
#' 
#' The parameter "penalty", which is p, should be ranging from 0.01 to 0.99 due to the iteration's convergent error.
#' The objective likelihood function to be optimized is,
#' \deqn{\sum_{k=1}^{n}\sum_{i<j}\exp(-\alpha t_{k})\cdot(y_{ij}(\tau h_{ij}^{t_{k}}+\mu_{i}-\mu_{j})-\log(1+\exp(\tau h_{ij}^{t_{k}}+\mu_{i}-\mu_{j})))}
#' With the Lasso constraint,
#' \deqn{\sum_{i<j}w_{ij}\left|\mu_{i}-\mu_{j}\right|\leq s}
#' where n is the number of matches, \eqn{\alpha} is the exponential decay rate, \eqn{\tau} is the home parameter and 
#' \eqn{y_{ij}} takes 0 if i is defeated by j, 1 otherwise. \eqn{\mu_{i}} is the team i's ability score and penalty is 1-s/max(s).
#'
#' summary() function can be applied to view the outputs.
#' @return The list with class "BTF" contains estimated abilities and other parameters.
#' \item{ability}{Estimated ability scores}
#' \item{df}{Degree of freedom (number of distinct \eqn{\mu})}
#' \item{penalty}{Amount of Lasso Penalty}
#' \item{decay.rate}{Exponential decay rate}
#' \item{lambda}{Corresponding Lasso lambda given penalty rate}
#' @seealso \code{\link{BTdataframe}} for dataframe initialization
#' @references 
#' Masarotto, G. and Varin, C.(2012) The Ranking Lasso and its Application to Sport Tournaments. 
#' *The Annals of Applied Statistics* **6** 1949--1970.
#' 
#' Zou, H. (2006) The adaptive lasso and its oracle properties. 
#' *J.Amer.Statist.Assoc* **101** 1418--1429.
#' @examples
#' ##Initializing Dataframe
#' x <- BTdataframe(NFL2010)
#' 
#'\dontrun{
#' ##BTdecayLasso run with exponential decay rate 0.005 and Lasso penaty 0.5
#' y <- BTdecayLassoF(x$dataframe, x$ability, 0.5, decay.rate = 0.005, 
#'                    fixed = x$worstTeam)
#' summary(y)
#' }
#' 
#' @export

BTdecayLassoF <- function(dataframe, ability, penalty, decay.rate = 0, fixed = 1, thersh = 1e-5, max = 100, iter = 100) {
  
  
  if (penalty > 1 || penalty < 0) {
    stop("Please provide a penalty ranging from 0 to 1")
  }
  
  df <- dataframe
  n <- nrow(ability) - 1
  df[, 5] <- df[, 5] - df[1, 5]
  p <- 1 - penalty
  
  if(!(fixed %in% seq(1, n, 1))){
    stop("The fixed team's index must be an integer index of one of all teams")
  }
  
  if (p > 0.99) {
    BT <- BTdecay(df, ability, decay.rate = decay.rate, fixed = fixed, iter = iter)
    ability <- BT$ability
    s <- list(ability = round(ability, -log10(thersh)), df = n, penalty = 0, decay.rate = decay.rate, lambda = 0)
  } else {
    if (p < 0.01) {
      s <- list(ability = ability, df = 1, penalty = 1, decay.rate = decay.rate, lambda = Inf)
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
    }
  }
  class(s) <- "BTF"
  s
}


