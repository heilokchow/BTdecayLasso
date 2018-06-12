#' Bradley-Terry Model with Exponential Decayed weighted likelihood and weighted Lasso with a given Lasso lambda and Lasso weight
#' 
#' @description 
#' Normal Adaptive Lasso's weight can be determined using funcion BTLasso.weight.
#' If you want to use equal weight, when there are n teams, you can define a n*n upper triangular matrix with diagonal is 0 and other entries of upper-right part are 1
#' @param dataframe Matrix with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' It can be generated using function BTdataframe.
#' @param ability A column vector of teams ability, the last row is the home parameter.
#' The row number is consistent with the team's index shown in dataframe.
#' It can be generated using function BTdataframe.
#' @param lambda The amount of Lasso penalty induced.
#' @param weight Weight for Lasso penalty on different abilities
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
#' @return
#' \item{ability}{Estimated ability scores}
#' \item{df}{Degree of freedom (number of distinct \eqn{\mu})}
#' \item{penalty}{s/max(s)}
#' @examples
#' ##Initializing Dataframe
#' x <- BTdataframe(NFL2010)
#' 
#' ##Defining Adaptive weight
#' w1 <- BTLasso.weight(x$df, x$ability, 0.1, decay.rate = 0.005)
#' 
#' ##BTdecayLasso.step2 run with exponential decay rate 0.005 and lambda 0.1
#' y1 <- BTdecayLasso.step2(x$df, x$ability, 0.1, w1, decay.rate = 0.005, fixed = x$worstTeam)
#' 
#' ##Defining equal weight
#' n <- nrow(x$ability) - 1
#' w2 <- matrix(1, nrow = n, ncol = n)
#' w2[lower.tri(w2, diag = TRUE)] <- 0
#' 
#' ##BTdecayLasso.step2 run with exponential decay rate 0.005 and lambda 0.1
#' y2 <- BTdecayLasso.step2(x$df, x$ability, 0.1, w2, decay.rate = 0.005, fixed = x$worstTeam)
#' @export
BTdecayLasso.step2 <- function(dataframe, ability, lambda, weight, decay.rate = 0, fixed = 1, thersh = 1e-5, iter = 100) {
  u <- decay.rate
  n <- nrow(ability) - 1
  theta <- matrix(0, nrow = n, ncol = n) 
  Lagrangian <- matrix(0, nrow = n, ncol = n) 
  ability[, 1] <- 0
  
  stop <- 0
  j <- 1
  v <- 10
  while (stop==0) {
    ability <- BTdecayLasso.step1(dataframe, ability, weight, Lagrangian, theta, v, lambda, 
                                  decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
    theta <- BTtheta(ability, weight, Lagrangian, v, lambda)
    Lagrangian0 <- BTLagrangian(Lagrangian, ability, theta, v)
    k <- sum(abs(Lagrangian0 - Lagrangian))
    if (k < thersh) {
      stop <- 1
    } else {
      Lagrangian <- Lagrangian0
      v <- max(Lagrangian^2)
    }
    s <- penaltyAmount(ability, weight)
  }
  cat(s,'\n')
  
  ability0 <- ability
  ability0[, 1] <- 0
  BT <- BTdecay(dataframe, ability0, decay.rate = decay.rate, fixed = fixed, iter = iter)
  ability0 <- BT$ability
  s0 <- penaltyAmount(ability0, weight)
  p <- s/s0
  
  degree <- round(ability[1:n, 1], -log10(thersh))
  degree <- length(unique(degree))
  
  output <- list(ability = ability, df = degree, penalty = p)
  output
}
