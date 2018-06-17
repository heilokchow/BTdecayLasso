#' Bradley-Terry Model with Exponential Decayed weighted likelihood and weighted Lasso with AIC or BIC criteria
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
#' It can be generated using function BTdataframe.
#' @param weight Weight for Lasso penalty on different abilities
#' @param Lambda A sequence of Lambda
#' @param criteria "AIC" or "BIC"
#' @param type "HYBRID" or "LASSO"
#' @param decay.rate The exponential decay rate. Usually ranging from (0, 0.1), A larger decay rate weights more
#' importance to most recent matches and the estimated parameters reflect more on recent behaviour.
#' @param fixed A teams index whose ability will be fixed as 0 (usually the team loss most which can be
#' generated using function BTdataframe).
#' @param thersh Thershold for convergency
#' @param max Maximum weight for w_{ij} (weight used for Adaptive Lasso)
#' @param iter Number of iterations used in L-BFGS-B algorithm.
#' @details Model selection through AIC or BIC method
#' @return
#' \item{Score}{Lowest AIC or BIC score}
#' \item{Optimal.degree}{The degree of freedom where lowest AIC or BIC score is achieved}
#' \item{Optimal.ability}{The ability where lowest AIC or BIC score is achieved}
#' \item{ability}{Matrix contains all abilities computed in this algorithm}
#' \item{Optimal.lambda}{The lambda where lowest score is attained}
#' \item{Optimal.penalty}{The penalty (s/\eqn{\max(s)}) where lowest score is attained}
#' \item{type}{Type of model selection method}
#' @examples 
#' ##Initializing Dataframe
#' x <- BTdataframe(NFL2010)
#' 
#' ##Model selection through AIC
#' w <- BTLasso.weight(NFL$df, NFL$ability)
#' z <- BTdecayLassoC(NFL$df, NFL$ability, w, criteria = "AIC", type = "LASSO")
#' @export

BTdecayLassoC <- function(dataframe, ability, weight = NULL, lambda = NULL, criteria = "AIC", type = "HYBRID", decay.rate = 0, 
                          fixed = 1, thersh = 1e-5, iter = 100, max = 100) {
  Lp <- BTdecayLasso(dataframe, ability, lambda = lambda, weight = weight, path = TRUE, decay.rate = decay.rate,
                     fixed = fixed, thersh = thersh, max = max, iter = iter)
  
  n1 <- nrow(dataframe)
  if (criteria == "AIC") {
    mul <- 2
  } else if (criteria == "BIC") {
    mul <- log(n1)
  } else {
    stop("criteria should either be AIC or BIC")
  }
  y <- Lp$df.path * mul
  
  x <- 2 * Lp$HYBRID.likelihood + y
  ind <- which.min(x)
  ind <- ind[length(ind)]
  dg <- Lp$df.path[ind]
  
  if (type == "HYBRID") {
    output <- list(Score = min(x), Optimal.degree = dg, Optimal.ability = Lp$ability.path[ind], 
                   Optimal.lambda = Lp$Lambda.path[ind], Optimal.penalty = Lp$penalty.path[ind], type = type)
  } else if (type == "LASSO") {
    m0 <- Lp$Lambda.path[ind]
    m1 <- Lp$Lambda.path[ind + 1]
    k <- m0 - m1
    while (k > thersh * 1e2) {
      k <- k * 0.5
      m1 <- m0 - k
      BT <- BTdecayLasso.step2(dataframe, ability, m1, weight, decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
      if (dg == BT$df) {
        m0 <- m1
      } 
    }
    l <- BTLikelihood(dataframe, BT$ability, decay.rate = decay.rate)
    output <- list(Score = 2 * l + dg * mul, Optimal.degree = dg, Optimal.ability = BT$ability,
                   Optimal.lambda = m1, Optimal.penalty = BT$penalty, type = type)
  } else {
    stop("Please provide a selection type HYBRID or LASSO")
  }
  output
}

