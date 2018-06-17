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
#' @param criteria AIC or BIC
#' @param decay.rate The exponential decay rate. Usually ranging from (0, 0.1), A larger decay rate weights more
#' importance to most recent matches and the estimated parameters reflect more on recent behaviour.
#' @param fixed A teams index whose ability will be fixed as 0 (usually the team loss most which can be
#' generated using function BTdataframe).
#' @param thersh Thershold for convergency
#' @param max Maximum weight for w_{ij} (weight used for Adaptive Lasso)
#' @param iter Number of iterations used in L-BFGS-B algorithm.
#' @details Model selection through AIC or BIC method
#' @return
#' \item{Optimal.Likelihood}{Lowest AIC or BIC score}
#' \item{Optimal.degree}{The degree of freedom where lowest AIC or BIC score is achieved}
#' \item{Optimal.ability}{The ability where lowest AIC or BIC score is achieved}
#' \item{Likelihood}{Sequence contains all likelihood computed in this algorithm}
#' \item{degree}{Sequence contains all degrees computed in this algorithm}
#' \item{ability}{Matrix contains all abilities computed in this algorithm}
#' \item{lambda.min}{The lambda where lowest score is attained}
#' \item{penalty.min}{The penalty (s/\eqn{\max(s)}) where lowest score is attained}
#' @examples 
#' ##Initializing Dataframe
#' x <- BTdataframe(NFL2010)
#' 
#' ##Model selection through AIC
#' w <- BTLasso.weight(NFL$df, NFL$ability)
#' z <- BTdecayLassoC(NFL$df, NFL$ability, w, criteria = "AIC")
#' @export

BTdecayLassoC <- function(dataframe, ability, weight = NULL, lambda = NULL, criteria = "AIC", decay.rate = 0, fixed = 1, thersh = 1e-5, iter = 100) {
  Lp <- BTdecayLasso(dataframe, ability, lambda = lambda, weight = weight, path = TRUE, decay.rate = decay.rate,
                     fixed = fixed, thersh = thersh, max = max, iter = iter)
  
  if (criteria == "AIC") {
    x <- 2* Lp$df.path + degree * 2
  } else if (criteria == "BIC") {
    s0 <- 2 * s0 + degree * log(n1)
  } else {
    stop("Please specify AIC or BIC for criteria")
  }
  
  flag <- 0
  if (is.null(Lambda)) {
    k <- 0.25
    m <- seq(-6, -0.5, k)
    Lambda <- exp(m)
  } else {
    flag <- 1
  }
  
  n2 <- length(Lambda)
  ability0 <- ability[, -1]
  s <- c()
  n1 <- nrow(dataframe)
  d <- c()
  
  for (i in 1:n2) {
    BT <- BTdecayLasso.step2(dataframe, ability, Lambda[i], weight, decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
    degree <- BT$df
    s0 <- BTLikelihood(dataframe, BT$ability, decay.rate = decay.rate)
    
    if (criteria == "AIC") {
      s0 <- 2* s0 + degree * 2
    } else if (criteria == "BIC") {
      s0 <- 2 * s0 + degree * log(n1)
    } else {
      stop("Please specify AIC or BIC for criteria")
    }
    
    s <- c(s, s0)
    ability0 <- cbind(ability0, BT$ability)
    d <- c(d, degree)
  }
  
  j <- which.min(s)
  
  m0 <- m[j - 1]
  m1 <- m[j]
  s1 <- s[j]
  
  if (flag == 1) {
    k <- m1 - m0
  }
  
  while (k > thersh * 100) {
     k <- 0.5 * k
     m2 <- m1 - k
     BT <- BTdecayLasso.step2(dataframe, ability, exp(m2), weight, decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
     s0 <- BTLikelihood(dataframe, BT$ability, decay.rate = decay.rate)
     if (criteria == "AIC") {
       s0 <- 2* s0 + degree * 2
     } else if (criteria == "BIC") {
       s0 <- 2 * s0 + degree * log(n1)
     } else {
       stop("Please specify AIC or BIC for criteria")
     }
     degree <- BT$df
     
     if (s0 < s1) {
       m1 <- m2
       s1 <- s0
     }
     ability0 <- cbind(ability0, BT$ability)
     
     s <- c(s, s0)
     d <- c(d, degree)
  }
  
  output <- list(Optimal.Likelood = s0, Optimal.degree = degree, Optimal.ability = BT$ability,
                 Likelihood = s, degree = d, ability = ability0, lambda.min = m2, penalty.min = BT$penalty)
  output
}

