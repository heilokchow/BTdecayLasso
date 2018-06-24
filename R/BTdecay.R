#' Bradley-Terry Model with Exponential Decayed weighted likelihood
#'
#' Exponential decay rate is applied to the likelihood function to achieve a better track of current abilities. When "decay.rate" is setting as 0,
#' this is a standard Bradley-Terry Model whose estimated parameters are equivalent to package "BradleyTerry2".
#' Further detailed description is attached in \code{\link{BTdecayLasso}}.
#'
#' @param dataframe Generated using \code{\link{BTdataframe}} given raw data.
#' @param ability A column vector of teams ability, the last row is the home parameter.
#' The row number is consistent with the team's index shown in dataframe. It can be generated using \code{\link{BTdataframe}} given raw data.
#' @param decay.rate The exponential decay rate. Usually ranging from (0, 0.01), A larger decay rate weights more
#' importance to most recent matches and the estimated parameters reflect more on recent behaviour.
#' @param fixed A teams index whose ability will be fixed as 0. The worstTeam's index
#' can be generated using \code{\link{BTdataframe}} given raw data.
#' @param iter Number of iterations used in L-BFGS-B algorithm.
#' @details
#' The standard Bradley-Terry Model defines the winning probability of i against j,
#' \deqn{P(Y_{ij}=1)=\frac{\exp(\tau h_{ij}^{t_{k}}+\mu_{i}-\mu_{j})}{1+\exp(\tau h_{ij}^{t_{k}}+\mu_{i}-\mu_{j})}}
#' \eqn{\tau} is the home parameter and \eqn{\mu_{i}} is the team i's ability score. \eqn{h_{ij}} takes 1 if team i is at home, -1 otherwise.
#' Given, a complete tournament's result. The objective likelihood function with an exponential decay rate is,
#' \deqn{\sum_{k=1}^{n}\sum_{i<j}\exp(-\alpha t_{k})\cdot(y_{ij}(\tau h_{ij}^{t_{k}}+\mu_{i}-\mu_{j})-\log(1+\exp(\tau h_{ij}^{t_{k}}+\mu_{i}-\mu_{j})))}
#' where n is the number of matches, \eqn{\alpha} is the exponential decay rate and \eqn{y_{ij}} takes 0 if i is defeated by j, 1 otherwise. \eqn{t_{k}} is
#' the time lag (time until now). 
#' This likelihood function is optimized using L-BFGS-B method with package \bold{optimr} and summary() function with S3 method can be applied to view the outputs.
#' @return List with class "BT" contains estimated abilities and convergent code, 0 stands for convergence reaches,
#' 1 stands for convergence not reaches. If 1 is returned, we suggest that decay rate should be set lower.
#' Bradley-Terry model fails to model the situation when a team wins or loses in all matches.
#' If a high decay rate is considered, a team who only loses or wins 1 matches long time ago will also causes the same problem.
#' \item{ability}{Estimated ability scores}
#' \item{convergence}{0 stands for convergent, 1 stands for not convergent}
#' \item{decay.rate}{Decay rate of this model}
#' @examples 
#' ##Initializing Dataframe
#' x <- BTdataframe(NFL2010)
#' 
#' ##Standard Bradley-Terry Model optimization
#' y <- BTdecay(x$dataframe, x$ability, decay.rate = 0, fixed = x$worstTeam)
#' summary(y)
#' 
#' ##Dynamic approximation of current ability scores using exponential decayed likelihood.
#' ##If we take decay.rate = 0.005
#' ##Match happens one month before will weight exp(-0.15)=0.86 on log-likelihood function
#' z <- BTdecay(x$dataframe, x$ability, decay.rate = 0.005, fixed = x$worstTeam)
#' summary(z)
#' @import optimr
#' @export

BTdecay <- function(dataframe, ability, decay.rate = 0, fixed = 1, iter = 100){
  
  
  ## Initialize the parameters
  df <- as.matrix(dataframe)
  u <- decay.rate
  n1 <- nrow(df)
  n <- nrow(ability) - 1
  counts <- matrix(0, nrow = n, ncol = 2)
  
  if(!(fixed %in% seq(1, n, 1))){
    stop("The fixed team's index must be an integer index of one of all teams")
  }
  
  ## Check the validity of standard Bradley-Terry Model
  for (i in 1:n1) {
    a1 <- df[i, 1]
    a2 <- df[i, 2]
    counts[a1, 1] <- counts[a1, 1] + df[i, 3]
    counts[a2, 1] <- counts[a2, 1] + df[i, 4]
    counts[a1, 2] <- counts[a1, 2] + df[i, 3] + df[i, 4]
    counts[a2, 2] <- counts[a2, 2] + df[i, 3] + df[i, 4]
  }
  
  win <- c()
  loss <- c()
  for (i in 1:n) {
    if (counts[i, 1] == counts[i, 2]) {
      win <- c(win , i)
    } else if (counts[i, 1] == 0) {
      loss <- c(loss, i)
    }
  }
  if (!is.null(win)) stop('Bradley-Terry Model cannot deal with the case if there exists team who wins all matches')
  if (!is.null(loss)) stop('Bradley-Terry Model cannot deal with the case if there exists team who loses all matches')
  
  ## Iterations of the estimation of ability scores
  fn <- function(ability){
    s <- 0
    for(i in 1:n1){
      a1 <- df[i, 1]
      a2 <- df[i, 2]
      C <- exp(-u * df[i, 5])
      x <- ability[a1]
      y <- ability[a2]
      p <- ability[n + 1] + x - y
      q <- exp(p)
      s <- s - (df[i, 3] * p - (df[i, 3] + df[i, 4]) * log(q + 1)) * C
    }
    s
  }
  
  gr <- function(ability){
    Grad <- rep(0, n + 1)
    for(i in 1:n1){
      a1 <- df[i, 1]
      a2 <- df[i, 2]
      C <- exp(-u * df[i, 5])
      x <- ability[a1]
      y <- ability[a2]
      p <- ability[n + 1] + x - y
      q <- exp(p)
      A <- -(df[i, 3] * (1/(q + 1)) + df[i, 4] * (-q/(q + 1))) * C
      Grad[a1] <- Grad[a1] + A
      Grad[a2] <- Grad[a2] - A
      Grad[n + 1] <- Grad[n + 1] + A
    }
    Grad
  }
  
  xa <- optimr::optimr(rep(0, n + 1), fn, gr = gr, method = "L-BFGS-B", control = list(maxit = iter))
  
  ability[, 1] <- xa$par - xa$par[fixed]
  ability[n + 1, 1] <- xa$par[n + 1]
  output <- list(ability = ability, convergence = xa$convergence, decay.rate = decay.rate)
  
  if(xa$convergence == 1){
    stop("Iterations diverge, please provide a smaller decay rate or more data")
  }
  
  class(output) <- "BT"
  output
}
