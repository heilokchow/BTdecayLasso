#' Compute the standard deviation of Bradley-Terry decay Lasso model by bootstrapping
#' 
#' Bootstrapping is done assuming that Maximum Likelihood's estimation reflects the true abilities.
#' Same level of Lasso penalty "lambda" should be applied in different simulation models for Lasso induced estimation.
#' 
#' @param dataframe Generated using \code{\link{BTdataframe}} given raw data.
#' @param ability A column vector of teams ability, the last row is the home parameter.
#' The row number is consistent with the team's index shown in dataframe. It can be generated using \code{\link{BTdataframe}} given raw data.
#' @param lambda The amount of Lasso penalty induced, only a single scalar is accepted in bootstrapping.
#' @param boot Amount of simulations.
#' @param weight Weight for Lasso penalty on different abilities.
#' @param decay.rate The exponential decay rate. Usually ranging from (0, 0.01), A larger decay rate weights more
#' importance to most recent matches and the estimated parameters reflect more on recent behaviour.
#' @param fixed A teams index whose ability will be fixed as 0. The worstTeam's index
#' can be generated using \code{\link{BTdataframe}} given raw data.
#' @param thersh Threshold for convergency
#' @param max Maximum weight for \eqn{w_{ij}} (weight used for Adaptive Lasso).
#' @param iter Number of iterations used in L-BFGS-B algorithm.
#' @details 100 times of simulation will be done by default, user can adjust the numbers of simulation by input of boot. However, bootstrapping process
#' is time consuming and usually 1000 time of simulations is enough to provide a stable result.
#' 
#' More detailed description of "lambda", "penalty" and "weight" are documented in \code{\link{BTdecayLasso}}.
#' 
#' summary() function follows S3 method can be applied to view the outputs.
#' @return A list with class "boot" contain Lasso and Hybrid Lasso's bootstrapping's mean and standard deviation.
#' \item{Lasso}{Lasso bootstrapping's result. A three column matrix where first column is the original dataset's estimation, the second column is bootstrapping mean and the last column is the
#' bootstrapping standard deviation}
#' \item{HYBRID.Lasso}{HYBRID Lasso bootstrapping's result. A three column matrix where the first column is the original dataset's estimation, the second column is bootstrapping mean and the last column is the
#' bootstrapping standard deviation}
#' @seealso \code{\link{BTdataframe}} for dataframe initialization, \code{\link{BTdecayLasso}} for detailed description
#' @references 
#' Masarotto, G. and Varin, C.(2012) The Ranking Lasso and its Application to Sport Tournaments. 
#' *The Annals of Applied Statistics* **6** 1949--1970.
#' 
#' Zou, H. (2006) The adaptive lasso and its oracle properties. 
#' *J.Amer.Statist.Assoc* **101** 1418--1429.
#' @examples
#' ##Initialize the dataframe and ability
#' x <- BTdataframe(NFL2010)
#' 
#' ##The following code runs the main results
#' ##But they will not be run in R CMD check since these iterations are time-consuming
#' \dontrun{
#' ##Run Lasso estimate for whole Lasso path
#' z <- BTdecayLasso(x$dataframe, x$ability, fixed = x$worstTeam)
#' 
#' ##Model selection using AIC with Lasso's likelihood
#' z1 <- BTdecayLassoC(x$dataframe, x$ability, model = z, 
#'                     criteria = "AIC", type = "LASSO", fixed = x$worstTeam)
#' 
#' ##Bootstrapping for model with lowest AIC score for 100 times.
#' ##Note that the decay.rate used in model selection should be consistent with
#' ##the one which is used in whole Lasso path's run (keep the same model)
#' z2 <- boot.BTdecayLasso(x$dataframe, x$ability, lambda = z1$Optimal.lambda, 
#'                         boot = 100, fixed = x$worstTeam)
#' }
#' 
#' @export
#' @import stats

boot.BTdecayLasso <- function(dataframe, ability, lambda, boot = 100, weight = NULL, decay.rate = 0, fixed = 1,
                              thersh = 1e-5, max = 100, iter = 100) {
  
  boot <- round(boot)
  if (boot < 2) {
    stop("Boot should be an integer greater than 1")
  }
  
  BT <- BTdecay(dataframe, ability, decay.rate = decay.rate, fixed = fixed, iter = iter)
  Tability <- BT$ability
  Bability <- Tability[, -1]
  Hability <- Tability[, -1]
  n1 <- nrow(dataframe)
  n <- nrow(ability) -1
  
  if(!(fixed %in% seq(1, n, 1))){
    stop("The fixed team's index must be an integer index of one of all teams")
  }
  
  y1 <- sapply(dataframe[, 1], function(x) Tability[x, 1])
  y2 <- sapply(dataframe[, 2], function(x) Tability[x, 1])
  t <- exp(Tability[n + 1] + y1 - y2)
  y <- t/(1 + t)
  
  set.seed(271)
  dataframe1 <- dataframe
  for (i in 1:boot) {
    
    stop <- 0
    while (stop == 0) {
      r <- stats::runif(n1)
      dataframe1[, 3] <- as.numeric(r < y)
      dataframe1[, 4] <- 1 - dataframe1[, 3]
             
      counts <- matrix(0, nrow = n, ncol = 2)
      for (i in 1:n1) {
        a1 <- dataframe1[i, 1]
        a2 <- dataframe1[i, 2]
        counts[a1, 1] <- counts[a1, 1] + dataframe1[i, 3]
        counts[a2, 1] <- counts[a2, 1] + dataframe1[i, 4]
        counts[a1, 2] <- counts[a1, 2] + dataframe1[i, 3] + dataframe1[i, 4]
        counts[a2, 2] <- counts[a2, 2] + dataframe1[i, 3] + dataframe1[i, 4]
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
      if (is.null(win) && is.null(loss)) {
        stop <- 1
      } 
    }
    
    BTL <- BTdecayLasso(dataframe1, ability, lambda = lambda, decay.rate = decay.rate, path = FALSE, fixed = fixed, 
                        thersh = thersh, max = max, iter = iter)
    Bability <- cbind(Bability, BTL$ability)
    Hability <- cbind(Hability, BTL$HYBRID.ability)
  }
  
  Bmean <- matrix(apply(Bability, 1, mean))
  Bsd <- matrix(apply(Bability, 1, stats::sd))
  Hmean <- matrix(apply(Hability, 1, mean))
  Hsd <- matrix(apply(Hability, 1, stats::sd))
  BTL <- BTdecayLasso(dataframe, ability, lambda = lambda, decay.rate = decay.rate, path = FALSE, fixed = fixed, 
                      thersh = thersh, max = max, iter = iter)
  Bori <- BTL$ability
  Hori <- BTL$HYBRID.ability
  
  Bout <- cbind(Bori, Bmean, Bsd)
  colnames(Bout) <- c("Original", "Est.Mean", "Est.std")
  Hout <- cbind(Hori, Hmean, Hsd)
  colnames(Hout) <- c("Original", "Est.Mean", "Est.std")
  output <- list(penalty = BTL$penalty, Lasso = Bout, HYBRID.Lasso = Hout, decay.rate = decay.rate)
  class(output) <- "boot"
  output
}

