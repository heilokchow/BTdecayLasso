#' Compute the standard deviation of Bradley-Terry decay Lasso model by bootstrapping
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
#' @param lambda The amount of Lasso penalty induced, only a single scalar is accepted in bootstrapping.
#' @param boot Amount of simulations.
#' @param weight Weight for Lasso penalty on different abilities.
#' @param decay.rate The exponential decay rate. Usually ranging from (0, 0.1), A larger decay rate weights more
#' importance to most recent matches and the estimated parameters reflect more on recent behaviour.
#' @param fixed A teams index whose ability will be fixed as 0 (usually the team loss most which can be
#' generated using function BTdataframe).
#' @param thersh Thershold for convergency
#' @param max Maximum weight for w_{ij} (weight used for Adaptive Lasso)
#' @param iter Number of iterations used in L-BFGS-B algorithm.
#' @details 100 times of simulation will be done by default, user can adjust the numbers of simulation by input of boot. However, bootstrapping process
#' is time consuming and in order to produce stable result, 1000 times of simulations is enough.
#' @return A list contain Lasso and Hybrid Lasso's bootstrapping's mean and standard deviation.
#' \item{Lasso}{Lasso bootstrapping's result. A three column matrix where first column is the original dataset's estimation, the second column is bootstrapping mean and the last column is the
#' bootstrapping standard deviation}
#' \item{HYBRID.Lasso}{HYBRID Lasso bootstrapping's result. A three column matrix where the first column is the original dataset's estimation, the second column is bootstrapping mean and the last column is the
#' bootstrapping standard deviation}
#' @examples
#' ##Initialize the dataframe and ability
#' NFL <- BTdataframe(NFL2010)
#' 
#' ##Run Lasso estimate for whole Lasso path
#' z <- BTdecayLasso(NFL$dataframe, NFL$ability, fixed = NFL$worstTeam)
#' 
#' ##Model selection using AIC with Lasso's likelihood
#' z1 <- BTdecayLassoC(NFL$dataframe, NFL$ability, model = z, criteria = "AIC", type = "LASSO", fixed = NFL$worstTeam)
#' 
#' ##Bootstrapping for model with lowest AIC score for 100 times.
#' z2 <- boot.BTdecayLasso(NFL$dataframe, NFL$ability, lambda = z1$Optimal.lambda, boot = 100, fixed = NFL$worstTeam)
#' @export

boot.BTdecayLasso <- function(dataframe, ability, lambda, boot = 100, weight = NULL, decay.rate = 0, fixed = 1,
                              thersh = 1e-5, max = 100, iter = 100) {
  BT <- BTdecay(dataframe, ability, decay.rate = decay.rate, fixed = fixed, iter = iter)
  Tability <- BT$ability
  Bability <- Tability[, -1]
  Hability <- Tability[, -1]
  n1 <- nrow(dataframe)
  n <- nrow(ability) -1
  
  y1 <- sapply(dataframe[, 1], function(x) Tability[x, 1])
  y2 <- sapply(dataframe[, 2], function(x) Tability[x, 1])
  t <- exp(Tability[n + 1] + y1 - y2)
  y <- t/(1 + t)
  
  set.seed(271)
  dataframe1 <- dataframe
  for (i in 1:boot) {
    
    stop <- 0
    while (stop == 0) {
      r <- runif(n1)
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
  Bsd <- matrix(apply(Bability, 1, sd))
  Hmean <- matrix(apply(Hability, 1, mean))
  Hsd <- matrix(apply(Hability, 1, sd))
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

