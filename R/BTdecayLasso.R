#' Bradley-Terry Model with Exponential Decayed weighted likelihood and Adaptive Lasso
#' 
#' @description 
#' Bradley-Terry model is applied for paired comparison data. Teams' ability score is estimated by maximizing log-likelihood function.
#' 
#' To achieve a better track of current abilities, we apply an exponential decay rate to weight the log-likelihood function.
#' The most current matches will weight more than previous matches. Parameter "decay.rate" in most functions of this package is used
#' to set the amount of exponential decay rate. decay.rate should be non-negative and the appropriate range of it depends on time scale in original dataframe.
#' (see \code{\link{BTdataframe}} and parameter "dataframe"'s definition of fifth column) For example,
#' a unit of week with a "decay.rate" 0.007 is equivalent to the unit of day with "decay.rate" 0.001. Usually, for sports matches,
#' if we take the unit of day, it's ranging from 0 to 0.01. The higher choice of "decay.rate", the better track of current teams' ability
#' with a side effect of higher variance.
#' 
#' If "decay.rate" is too large, for example "0.1" with a unit of day, \eqn{\exp(-0.7)} = 0.50. Only half weight will be add to the likelihood for matches played
#' one week ago and \eqn{\exp(-3.1)} = 0.05 suggests that previous matches took place one month ago will have little effect. Therefore, Only a few matches are
#' accounted for ability's estimation. It will lead to a very high variance and uncertainty. Since standard Bradley-Terry model
#' can not handle the case where there is a team who wins or loses all matches, such estimation may not provide convergent results. 
#' Thus, if our estimation provides divergent result, an error will be returned and we suggest user to chose a smaller "decay.rate"
#' or adding more match results into the same modeling period.
#' 
#' By default, the Adaptive Lasso is inplemented for variance reduction and team's grouping. Adaptive Lasso is proved to have good grouping property.
#' Apart from adaptive lasso, user can define own weight for different
#' Lasso constriant \eqn{\left|\mu_{i}-\mu_{j}\right|} where \eqn{\mu_{i}} is team i's ability.
#' 
#' Also by default, the whole Lasso path will be run. Similar to package "glmnet", user can provide their own choice of Lasso penalty "lambda" and determine whether the
#' whole Lasso path will be run (since such run is time-consuming). However, we suggest that if user is not familiar with the actual relationship among
#' lambda, the amount of penalty, the amount of shrinkage and grouping effect, a whole Lasso path should be run and selection of an
#' approperiate lambda is done by AIC or BIC criteria using \code{\link{BTdecayLassoC}} (since this model is time related, cross-validation method cannot be applied). Also, users can
#' use \code{\link{BTdecayLassoF}} to run with a specific Lasso penalty ranging from 0 to 1 (1 penalty means all estimators will shrink to 0).
#' 
#' Two sets of estimated abilities will be given, the biased Lasso estimation and the HYBRID Lasso's estimation.
#' HYBRID Lasso estimation solves the restricted Maximum Likelihood optimization based on the group determined by Lasso's estimation (Different team's ability will converges to
#' the same value if Lasso penalty is added and these teams' ability is setting to be equal as a restriction).
#' 
#' In addition, summary() using S3 method can be applied to view the outputs.
#' 
#' @param dataframe Generated using \code{\link{BTdataframe}} given raw data.
#' @param ability A column vector of teams ability, the last row is the home parameter.
#' The row number is consistent with the team's index shown in dataframe. It can be generated using \code{\link{BTdataframe}} given raw data.
#' @param lambda The amount of Lasso penalty induced. The input should be a positive scalar or a sequence.
#' @param weight Weight for Lasso penalty on different abilities.
#' @param path whether the whole Lasso path will be run (plot.BTdecayLasso is enabled only if path = TRUE)
#' @param decay.rate A non-negative exponential decay rate. Usually ranging from (0, 0.01), A larger decay rate weights more
#' importance to most recent matches and the estimated parameters reflect more on recent behaviour.
#' @param fixed A teams index whose ability will be fixed as 0. The worstTeam's index
#' can be generated using \code{\link{BTdataframe}} given raw data.
#' @param thersh Threshold for convergency used for Augmented Lagrangian Method.
#' @param max Maximum weight for w_{ij} (weight used for Adaptive Lasso)
#' @param iter Number of iterations used in L-BFGS-B algorithm.
#' @details
#' According to \code{\link{BTdecay}}, the objective likelihood function to be optimized is,
#' \deqn{\sum_{k=1}^{n}\sum_{i<j}\exp(-\alpha t_{k})\cdot(y_{ij}(\tau h_{ij}^{t_{k}}+\mu_{i}-\mu_{j})-\log(1+\exp(\tau h_{ij}^{t_{k}}+\mu_{i}-\mu_{j})))}
#' The Lasso constraint is given as,
#' \deqn{\sum_{i<j}w_{ij}\left|\mu_{i}-\mu_{j}\right|\leq s}
#' where \eqn{w_{ij}} are predefined weight. For Adaptive Lasso, \eqn{\left|w_{ij}=1/(\mu_{i}^{MLE}-\mu_{j}^{MLE})\right|}.
#' 
#' Maximize this constraint objective function is equivalent to minimizing the following equation,
#' \deqn{-l(\mu,\tau)+\lambda\sum_{i<j}w_{ij}|\mu_{i}-\mu_{j}|}
#' Where \eqn{-l(\mu,\tau)} is taking negative value of objective function above.  Increase "lambda" will decrease "s", their relationship is
#' monotone. Here, we define "penalty" as \eqn{1-s/\max(s)}. Thus, "lambda" and "penalty" has a positive correlation.
#' @return
#' \item{ability}{Estimated ability scores with user given lambda}
#' \item{likelihood}{Negative likelihood of objective function with user given lambda}
#' \item{df}{Degree of freedom with user given lambda(number of distinct \eqn{\mu})}
#' \item{penalty}{\eqn{s/max(s)} with user given lambda}
#' \item{Lambda}{User given lambda}
#' \item{ability.path}{if path = TRUE, estimated ability scores on whole Lasso path}
#' \item{likelihood.path}{if path = TRUE, negative likelihood of objective function on whole Lasso path}
#' \item{df.path}{if path = TRUE, degree of freedom on whole Lasso path(number of distinct \eqn{\mu})}
#' \item{penalty.path}{if path = TRUE, \eqn{s/max(s)} on whole Lasso path}
#' \item{Lambda.path}{if path = TRUE, Whole Lasso path}
#' \item{path}{Whether whole Lasso path will be run}
#' \item{HYBRID.ability.path}{If path = TRUE, the whole path of evolving of HYBRID ability}
#' \item{HYBRID.likelihood.path}{if path = TRUE, the whole path of HYBRID likelihood}
#' @seealso \code{\link{BTdataframe}} for dataframe initialization,
#' \code{\link{plot.swlasso}},  \code{\link{plot.wlasso}} are used for Lasso path plot if path = TRUE in this function's run
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
#' ##The following code runs the main results
#' ##But they will not be run in R CMD check since these iterations are time-consuming
#' ##Usually a single lambda's run will take 1-20 s
#' ##The whole Adaptive Lasso run will take 5-20 min
#' \dontrun{
#' ##BTdecayLasso run with exponential decay rate 0.005 and 
#' ##lambda 0.1 on whole lasso path using adaptive lasso
#' y1 <- BTdecayLasso(x$dataframe, x$ability, lambda = 0.1, 
#'                    decay.rate = 0.005, fixed = x$worstTeam)
#' summary(y1)
#' 
#' ##Defining equal weight
#' ##Note that comparing to Adaptive weight, the user defined weight may not be 
#' ##efficient in groupiing. Therefore, to run the whole Lasso path 
#' ##(evolving of distinct ability scores), it may take a much longer time. 
#' ##We recommend the user to apply the default setting,
#' ##where Adaptive Lasso will be run.
#' n <- nrow(x$ability) - 1
#' w2 <- matrix(1, nrow = n, ncol = n)
#' w2[lower.tri(w2, diag = TRUE)] <- 0
#' 
#' ##BTdecayLasso run with exponential decay rate 0.005 and with a specific lambda 0.1
#' y2 <- BTdecayLasso(x$dataframe, x$ability, lambda = 0.1, weight = w2, 
#'                    path = FALSE, decay.rate = 0.005, fixed = x$worstTeam)
#' 
#' ##BTdecayLasso run with exponential decay rate 0.005 and with a specific lambda 0.1
#' ##Time-consuming
#' y3 <- BTdecayLasso(x$dataframe, x$ability, lambda = 0.1, weight = w2, 
#'                    path = TRUE, decay.rate = 0.005, fixed = x$worstTeam)
#' summary(y2)
#' 
#' ##Plot the Lasso path (S3 method)
#' plot(y1)
#' plot(y3)
#' }
#' 
#' @export

BTdecayLasso <- function(dataframe, ability, lambda = NULL, weight = NULL, path = TRUE, decay.rate = 0, fixed = 1, thersh = 1e-5, max = 100, iter = 100) {
  
  
  u <- decay.rate
  n <- nrow(ability) - 1
  theta <- matrix(0, nrow = n, ncol = n) 
  Lagrangian <- matrix(0, nrow = n, ncol = n) 
  ability[, 1] <- 0
  k0 <- 0.0675
  
  if(!(fixed %in% seq(1, n, 1))){
    stop("The fixed team's index must be an integer index of one of all teams")
  }
  
  if (is.null(weight)) {
    weight <- BTLasso.weight(dataframe, ability, decay.rate = decay.rate, fixed = fixed, thersh = thersh, max = max, iter = iter)
  }
  
  v <- 1 
  ability0 <- ability[, -1]
  Hability0 <- ability[, -1]
  l <- c()
  hl <- c()
  p <- c()
  slambda <- c()
  
  BT <- BTdecay(dataframe, ability, decay.rate = decay.rate, fixed = fixed, iter = iter)
  ability1 <- BT$ability
  s1 <- penaltyAmount(ability1, weight)
  l1 <- BTLikelihood(dataframe, ability1, decay.rate = decay.rate)
  df <- c()
  j <- 1
  
  if (path == FALSE && is.null(lambda)) {
    stop("Please provide a sequence of lambda or enable lasso path")
  }
  
  if (path == TRUE) {
    lambda0 <- 1
    lambda1 <- exp(-0.5)
    degree0 <- 0
    
    while (degree0 < (n - 1)) {
      stop <- 0
      while (stop==0) {
        ability <- BTdecayLasso.step1(dataframe, ability, weight, Lagrangian, theta, v, lambda1, 
                                      decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
        theta <- BTtheta(ability, weight, Lagrangian, v, lambda1)
        Lagrangian0 <- BTLagrangian(Lagrangian, ability, theta, v)
        k <- sum(abs(Lagrangian0 - Lagrangian))
        if (k < thersh) {
          stop <- 1
        } else {
          Lagrangian <- Lagrangian0
          v <- max(Lagrangian^2)
        }
        s0 <- penaltyAmount(ability, weight)
      }
      
      cat(s0, '\n')
      p0 <- s0/s1
      ability0 <- cbind(ability0, ability)
      l0 <- BTLikelihood(dataframe, ability, decay.rate = decay.rate)
      l <- c(l, l0)
      p <- c(p, p0)
      
      degree1 <- round(ability[1:n, 1], -log10(thersh)-1)
      degree <- length(unique(degree1))
      df <- c(df, degree)
      
      map <- function(x){
        if (x == 1){
          return(1)
        } else {
          match <- which(degree1[1:(x - 1)] == degree1[x])
          if (length(match) == 0) {
            return(length(unique(degree1[1:x])))
          } else {
            return(length(unique(degree1[1:match[1]])))
          }
        }
      }
      
      dataframe1 <- dataframe
      dataframe1[, 1] <- sapply(dataframe[, 1], map)
      dataframe1[, 2] <- sapply(dataframe[, 2], map)
      Hability1 <- matrix(0, nrow = (degree + 1), ncol = 1)
      HBT <- BTdecay(dataframe1, Hability1, decay.rate = decay.rate, fixed = map(fixed), iter = iter)
      Hability1 <- HBT$ability
      Hability <- ability
      for (i in 1:n) {
        Hability[i, 1] <- Hability1[map(i), 1]
      }
      Hability[(n + 1), 1] <- Hability1[(degree + 1), 1]
      Hability0 <- cbind(Hability0, Hability)
      hl0 <- BTLikelihood(dataframe, Hability, decay.rate = decay.rate)
      hl <- c(hl, hl0)
      
      slambda <- c(slambda, lambda1)
      
      if (degree > (degree0 + 1) && abs(lambda0 - lambda1) > thersh) {
        lambda1 <- (lambda0 + lambda1)/2
      } else {
        lambda0 <- lambda1
        lambda1 <- lambda1 * exp(-k0)
        degree0 <- max(degree0, degree)
      }
    }
  }
  
  if (!is.null(lambda)){
    lambda <- sort(lambda, decreasing = TRUE)
    for (i in 1:length(lambda)) {
      stop <- 0
      while (stop==0) {
        ability <- BTdecayLasso.step1(dataframe, ability, weight, Lagrangian, theta, v, lambda[i], 
                                      decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
        theta <- BTtheta(ability, weight, Lagrangian, v, lambda[i])
        Lagrangian0 <- BTLagrangian(Lagrangian, ability, theta, v)
        k <- sum(abs(Lagrangian0 - Lagrangian))
        if (k < thersh) {
          stop <- 1
        } else {
          Lagrangian <- Lagrangian0
          v <- max(Lagrangian^2)
        }
        s0 <- penaltyAmount(ability, weight)
      }
      
      cat(s0, '\n')
      p0 <- s0/s1
      ability0 <- cbind(ability0, ability)
      l0 <- BTLikelihood(dataframe, ability, decay.rate = decay.rate)
      l <- c(l, l0)
      p <- c(p, p0)
      
      degree1 <- round(ability[1:n, 1], -log10(thersh)-1)
      degree <- length(unique(degree1))
      df <- c(df, degree)
      
      map <- function(x){
        if (x == 1){
          return(1)
        } else {
          match <- which(degree1[1:(x - 1)] == degree1[x])
          if (length(match) == 0) {
            return(length(unique(degree1[1:x])))
          } else {
            return(length(unique(degree1[1:match[1]])))
          }
        }
      }
      
      dataframe1 <- dataframe
      dataframe1[, 1] <- sapply(dataframe[, 1], map)
      dataframe1[, 2] <- sapply(dataframe[, 2], map)
      Hability1 <- matrix(0, nrow = (degree + 1), ncol = 1)
      HBT <- BTdecay(dataframe1, Hability1, decay.rate = decay.rate, fixed = map(fixed), iter = iter)
      Hability1 <- HBT$ability
      Hability <- ability
      for (i in 1:n) {
        Hability[i, 1] <- Hability1[map(i), 1]
      }
      Hability[(n + 1), 1] <- Hability1[(degree + 1), 1]
      Hability0 <- cbind(Hability0, Hability)
      hl0 <- BTLikelihood(dataframe, Hability, decay.rate = decay.rate)
      hl <- c(hl, hl0)
      
      slambda <- c(slambda, lambda[i])
    }
  }
  
  
  ability0 <- cbind(ability0, ability1)
  Hability0 <- cbind(Hability0, ability1)
  l <- c(l, l1)
  hl <- c(hl, l1)
  p <- c(p, 1)
  df <- c(df, n)
  
  if (is.null(lambda)) {
    output <- list(ability.path = ability0, likelihood.path = l, penalty.path = p, df.path = df, Lambda.path = c(slambda, 0), path = path,
                   HYBRID.ability.path = Hability0, HYBRID.likelihood.path = hl, decay.rate = decay.rate)
    class(output) <- "wlasso"
    
  } else { 
    n3 <- length(lambda)
    n4 <- length(slambda)
    if (path == FALSE) {
      output <- list(ability = as.matrix(ability0[, 1:n4]), likelihood = l[1:n4], penalty = p[1:n4], df = df[1:n4], Lambda = slambda, path = path,
                     HYBRID.ability = as.matrix(Hability0[, 1:n4]), HYBRID.likelihood = hl[1:n4], decay.rate = decay.rate)
      class(output) <- "slasso"
    } else {
      output <- list(ability = as.matrix(ability0[, (n4 - n3 + 1):n4]), likelihood = l[(n4 - n3 + 1):n4], penalty = p[(n4 - n3 + 1):n4], df = df[(n4 - n3 + 1):n4], Lambda = lambda,
                     ability.path = ability0, likelihood.path = l, penalty.path = p, df.path = df, Lambda.path = c(slambda, 0), path = path,
                     HYBRID.ability = as.matrix(Hability0[, (n4 - n3 + 1):n4]), HYBRID.likelihood = hl[(n4 - n3 + 1):n4], 
                     HYBRID.ability.path = Hability0, HYBRID.likelihood.path = hl,
                     decay.rate = decay.rate)
      class(output) <- "swlasso"
    }
  }
  output
}