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
#' @param weight Weight for Lasso penalty on different abilities.
#' @param path whether the whole Lasso path will be run (plot.BTdecayLasso is enabled only if path = TRUE)
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
#' This likelihood function is optimized using L-BFGS-B method with package \bold{optimr}. Without specifying path = FALSE, the whole lasso path will be run
#' and we can use plot.BTdecayLasso to do a lasso path plot
#' @return
#' \item{ability}{Estimated ability scores with user given lambda}
#' \item{likelihood}{Negative likelihood of objective function with user given lambda}
#' \item{df}{Degree of freedom with user given lambda(number of distinct \eqn{\mu})}
#' \item{penalty}{\eqn{s/max(s)} with user given lambda}
#' \item{Lambda}{User given lambda}
#' \item{ability.path}{Estimated ability scores on whole Lasso path}
#' \item{likelihood.path}{Negative likelihood of objective function on whole Lasso path}
#' \item{df.path}{Degree of freedom on whole Lasso path(number of distinct \eqn{\mu})}
#' \item{penalty.path}{\eqn{s/max(s)} on whole Lasso path}
#' \item{Lambda.path}{Whole Lasso path}
#' \item{path}{Whether whole Lasso path will be run}
#' @examples
#' ##Initializing Dataframe
#' x <- BTdataframe(NFL2010)
#' 
#' ##BTdecayLasso run with exponential decay rate 0.005 and lambda 0.1 on whole lasso path using adaptive lasso
#' y1 <- BTdecayLasso.step2(x$df, x$ability, lambda = 0.1, decay.rate = 0.005, fixed = x$worstTeam)
#' 
#' ##Defining equal weight
#' n <- nrow(x$ability) - 1
#' w2 <- matrix(1, nrow = n, ncol = n)
#' w2[lower.tri(w2, diag = TRUE)] <- 0
#' 
#' ##BTdecayLasso run with exponential decay rate 0.005 and with a specific lambda 0.1
#' y2 <- BTdecayLasso.step2(x$df, x$ability, lambda = 0.1, weight = w2, path = FALSE, decay.rate = 0.005, fixed = x$worstTeam)
#' @export

BTdecayLasso <- function(dataframe, ability, lambda = NULL, weight = NULL, path = TRUE, decay.rate = 0, fixed = 1, thersh = 1e-5, max = 100, iter = 100) {
  u <- decay.rate
  n <- nrow(ability) - 1
  theta <- matrix(0, nrow = n, ncol = n) 
  Lagrangian <- matrix(0, nrow = n, ncol = n) 
  ability[, 1] <- 0
  
  if (path == TRUE){
    Lambda <- exp(seq(-0.5, -6, -0.0625))
    if (!is.null(lambda)){
      Lambda <- sort(c(Lambda, lambda), decreasing = TRUE)
    }
  } else if (is.null(lambda)) {
    stop("Please provide a sequence of lambda")
  } else {
    Lambda <- sort(lambda, decreasing = TRUE)
  }
  n2 <- length(Lambda)
  
  if (is.null(weight)) {
    weight <- BTLasso.weight(dataframe, ability, decay.rate = decay.rate, fixed = fixed, thersh = thersh, max = max, iter = iter)
  }
  
  v <- 1 
  ability0 <- ability[, -1]
  Hability0 <- ability[, -1]
  l <- c()
  hl <- c()
  p <- c()
  
  BT <- BTdecay(dataframe, ability, decay.rate = decay.rate, fixed = fixed, iter = iter)
  ability1 <- BT$ability
  s1 <- penaltyAmount(ability1, weight)
  l1 <- BTLikelihood(dataframe, ability1, decay.rate = decay.rate)
  df <- c()
  
  for (i in 1:n2) {
    stop <- 0
    j <- 1
    
    while (stop==0) {
      ability <- BTdecayLasso.step1(dataframe, ability, weight, Lagrangian, theta, v, Lambda[i], 
                                    decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
      theta <- BTtheta(ability, weight, Lagrangian, v, Lambda[i])
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
  }
  
  ability0 <- cbind(ability0, ability1)
  Hability0 <- cbind(Hability0, ability1)
  l <- c(l, l1)
  hl <- c(hl, l1)
  p <- c(p, 1)
  df <- c(df, n)
  
  if (is.null(lambda)) {
    output <- list(ability.path = ability0, likelihood.path = l, penalty.path = p, df.path = df, Lamda.path = c(Lambda, 0), path = path,
                   HYBRID.ability = Hability0, HYBRID.likelihood = hl)
    
  } else { 
    if (path == FALSE) {
      Lambda2 <- sort(lambda)
      x2 <- sapply(Lambda2, function(x) which(Lambda == x))
      ability2 <- ability0[, x2]
      l2 <- l[x2]
      p2 <- p[x2]
      df2 <- df[x2]
      output <- list(ability = ability2, likelihood = l2, penalty = p2, df = df2, Lambda = Lambda2, path = path,
                     HYBRID.ability = Hability0, HYBRID.likelihood = hl)
    } else {
      Lambda2 <- sort(lambda)
      x2 <- sapply(Lambda2, function(x) which(Lambda == x))
      ability2 <- ability0[, x2]
      l2 <- l[x2]
      p2 <- p[x2]
      df2 <- df[x2]
      output <- list(ability = ability2, likelihood = l2, penalty = p2, df = df2, Lambda = Lambda2,
                     ability.path = ability0, likelihood.path = l, penalty.path = p, df.path = df, Lambda.path = c(Lambda, 0), path = path,
                     HYBRID.ability = Hability0, HYBRID.likelihood = hl)
    }
  }
  output
}