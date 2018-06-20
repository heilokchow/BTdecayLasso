#' @export

summary.slasso <- function(object, ...) {
  y <- object$ability
  z <- object$HYBRID.ability
  for (i in 1:ncol(y)) {
    y[, i] <- round(y[, i], 3)
    z[, i] <- round(z[, i], 3)
  }
  cat("penalty:\n\n")
  cat(object$penalty)
  cat("\n\ndecay.rate:\n\n")
  cat(object$decay.rate)
  cat("\n\nLasso Estimates:\n\n")
  print(y)
  cat("\nHYBRID Lasso Estimates:\n\n")
  print(z)
}

#' @export

summary.swlasso <- function(object, ...) {
  y <- object$ability
  z <- object$HYBRID.ability
  for (i in 1:ncol(y)) {
    y[, i] <- round(y[, i], 3)
    z[, i] <- round(z[, i], 3)
  }
  cat("penalty:\n\n")
  cat(object$penalty)
  cat("\n\ndecay.rate:\n\n")
  cat(object$decay.rate)
  cat("\n\nLasso Estimates:\n\n")
  print(y)
  cat("\n")
  cat("HYBRID Lasso Estimates:\n\n")
  print(z)
}

#' @export

summary.wlasso <- function(object, ...) {
  cat("Degree Path:\n\n")
  y <- data.frame(Run = seq(1, length(object$df.path), 1), Degree = object$df.path)
  print(y)
}


#' @export
 
summary.BTC <- function(object, ...) {
  y <- object$Optimal.ability
  for (i in 1:ncol(y)) {
    y[, i] <- round(y[, i], 3)
  }
  cat("type:\n\n")
  cat(object$type)
  cat("\n\npenalty:\n\n")
  cat(object$Optimal.penalty)
  cat("\n\ndegree:\n\n")
  cat(object$Optimal.degree)
  cat("\n\ndecay.rate:\n\n")
  cat(object$decay.rate)
  cat("\n\nEstimates:\n\n")
  print(y)
}

#' @export

summary.boot <- function(object, ...) {
  y <- object$Lasso
  z <- object$HYBRID.Lasso
  for (i in 1:ncol(y)) {
    y[, i] <- round(y[, i], 3)
    z[, i] <- round(z[, i], 3)
  }
  cat("penalty:\n\n")
  cat(object$penalty)
  cat("\n\ndecay.rate:\n\n")
  cat(object$decay.rate)
  cat("\n\nLasso Estimates:\n\n")
  print(y)
  cat("\n")
  cat("HYBRID Lasso Estimates:\n\n")
  print(z)
}

#' @export

summary.BT <- function(object, ...) {
  y <- object$ability
  for (i in 1:ncol(y)) {
    y[, i] <- round(y[, i], 3)
  }
  cat("decay.rate:\n\n")
  cat(object$decay.rate)
  cat("\n\nBradley-Terry Estimates:\n\n")
  print(y)
}

#' @export

summary.BTF <- function(object, ...) {
  y <- object$ability
  for (i in 1:ncol(y)) {
    y[, i] <- round(y[, i], 3)
  }
  cat("penalty:\n\n")
  cat(object$penalty)
  cat("\n\ndecay.rate:\n\n")
  cat(object$decay.rate)
  cat("\n\nBradley-Terry Lasso Estimates:\n\n")
  print(y)
  cat("\nlambda:\n\n")
  cat(object$lambda)
}