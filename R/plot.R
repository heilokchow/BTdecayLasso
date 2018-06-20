#' Plot the Lasso path
#' 
#' @description Plot the whole lasso path run by BTdecayLasso() with given lambda and path = TRUE
#' @usage
#' ##S3 method for class "swlasso"
#' @param x Object with class "swlasso"
#' @param ... Further arguments pass to or from other methods
#' @export
#' @import ggplot2

plot.swlasso <- function(x, ...) {
  n <- nrow(x$ability.path) - 1
  df1 <- data.frame(ability = x$ability.path[1:n, 1], team = seq(1, n, 1), penalty = x$penalty.path[1])
  for (i in 1:(length(x$likelihood.path) - 1)) {
    df1 <- rbind(df1, data.frame(ability = x$ability.path[1:n, (i + 1)], team = seq(1, n, 1), penalty = x$penalty.path[i + 1]))
  }
  ggplot(df1, aes(x = penalty, y = ability, color = team)) + geom_line(aes(group = team))
  
}

#' Plot the Lasso path
#' 
#' @description Plot the whole lasso path run by BTdecayLasso() with lambda = NULL and path = TRUE
#' @usage
#' ##S3 method for class "wlasso"
#' @param x Object with class "wlasso"
#' @param ... Further arguments pass to or from other methods
#' @export
#' @import ggplot2

plot.wlasso <- function(x, ...) {
  n <- nrow(x$ability.path) - 1
  df1 <- data.frame(ability = x$ability.path[1:n, 1], team = seq(1, n, 1), penalty = x$penalty.path[1])
  for (i in 1:(length(x$likelihood.path) - 1)) {
    df1 <- rbind(df1, data.frame(ability = x$ability.path[1:n, (i + 1)], team = seq(1, n, 1), penalty = x$penalty.path[i + 1]))
  }
  ggplot(df1, aes(x = penalty, y = ability, color = team)) + geom_line(aes(group = team))
  
}