#' @export
BTdataframe <- function(dataframe){
  
  
  ## Intitialize dataframe of match results for BTdecayLasso function
  df1 <- dataframe
  df <- matrix(ncol = ncol(df1), nrow = nrow(df1))
  team <- as.character(unique(df1[, 1]))
  n <- length(team)
  
  for(i in 1:n){
    df[df1 == team[i]] <- i
  }
  
  df[, 3] <- df1[, 3]
  df[, 4] <- df1[, 4]
  df[, 5] <- df1[, 5] - df1[1, 5]
  
  ## Determine the team who has the worst performance
  k0 <- 1000
  i0 <- 1
  
  for (i in 1:n) {
    k1 <- sum(df[df[, 1] == i, 3]) - sum(df[df[, 1] == i, 4]) + sum(df[df[, 2] == i, 4]) - sum(df[df[, 2] == i, 3])
    if(k0 > k1){
      i0 <- i
      k0 <- k1
    }
  }
  
  ## Intitialize ablity vector for BTdecayLasso function
  ability <- matrix(0, ncol = 1, nrow = (n + 1))
  colnames(ability) <- c("score")
  rownames(ability) <- c(team, "at.home")
  
  output <- list(df = df, ability = ability, worstTeam = i0)
  output
}