#' Dataframe initialization
#' 
#' @details Initial the raw dataframe and return an un-estimated ability vector and the worst
#' team who loses most.
#' @param dataframe Raw dataframe input, an example data "NFL2010" is attached in package for reference
#' The raw data is a dataframe with 5 columns. First column is home teams.
#' Second column is away teams.
#' Third column is the number of wins of home teams (if home team defeats away team, record 1 here, 0 otherwise).
#' Fourth column is the number of wins of away teams (if home team defeats away team, record 0 here, 1 otherwise).
#' Fifth column is a scalar of time when the match is played until now (Time lag). Any time scale can be used here.
#' "NFL2010" applies the unit of day.
#' @param home Whether home effect will be considered, the default is TRUE.
#' @details Note that even if the tournament does not have any home team or away team, you can still provide the match results
#' according to the description above regardless of who is at home and who is away. By selecting the home = FALSE,
#' We duplicate the dataset, switch the home, away teams and also the home, away match results. Then this dataset will
#' be attached to the original dataset and all home and away win's number will be divided by 2. MLE estimation of home effect
#' is proved to be an exact 0.
#' 
#' The elimination of home effect by duplicating the original dataset will be less efficient than eliminating the home parameter
#' directly in iterations. Since most games such as football, basketball have home effect and this method provides an idea of
#' handling the case where some games have home effect and some games are played on neutral place, this method is applied here.
#' @return 
#' \item{dataframe}{dataframe for Bradley-Terry run}
#' \item{ability}{Initial ability vector for iterations}
#' \item{worstTeam}{The worst team whose ability can be set as 0 during any model's run}
#' @export
BTdataframe <- function(dataframe, home =TRUE) {
  
  
  ## Intitialize dataframe of match results for BTdecayLasso function
  df1 <- as.data.frame(dataframe)
  df <- matrix(ncol = ncol(df1), nrow = nrow(df1))
  team <- as.matrix(unique(df1[, 1]))
  n <- length(team)
  
  for(i in 1:n){
    df[df1 == team[i, 1]] <- i
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
  
  if (home == FALSE) {
    df2 <- df
    df2[, 1] <- df[, 2]
    df2[, 2] <- df[, 1]
    df2[, 3] <- df[, 4]
    df2[, 4] <- df[, 3]
    df <- rbind(df, df2)
    df[, 3] <- df[, 3]/2
    df[, 4] <- df[, 4]/2
  }
  
  output <- list(dataframe = df, ability = ability, worstTeam = i0)
  output
}