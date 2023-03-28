
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BTdecayLasso

[![CRAN](http://www.r-pkg.org/badges/version/BTdecayLasso)](https://cran.r-project.org/package=BTdecayLasso)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/BTdecayLasso)](https://cran.r-project.org/package=BTdecayLasso)
[![Last Commit](https://img.shields.io/github/last-commit/heilokchow/BTdecayLasso)](https://github.com/heilokchow/BTdecayLasso)

Bradley-Terry model is used for ranking in sports tournament. Given the
standard Bradley-Terry model, we use an exponential decay rate to weight
its log-likelihood function and apply Lasso penalty to achieve a
variance reduction and team grouping.

## Installation

You can install BTdecayLasso from github with:

``` r
# install.packages("devtools")
devtools::install_github("heilokchow/BTdecayLasso")
```

## Example

This is a basic example which shows you how to solve a common problem:

First, given raw datasets (five columns are home teams, away teams, home
wins, away wins, time until now), we convert this dataset into a
dataframe which can be used for other function’s input.

``` r
NFL <- BTdataframe(NFL2010)
```

Then, we comput the whole Lasso path for further analysis’s use. In this
example, to track the dynamically changing abilities, we set
‘decay.rate’ to be 0.005. A higher decay rate will give more unbiased
results for current abilites’ estimation with a side effect of higher
variance.

``` r
BTM <- BTdecayLasso(NFL$dataframe, NFL$ability, decay.rate = 0.005, fixed = NFL$worstTeam)
```

We can use ‘plot’ function to view the whole Lasso path.

``` r
plot(BTM)
```

The optimal model is selected using AIC criteria on HYBRID Lasso’s run
here.

``` r
BTO <- BTdecayLassoC(NFL$dataframe, NFL$ability, decay.rate = 0.005, fixed = NFL$worstTeam,
                     model = BTM, criteria = "AIC", type = "HYBRID")
summary(BTO)
```

Finally, we use bootstrapping to obtain the standard deviation of this
choosen model with 100 times of simulation.

``` r
BT <- boot.BTdecayLasso(NFL$dataframe, NFL$ability, BTO$Optimal.lambda, decay.rate = 0.005, 
                        fixed = NFL$worstTeam, boot = 100)
summary(BT)
```
