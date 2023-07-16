## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
library(conf)
#  display the probability mass function

km.pmf(n=3, h = 1/2, perc = 0.5)


## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
#  display the probability mass function with lower failure rate (higher censoring rate)

km.pmf(n=3, h = 1/3, perc = 0.75, sep = FALSE, xfrac = FALSE, cex.lollipop = 2)


## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
#  display the probability mass function with an earlier time of interest

km.pmf(n=3, h = 1/3, perc = 0.35, sep = FALSE, xfrac = FALSE, cex.lollipop = 2)


## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
#  display the probability mass function with a higher probability of failure (lower censoring rate)

km.pmf(n=3, h = 2/3, perc = 0.75, sep = FALSE, xfrac = FALSE, cex.lollipop = 2)


## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
#  display the probability mass function with a sample size of 8
km.pmf(8, 1/2, 0.75, sep = FALSE, xfrac = FALSE, cex.lollipop = 0.01)

