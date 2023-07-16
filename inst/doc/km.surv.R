## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
library(conf)
#  display the probability mass functions at various times of interest

km.surv(n = 3, h = 1/2)


## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
#  display the probability mass function at various times of interest
#  with the expected values in red connected with lines

km.surv(n = 4, h = 1/3, lambda = 100, ev = TRUE, line = TRUE)
title("High Censoring Rate")


## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
#  display the probability mass function at various times of interest
#  with the expected values in red connected with lines

km.surv(n = 4, h = 2/3, lambda = 100, ev = TRUE, line = TRUE)
title("High Failure Rate")

## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
#  display the probability mass function at various times of interest
km.surv(n = 7, h = 3/4, lambda = 50, graydots = TRUE, xfrac = FALSE)

## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
#  display the probability mass function at various times of interest
km.surv(n = 7, h = 3/4, lambda = 50, graydots = TRUE, xfrac = FALSE, gray.outline = FALSE)

## ---- fig.width = 6, fig.height = 5, fig.show = 'hold'------------------------
#  display the probability mass function at various times of interest
km.surv(n = 5, h = 5/8, lambda = 30, graydots = TRUE, ev = TRUE, gray.outline = FALSE, gray.cex = 1.25)

