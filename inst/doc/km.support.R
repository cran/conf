## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.width = 4, fig.height = 4, fig.show = 'hold'------------------------
library(conf)
#  display unsorted numerators and denominators of support values for n = 4
n = 4
s = km.support(n)
s
#  display sorted support values for n = 4 as decimals
sort(s$num / s$den)
#  display sorted support values for n = 4 as exact fractions
i <- order(s$num / s$den)
m <- length(s$num)
f <- ""
for (j in i[2:(m - 1)]) f <- paste(f, s$num[j], "/", s$den[j], ", ", sep = "")
cat(paste("The ", m, " support values for n = ", n, " are: 0, ", f, "1.\n", sep = ""))


