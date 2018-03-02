## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.width = 4, fig.height = 4, fig.show = 'hold'-------------------
library(conf)
ballbearing <- c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.48, 51.84,
                 51.96, 54.12, 55.56, 67.80, 68.64, 68.64, 68.88, 84.12,
                 93.12, 98.64, 105.12, 105.84, 127.92, 128.04, 173.40)
crplot(dataset = ballbearing, alpha = 0.05, distn = "weibull") 

## ---- fig.width = 4, fig.height = 4, fig.show = 'hold'-------------------
crplot(dataset = ballbearing, alpha = 0.05, distn = "weibull", 
       pts = FALSE, sf = c(2, 4), origin = TRUE)

## ---- fig.height = 4, fig.width = 4, fig.keep = 'none'-------------------
x <- crplot(dataset = ballbearing, alpha = 0.05, distn = "weibull", 
       pts = FALSE, sf = c(2, 4), origin = TRUE, info = TRUE)
names(x)

## ---- fig.height = 4, fig.width = 3.5, fig.show = 'asis'-----------------
# with confidence region data stored in x, it is now available for custom graphics
plot(x$kappa, x$lambda, type = 'l', lty = 5, xlim = c(0, 3), ylim = c(0, 0.0163),
     xlab = expression(kappa), ylab = expression(lambda))
polygon(x$kappa, x$lambda, col = "gray80", border = NA)

## ---- fig.show = 'hold'--------------------------------------------------
# record MLE values (previously output to screen when info = TRUE) & reproduce CR plot
kappa.hat <- 2.10206
lambda.hat <- 0.01221
crplot(dataset = ballbearing, alpha = 0.05, distn = "weibull", xlab = expression(kappa), 
       ylab = expression(lambda), main = paste0("Confidence Region"),
       pts = FALSE, mlelab = FALSE, sf = c(2, 4), origin = TRUE)

# 1st analysis plot: overlay phi angles (as line segments) on the confidence region plot
par(xpd = FALSE)
segments(rep(kappa.hat, length(x$phi)), rep(lambda.hat, length(x$phi)),
         x1 = kappa.hat * cos(x$phi) + kappa.hat, y1 = kappa.hat * sin(x$phi) + lambda.hat, lwd = 0.2)

# 2nd analysis plot: CDF of phi angles reveals most are very near 0 (2pi) and pi
plot.ecdf(x$phi, pch = 20, cex = 0.1, axes = FALSE, xlab = expression(phi), ylab = "", main = "CDF")
axis(side = 1, at = round(c(0, pi, 2*pi), 2))
axis(side = 2, at = c(0, 1), las = 2)

## ---- fig.show = 'hold'--------------------------------------------------
crplot(dataset = ballbearing, alpha = 0.05, distn = "weibull", 
       maxdeg = 3, main = "maxdeg = 3")
crplot(dataset = ballbearing, alpha = 0.05, distn = "weibull", 
       maxdeg = 3, main = "maxdeg = 3 (pts hidden)", pts = FALSE)
crplot(dataset = ballbearing, alpha = 0.05, distn = "weibull", 
       maxdeg = 20, main = "maxdeg = 20")
crplot(dataset = ballbearing, alpha = 0.05, distn = "weibull", 
       maxdeg = 20, main = "maxdeg = 20 (pts hidden)", pts = FALSE)

## ---- fig.show = 'hold'--------------------------------------------------
crplot(dataset = ballbearing, alpha = 0.05, distn = "invgauss", sf = c(2, 2),
       main = "default; heuristic = 1")
crplot(dataset = ballbearing, alpha = 0.05, distn = "invgauss", sf = c(2, 2), 
       heuristic = 0, ellipse_n = 100, main = "heuristic = 0")

## ---- fig.show = 'hold'--------------------------------------------------
crplot(dataset = ballbearing, alpha = 0.05, distn = "invgauss", sf = c(2, 2),
       heuristic = 0, ellipse_n = 40, main = "combination: step 1")
crplot(dataset = ballbearing, alpha = 0.05, distn = "invgauss", sf = c(2, 2), 
       maxdeg = 10, ellipse_n = 40, main = "combination: step 2")

## ---- fig.show = 'hold'--------------------------------------------------
crplot(dataset = ballbearing, alpha = 0.05, distn = "invgauss", sf = c(2, 2), 
       maxdeg = 10, main = "default (heuristic = 1)")

