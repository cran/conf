#' @export
###############################################################################
# R function to plot the actual coverage of several confidence interval
# procedures for a given sample size n.
###############################################################################
binomTestCoveragePlot <- function(n,
                              alpha = 0.05,
                              intervalType = "Clopper-Pearson",
                              plo = 0,
                              phi = 1,
                              clo = 1 - 2 * alpha,
                              chi = 1,
                              points = 5 + floor(250 / n),
                              showTrueCoverage = TRUE,
                              gridCurves = FALSE)
{
  #############################################################################

  # first, some parameter checking...

  if (!is.numeric(n) || length(n) != 1 || floor(n) != n  || n <= 0)
    stop("'n' must be a positive integer")
  if (missing(n))
    stop("argument 'n' is missing, with no default")

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a numeric value in (0, 1)")

  if (intervalType != "Wald" && intervalType != "Wilson-Score" && intervalType != "Clopper-Pearson"
      && intervalType != "Jeffreys" && intervalType != "Agresti-Coull" && intervalType != "Arcsine"
      && intervalType != "Blaker")
    stop ("'intervalType' is invalid. 'intervalType' must be one of the following: 'Clopper-Pearson',
          'Wald', 'Wilson-Score', 'Jeffreys', 'Agresti-Coull', 'Arcsine', 'Blaker'")

  if (!is.numeric(plo) || length(plo) != 1)
    stop("'plo' must be a numeric value")

  if (!is.numeric(phi) || length(phi) != 1)
    stop("'phi' must be a numeric value")

  if (plo >= phi)
    stop("'plo' must be less than 'phi'")

  if (!is.numeric(clo) || length(clo) != 1)
    stop("'clo' must be a numeric value")

  if (!is.numeric(chi) || length(chi) != 1)
    stop("'chi' must be a numeric value")

  if (clo >= chi)
    stop("'clo' must be less than 'chi'")

  if (!is.numeric(points) || length(points) != 1 || floor(points) != points || points <= 1)
    stop("'points' must be a positive integer greater than 1")

  if (!is.logical(showTrueCoverage) || length(showTrueCoverage) != 1)
    stop("'showTrueCoverage' must be a single logical value")

  if (!is.logical(gridCurves) || length(gridCurves) != 1)
    stop("'gridCurves' must be a single logical value")

 ##################################################################################################

  plot(NULL, xlab = "p", ylab = "coverage", xlim = c(plo, phi), ylim = c(clo, chi),
       main = paste(intervalType, "coverage for n = ", n), las = 1)
  if (showTrueCoverage) abline(h = 1 - alpha, col = "red")
  if (gridCurves) {
    p = seq(0, 1, by = 0.001)
    for (x in 0:n) {
      for (y in x:n) {
        coverage = rep(0, length(p))
        for (z in x:y) coverage = coverage + dbinom(z, n, p)
        lines(p, coverage, col = "gray", lwd = 0.4)
      }
    }
  }

  if (intervalType == "Clopper-Pearson") {
    l <- rep(0, n + 1)
    u <- rep(1, n + 1)
    for (x in 1:n)
      l[x + 1] <- 1 / (1 + (n - x + 1) / (x * qf(alpha / 2, 2 * x, 2 * (n - x + 1))))
    for (x in 0:(n - 1))
      u[x + 1] <- 1 / (1 + (n - x) / ((x + 1) * qf(1 - alpha / 2, 2 * (x + 1), 2 * (n - x))))
    b <- sort(c(l, u))
    for (i in 1:(length(b) - 1)) {
      p <- seq(b[i], b[i + 1], length = points)
      coverage <- rep(0, points)
      mid <- (b[i] + b[i + 1]) / 2
      for (x in 0:n)
        if (l[x + 1] < mid && mid < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
      lines(p, coverage, type = "l")
    }
  }

  if (intervalType == "Wald") {
    l <- rep(0, n + 1)
    u <- rep(0, n + 1)
    crit <- qnorm(1 - alpha / 2)
    for (x in 0:n) {
      shat <- x / n
      l[x + 1] <- shat - crit * sqrt(shat * (1 - shat) / n)
      if (l[x + 1] < 0) l[x + 1] <- 0
      u[x + 1] <- shat + crit * sqrt(shat * (1 - shat) / n)
      if (u[x + 1] > 1) u[x + 1] <- 1
    }
    b <- sort(c(l[l >= 0], u[u <= 1]))
    b <- unique(b)
    for (i in 1:(length(b) - 1)) {
      p <- seq(b[i], b[i + 1], length = points)
      coverage <- rep(0, points)
      mid <- (b[i] + b[i + 1]) / 2
      for (x in 0:n)
        if (l[x + 1] < mid && mid < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
      lines(p, coverage, type = "l")
    }
  }

  if (intervalType == "Wilson-Score") {
    l <- rep(0, n + 1)
    u <- rep(1, n + 1)
    crit <- qnorm(1 - alpha / 2)
    for (x in 1:n)
      l[x + 1] <- ((x / n) + (crit ^ 2 / (2 * n)) - crit * sqrt(((x / n) * (1 - (x / n))) / n + crit ^ 2 / (4 * n ^ 2))) / (1 + crit ^ 2 / n)
    for (x in 0:(n - 1))
      u[x + 1] <- ((x / n) + (crit ^ 2 / (2 * n)) + crit * sqrt(((x / n) * (1 - (x / n))) / n + crit ^ 2 / (4 * n ^ 2))) / (1 + crit ^ 2 / n)
    b <- sort(c(l, u))
    for (i in 1:(length(b) - 1)) {
      p <- seq(b[i], b[i + 1], length = points)
      coverage <- rep(0, points)
      mid <- (b[i] + b[i + 1]) / 2
      for (x in 0:n)
        if (l[x + 1] < mid && mid < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
      lines(p, coverage, type = "l")
    }
  }

  if (intervalType == "Jeffreys") {
    l <- rep(0, n + 1)
    u <- rep(1, n + 1)
    for (x in 1:n)
      l[x + 1] <- qbeta(alpha / 2, x + 1 / 2, n - x + 1 / 2)
    for (x in 0:(n - 1))
      u[x + 1] <- qbeta(1 - alpha / 2, x + 1 / 2, n - x + 1 / 2)
    b <- sort(c(l, u))
    for (i in 1:(length(b) - 1)) {
      p <- seq(b[i], b[i + 1], length = points)
      coverage <- rep(0, points)
      mid <- (b[i] + b[i + 1]) / 2
      for (x in 0:n)
        if (l[x + 1] < mid && mid < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
      lines(p, coverage, type = "l")
    }
  }

  if (intervalType == "Agresti-Coull") {
    l <- rep(0, n + 1)
    u <- rep(0, n + 1)
    crit <- qnorm(1 - alpha / 2)
    nhat <- n + crit ^ 2
    for (x in 0:n) {
      phat <- (x + (1 / 2) * crit ^ 2) / nhat
      l[x + 1] <- phat - crit * sqrt((phat * (1 - phat)) / nhat)
      u[x + 1] <- phat + crit * sqrt((phat * (1 - phat)) / nhat)
    }
    l[l < 0] <- 0
    u[u > 1] <- 1
    b <- sort(c(l, u))
    b <- unique(b)
    for (i in 1:(length(b) - 1)) {
      p <- seq(b[i], b[i + 1], length = points)
      coverage <- rep(0, points)
      mid <- (b[i] + b[i + 1]) / 2
      for (x in 0:n)
         if (l[x + 1] < mid && mid < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
      lines(p, coverage, type = "l")
    }
  }

  if (intervalType == "Arcsine") {
    l <- rep(0, n + 1)
    u <- rep(1, n + 1)
    crit <- qnorm(1 - alpha / 2)
    for (x in 1:n) {
      phat <- (x + 0.375) / (n + 0.75)
      l[x + 1] <- sin(asin(sqrt(phat)) - (crit /  (2 * sqrt(n)))) ^ 2
    }
    for (x in 0:(n - 1)){
      phat <- (x + 0.375) / (n + 0.75)
      u[x + 1] <- sin(asin(sqrt(phat)) + (crit /  (2 * sqrt(n)))) ^ 2
    }
    b <- sort(c(l[l >= 0], u[u <= 1]))
    b <- unique(b)
    for (i in 1:(length(b) - 1)) {
      p <- seq(b[i], b[i + 1], length = points)
      coverage <- rep(0, points)
      mid <- (b[i] + b[i + 1]) / 2
      for (x in 0:n)
        if (l[x + 1] < mid && mid < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
      lines(p, coverage, type = "l")
    }
  }

  if (intervalType == "Blaker") {
    l <- rep(0, n + 1)
    u <- rep(1, n + 1)
    for (x in 1:n) {
      l[x + 1] <- binomTest(n, x, alpha = alpha, intervalType = "Blaker")[1]
    }
    for (x in 1:n) {
      u[x] <- binomTest(n, x - 1, alpha = alpha, intervalType = "Blaker")[2]
    }
    b <- sort(c(l, u))
    for (i in 1:(length(b) - 1)) {
      p <- seq(b[i], b[i + 1], length = points)
      coverage <- rep(0, points)
      mid <- (b[i] + b[i + 1]) / 2
      for (x in 0:n)
        if (l[x + 1] < mid && mid < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
      lines(p, coverage, type = "l")
    }
  }
}
