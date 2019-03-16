#' @export
###############################################################################
# R function to calculate lower and upper limits of several confidence interval
# procedures for a given sample size n and number of successes x.
###############################################################################
binomTest <- function(n, x,
                      alpha = 0.05,
                      intervalType = "Clopper-Pearson")
{
  #############################################################################

  # first, some parameter checking...

  if (!is.numeric(n) || length(n) != 1 || floor(n) != n  || n <= 0)
    stop("'n' must be a positive integer")
  if (missing(n))
    stop("argument 'n' is missing, with no default")

  if (!is.numeric(x) || length(x) != 1 || floor(x) != x  || x < 0 || x > n)
    stop("'x' must be a positive integer between 0 and n")
  if (missing(x))
    stop("argument 'x' is missing, with no default")

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a numeric value in (0, 1)")

  if (intervalType != "Wald" && intervalType != "Wilson-Score" && intervalType != "Clopper-Pearson"
      && intervalType != "Jeffreys" && intervalType != "Agresti-Coull" && intervalType != "Arcsine"
      && intervalType != "Blaker")
    stop ("'intervalType' is invalid. 'intervalType' must be one of the following: 'Clopper-Pearson',
          'Wald', 'Wilson-Score', 'Jeffreys', 'Agresti-Coull', 'Arcsine', 'Blaker'")

  ##################################################################################################

  if (intervalType == "Clopper-Pearson") {
      l <- qbeta(alpha / 2, x, n - x + 1)
      u <- qbeta(1 - alpha / 2, x + 1, n - x)
      return(c(l, u))
  }

  if (intervalType == "Wald") {
      crit <- qnorm(1 - alpha / 2)
      shat <- x / n
      l <- shat - crit * sqrt(shat * (1 - shat) / n)
      u <- shat + crit * sqrt(shat * (1 - shat) / n)
      if (l < 0) l = 0
      if (u > 1) u = 1
      return(c(l, u))
  }

  if (intervalType == "Wilson-Score") {
      crit <- qnorm(1 - alpha / 2)
      l <- ((x / n) + (crit ^ 2 / (2 * n)) - crit * sqrt(((x / n) * (1 - (x / n))) / n + crit ^ 2 / (4 * n ^ 2))) / (1 + crit ^ 2 / n)
      u <- ((x / n) + (crit ^ 2 / (2 * n)) + crit * sqrt(((x / n) * (1 - (x / n))) / n + crit ^ 2 / (4 * n ^ 2))) / (1 + crit ^ 2 / n)
      if (x == 0) l <- 0
      if (x == n) u <- 1
      return(c(l, u))
  }

  if (intervalType == "Jeffreys") {
    l <- qbeta(alpha / 2, x + 1 / 2, n - x + 1 / 2)
    u <- qbeta(1 - alpha / 2, x + 1 / 2, n - x + 1 / 2)
    if (x == 0) l <- 0
    if (x == n) u <- 1
    return(c(l, u))
  }

  if (intervalType == "Agresti-Coull") {
    crit <- qnorm(1 - alpha / 2)
    nhat <- n + crit ^ 2
    phat <- (x + (1 / 2) * crit ^ 2) / nhat
    l <- phat - crit * sqrt((phat * (1 - phat)) / nhat)
    u <- phat + crit * sqrt((phat * (1 - phat)) / nhat)
    if (l < 0) l <- 0
    if (u > 1) u <- 1
    return(c(l, u))
  }

  if (intervalType == "Arcsine") {
    crit <- qnorm(1 - alpha / 2)
    phat <- (x + 0.375) / (n + 0.75)
    l <- sin(asin(sqrt(phat)) - (crit / (2 * sqrt(n)))) ^ 2
    u <- sin(asin(sqrt(phat)) + (crit / (2 * sqrt(n)))) ^ 2
    if (x == 0) l <- 0
    if (x == n) u <- 1
    return(c(l, u))
  }

  if (intervalType == "Blaker") {
    acceptbin <- function(x,n,p) {
      p1 <- 1 - pbinom(x - 1, n, p)
      p2 <- pbinom(x, n, p)
      a1 <- p1 + pbinom(qbinom(p1, n, p) - 1, n, p)
      a2 <- p2 + 1 - pbinom(qbinom(1 - p2, n, p), n, p)
      return(min(a1, a2))
    }
    acceptinterval <- function(x, n, level, tolerance = 10^(-6)) {
      lower <- 0
      upper <- 1
      if (x != 0) {lower <- qbeta((1 - level) / 2, x, n - x + 1)
        while (acceptbin(x, n, lower) < (1 - level)) lower <- lower + tolerance
      }
      if (x != n) {upper <- qbeta(1 - (1 - level) / 2, x + 1, n - x)
        while (acceptbin(x, n, upper) < (1 - level)) upper <- upper - tolerance
      }
      return(c(lower, upper))
    }
    if (x == 0) {
      l <- 0
      u <- acceptinterval(x, n, 1 - alpha)[2]
    }
    else {
      if (x == n) {
      l <- acceptinterval(x, n, 1 - alpha)[1]
      u <- 1
    }
      else {
      l <- acceptinterval(x, n, 1 - alpha)[1]
      u <- acceptinterval(x, n, 1 - alpha)[2]
      }
    }
    return(c(l, u))
  }
}
