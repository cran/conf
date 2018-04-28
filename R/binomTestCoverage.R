#' @export
###############################################################################
# R function to calculate the actual coverage of several confidence
# interval procedures for a given sample size n and binomial parameter p.
###############################################################################
binomTestCoverage <- function(n, p,
                              alpha = 0.05,
                              intervalType = "Clopper-Pearson")
{
  #############################################################################

  # first, some parameter checking...

  if (!is.numeric(n) || length(n) != 1 || floor(n) != n  || n <= 0)
    stop("'n' must be a positive integer")
  if (missing(n))
    stop("argument 'n' is missing, with no default")

  if (!is.numeric(p) || length(p) != 1 || p < 0 || p > 1)
    stop("'p' must be a numeric value in [0, 1]")
  if (missing(p))
    stop("argument 'p' is missing, with no default")

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a numeric value in (0, 1)")

  if (intervalType != "Wald" && intervalType != "Wilson-Score" && intervalType != "Clopper-Pearson"
      && intervalType != "Jeffreys" && intervalType != "Agresti-Coull" && intervalType != "Arcsine"
      && intervalType != "Blaker")
    stop ("'intervalType' is invalid. 'intervalType' must be one of the following: 'Clopper-Pearson',
          'Wald', 'Wilson-Score', 'Jeffreys', 'Agresti-Coull', 'Arcsine', 'Blaker'")


  ##################################################################################################

  if (intervalType == "Clopper-Pearson") {
    l <- rep(0, n + 1)
    u <- rep(1, n + 1)
    for (x in 1:n)
      l[x + 1] <- 1 / (1 + (n - x + 1) / (x * qf(alpha / 2, 2 * x, 2 * (n - x + 1))))
    for (x in 0:(n - 1))
      u[x + 1] <- 1 / (1 + (n - x) / ((x + 1) * qf(1 - alpha / 2, 2 * (x + 1), 2 * (n - x))))
    coverage <- 0
    for (x in 0:n) {
      if (l[x + 1] < p && p < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
    }
    if (p == 0 || p == 1) coverage <- 1
    return(coverage)
  }

  if (intervalType == "Wald") {
    l <- rep(0, n + 1)
    u <- rep(0, n + 1)
    crit <- qnorm(1 - alpha / 2)
    for (x in 0:n) {
      shat <- x / n
      l[x + 1] <- shat - crit * sqrt(shat * (1 - shat) / n)
      u[x + 1] <- shat + crit * sqrt(shat * (1 - shat) / n)
    }
    coverage <- 0
    for (x in 0:n) {
      if (l[x + 1] < p && p < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
    }
    if (p == 0 || p == 1) coverage <- 0
    return(coverage)
  }

  if (intervalType == "Wilson-Score") {
    l <- rep(0, n + 1)
    u <- rep(1, n + 1)
    crit <- qnorm(1 - alpha / 2)
    for (x in 1:n)
      l[x + 1] <- ((x / n) + (crit ^ 2 / (2 * n)) - crit * sqrt(((x / n) * (1 - (x / n))) / n + crit ^ 2 / (4 * n ^ 2))) / (1 + crit ^ 2 / n)
    for (x in 0:(n - 1))
      u[x + 1] <- ((x / n) + (crit ^ 2 / (2 * n)) + crit * sqrt(((x / n) * (1 - (x / n))) / n + crit ^ 2 / (4 * n ^ 2))) / (1 + crit ^ 2 / n)
    coverage <- 0
    for (x in 0:n) {
      if (l[x + 1] < p && p < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
    }
    if (p == 0 || p == 1) coverage <- 1
    return(coverage)
  }

  if (intervalType == "Jeffreys") {
    l <- rep(0, n + 1)
    u <- rep(1, n + 1)
    crit <- qnorm(1 - alpha / 2)
    for (x in 1:n)
      l[x + 1] <- qbeta(alpha / 2, x + 1/2, n - x + 1/2)
    for (x in 0:(n - 1))
      u[x + 1] <- qbeta(1 - alpha / 2, x + 1/2, n - x + 1/2)
    coverage <- 0
    for (x in 0:n) {
      if (l[x + 1] < p && p < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
    }
    if (p == 0 || p == 1) coverage <- 1
    return(coverage)
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
    coverage <- 0
    for (x in 0:n) {
      if (l[x + 1] < p && p < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
    }
    if (p == 0 || p == 1) coverage <- 1
    return(coverage)
  }

  if (intervalType == "Arcsine") {
    l <- rep(0, n + 1)
    u <- rep(0, n + 1)
    crit <- qnorm(1 - alpha / 2)
    for (x in 0:n) {
      phat <- (x + 0.375) / (n + 0.75)
      l[x + 1] <- sin(asin(sqrt(phat)) - (crit /  (2 * sqrt(n)))) ^ 2
      u[x + 1] <- sin(asin(sqrt(phat)) + (crit /  (2 * sqrt(n)))) ^ 2
    }
    coverage <- 0
    for (x in 0:n) {
      if (l[x + 1] < p && p < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
    }
    if (p == 0 || p == 1) coverage <- 0
    return(coverage)
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
    coverage <- 0
    for (x in 0:n) {
      if (l[x + 1] < p && p < u[x + 1]) coverage <- coverage + dbinom(x, n, p)
    }
    if (p == 0 || p == 1) coverage <- 1
    return(coverage)
  }
}
