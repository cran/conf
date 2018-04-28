#' @export
###############################################################################
# R function to calculate lower and upper limits of ensemble confidence interval
# procedures for a given sample size n and number of successes x.
###############################################################################
binomTestEnsemble <- function(n, x,
                              alpha = 0.05,
                              CP = TRUE,
                              WS = TRUE,
                              JF = TRUE,
                              AC = TRUE,
                              AR = TRUE)
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

  if (!is.logical(CP) || length(CP) != 1)
    stop("'CP' must be a single logical value")

  if (!is.logical(WS) || length(WS) != 1)
    stop("'WS' must be a single logical value")

  if (!is.logical(JF) || length(JF) != 1)
    stop("'JF' must be a single logical value")

  if (!is.logical(AC) || length(AC) != 1)
    stop("'AC' must be a single logical value")

  if (!is.logical(AR) || length(AR) != 1)
    stop("'AR' must be a single logical value")

  if (sum(c(CP, WS, JF, AC, AR)) == 0)
    stop("at least one of 'CP', 'WS', 'JF', 'AC', 'AR' must be TRUE")

  #############################################################################

  combination <- c(CP, WS, JF, AC, AR)
  len = length(which(combination == 1))
  phat <- x / n

  cp <- binomTest(n, x, alpha = alpha)
  ws <- binomTest(n, x, intervalType = "Wilson-Score", alpha = alpha)
  jf <- binomTest(n, x, intervalType = "Jeffreys", alpha = alpha)
  ac <- binomTest(n, x, intervalType = "Agresti-Coull", alpha = alpha)
  ar <- binomTest(n, x, intervalType = "Arcsine", alpha = alpha)

  accp <- binomTestCoverage(n, phat, alpha = alpha)
  acws <- binomTestCoverage(n, phat, intervalType = "Wilson-Score", alpha = alpha)
  acjf <- binomTestCoverage(n, phat, intervalType = "Jeffreys", alpha = alpha)
  acac <- binomTestCoverage(n, phat, intervalType = "Agresti-Coull", alpha = alpha)
  acar <- binomTestCoverage(n, phat, intervalType = "Arcsine", alpha = alpha)

  actcov <- c(accp,acws,acjf,acac,acar)
  actcov <- actcov[which(combination == 1)]

  CIleft <- c(cp[1], ws[1], jf[1], ac[1], ar[1])
  CIleft <- CIleft[which(combination == 1)]
  xsort <- sort(CIleft)

  indexhi <- which (actcov >= 1 - alpha)
  leftxhi <- mean(CIleft[indexhi])
  leftyhi <- mean(actcov[indexhi])
  indexlo <- which (actcov < 1 - alpha)
  leftxlo <- mean(CIleft[indexlo])
  leftylo <- mean(actcov[indexlo])

  slope <- (leftyhi - leftylo) / (leftxhi - leftxlo)
  eleft <- (1 - alpha + slope * leftxlo - leftylo) / slope

  if (phat == 0) eleft <- 0
  if (length(unique(actcov)) == 1){
    if (mean(actcov) >= 1 - alpha) {
      eleft <- xsort[len]
    }
    if (mean(actcov) < 1 - alpha) {
      eleft <- xsort[1]
    }
  }
  if (min(actcov) >= 1 - alpha) {
    eleft <- xsort[len]
  }
  if (max(actcov) < 1 - alpha) {
    eleft <- xsort[1]
  }
  if (eleft < 0) eleft <- 0
  if (eleft > 1) eleft <- 1

  CIright <- c(cp[2], ws[2], jf[2], ac[2], ar[2])
  CIright <- CIright[which(combination == 1)]
  xrsort <- sort(CIright)
  indexhi <- which (actcov >= 1 - alpha)
  rightxhi <- mean(CIright[indexhi])
  rightyhi <- mean(actcov[indexhi])
  indexlo <- which(actcov < 1 - alpha)
  rightxlo <- mean(CIright[indexlo])
  rightylo <- mean(actcov[indexlo])
  slope <- (rightyhi - rightylo) / (rightxhi - rightxlo)
  eright <- (1 - alpha + slope * rightxlo - rightylo) / slope

  if (phat == 1) eright <- 1
  if (length(unique(actcov)) == 1) {
    if (mean(actcov) >= 1 - alpha) {
      eright <- xrsort[1]
    }
    if (mean(actcov) < 1 - alpha) {
      eright <- xrsort[len]
    }
  }
  if (min(actcov) >= 1 - alpha) {
    eright <- xrsort[1]
  }
  if (max(actcov) < 1 - alpha) {
    eright <- xrsort[len]
  }

  if (eright < 0) eright <- 0
  if (eright > 1) eright <- 1

  return(c(eleft,eright))
}
