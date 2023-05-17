#' @export
######################################################################
# R function to calculate the support values of the Kaplan-Meier
# product-limit estimator for a given sample size n
######################################################################
km.support <- function(n)
{
######################################################################

  # first, some parameter checking...

  if (!is.numeric(n) || length(n) != 1 || floor(n) != n  || n <= 0)
    stop("'n' must be a positive integer")
  if (missing(n))
    stop("argument 'n' is missing, with no default")

######################################################################

  # next user warnings based on the value of n
  if (n > 35) stop("modify code for n > 35 because of CPU restrictions")
  if (n > 22 && n <= 27) cat("please wait, this might take a few seconds\n")
  if (n > 27) cat("please wait, this might take a few minutes\n")

######################################################################

# memory.limit(size = 160000)
  nn <- 1
  numer <- 1
  denom <- 1
  if (n == 1) return(list(num = c(0, numer), den = c(1, denom)))
  for (nn in 2:n) {
    numer.new <- (nn - 1) * numer
    denom.new <- nn * denom
    for (i in 1:length(numer.new)) {
      numerator <- numer.new[i]
      denominator <- denom.new[i]
      remainder <- -1
      while (remainder != 0) {
        remainder <- numerator %% denominator
        numerator <- denominator
        if (remainder != 0) denominator <- remainder
      }
      numer.new[i] <- numer.new[i] / denominator
      denom.new[i] <- denom.new[i] / denominator
    }
    numer <- c(numer, numer.new)
    denom <- c(denom, denom.new)
    temp  <- complex(real = numer, imaginary = denom)
    temp  <- unique(temp)
    numer <- Re(temp)
    denom <- Im(temp)
  }
  return(list(num = c(0, numer), den = c(1, denom)))
}
