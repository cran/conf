#' @export
######################################################################
# R function to generate a matrix that contains all possible delta
# combinations with n items on test, as well as the associated
# support values as numeric values and fractions
######################################################################
km.outcomes <- function(n)
{
######################################################################
  
  # first, some parameter checking...
  
  if (!is.numeric(n) || length(n) != 1 || floor(n) != n  || n <= 0)
    stop("'n' must be a positive integer")
  if (missing(n))
    stop("argument 'n' is missing, with no default")
  
######################################################################
  
  # next user warnings based on the value of n

  if (n > 24) stop("modify code for n > 24 because of CPU restrictions")
  if (n > 16 && n <= 21) cat("please wait, this might take a few seconds\n")
  if (n > 21) cat("please wait, this might take a few minutes\n")
  
######################################################################
  
  # memory.limit(size = 160000) 
  a <- matrix(c(rep(-1, n), rep(1, 3)), 1, n + 3)
  for (i in 1:n) {
    b <- as.matrix(expand.grid(replicate(i, 0:1, simplify = FALSE)))
    c <- matrix(NA, nrow(b), n + 3 - i)
    b <- cbind(b, c)
    a <- rbind(a, b)
  }
  numerator <- (n - 1):0
  denominator <- n:1
  for (k in 2:nrow(a)) {
    if (is.na(a[k, n]) || a[k, n] == 1) {
      i <- which(!is.na(a[k, 1:n]) & a[k, 1:n] == 1)
      numer <- prod(numerator[i])
      denom <- prod(denominator[i])
      a[k, (n + 1):(n + 3)] <- c(numer / denom, numer, denom) 
    }
  }
  loc <- c(0, rep(1:n, 2 ^ (1:n)))
  a <- cbind(loc, a)
  colnames(a) <- c("l", paste("d", c(1:n), sep = ""), "S(t)", "num", "den")
  return(a)
}
