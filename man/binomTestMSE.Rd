\name{binomTestMSE}
\alias{binomTestMSE}

\title{RMSE-Minimizing Confidence Intervals for Binomial Proportions}
\description{
  Generates lower and upper confidence interval limits for a binomial proportion that minimizes
  the root mean square error (RMSE) of the actual coverage function.
}
\usage{
  binomTestMSE(n, x,
               alpha = 0.05,
               smooth = 1,
               showRMSE = TRUE,
               showAll = FALSE)
}

\arguments{
  \item{n}{sample size}
  \item{x}{number of successes}
  \item{alpha}{significance level for confidence interval}
  \item{smooth}{smoothness index}
  \item{showRMSE}{a logical variable indicating whether to show the value of RSME}
  \item{showAll}{a logical variable indicating whether to show confidence intervals of all possible number of successes}
}

\details{
  Generates lower and upper confidence interval limits for a binomial proportion for
  \itemize{
    \item various sample sizes,
    \item various numbers of successes.
  }

  When the \code{binomTestMSE} function is called, it returns a two-element vector in which
  \itemize{
    \item the first element is the lower bound of the RMSE-minimizing confidence interval, and
    \item the second element is the upper bound of the RMSE-minimizing confidence interval.
  }

  An RMSE-minimizing two-sided 100 * (1 - alpha) percent confidence interval
  for p is constructed from a random sample of size n from a Bernoulli(p)
  population.  The parameter \code{x} gives the number of successes in the n
  mutually independent Bernoulli trials.  For n <= 15, all possible jumps
  between acceptance curves associated with the actual coverage function are
  enumerated based on their one-to-one relationship with the symmetric Dyck
  paths.  For each sequence of jumps between acceptance curves, the confidence
  interval bounds that are returned are associated with discontinuities in the
  actual coverage function that together result in the lowest possible RMSE.  A
  set of smoothness constraints that build on four existing non-conservative
  confidence intervals (Wilson-score, Jeffreys, Arcsine, and Agresti-Coull) is
  used if the smoothness index \code{smooth} is set to one.  These constraints
  ensure that the RMSE-confidence interval achieves smoothness, a preferable
  property of the binomial confidence interval that is related to lower bound
  differences for adjacent values of \code{x}.  There is a trade-off between
  the RMSE and the smoothness.  For n > 100, smoothness is required.  The RMSE
  usually increases if the smoothness constraints are used.  For n > 15, only
  the symmetric Dyck paths associated with the Wilson–score, Jeffreys, Arcsine,
  and Agresti–Coull confidence interval procedures are used instead of
  enumerating because the computation time increases in a factorial fashion in
  n.  The minimal RMSE is not guaranteed for n > 15 because another symmetric
  Dyck path other than those associated with the four existing confidence
  interval procedures might prove to be optimal.  However, this procedure does
  ensure a lower RMSE than any of the four existing confidence intervals for
  all n.

}

\references{
  Feng, K., Sasinowska, H., Leemis, L. (2022), "RMSE-Minimizing Confidence Interval
  for the Binomial Parameter", Computational Statistics, 37 (4), 2022, 1855-1885.
}

\author{Kexin Feng (\email{kfeng@caltech.edu}),
  Larry Leemis (\email{leemis@math.wm.edu}),
  Heather Sasinowska (\email{hdsasinowska@wm.edu})
}

\examples{
  binomTestMSE(10, 3)
}

\keyword{actual coverage function}
\keyword{binomial proportion}
\keyword{confidence interval}
\keyword{Dyck words}
\keyword{RMSE}
\concept{confidence intervals for binomial proportions}
