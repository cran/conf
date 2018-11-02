\name{binomTestCoverage}
\alias{binomTestCoverage}

\title{Actual Coverage Calculation for Binomial Proportions}
\description{
  Calculates the actual coverage of a confidence interval for a
  binomial proportion for a particular sample size \eqn{n} and
  a particular value of the probability of success \eqn{p} for several
  confidence interval procedures.
}
\usage{
  binomTestCoverage(n, p,
                    alpha = 0.05,
                    intervalType = "Clopper-Pearson")
}

\arguments{
  \item{n}{sample size}
  \item{p}{population probability of success}
  \item{alpha}{significance level for confidence interval}
  \item{intervalType}{type of confidence interval used; either "Clopper-Pearson", "Wald", "Wilson-Score", "Jeffreys", "Agresti-Coull", "Arcsine", or "Blaker"}
  }

\details{
  Calculates the actual coverage of a confidence interval
  procedure at a particular value of \eqn{p} for
  \itemize{
    \item various types of confidence intervals,
    \item various probabilities of success \eqn{p}, and
    \item various sample sizes \eqn{n}.
  }

  The actual coverage for a particular value of \eqn{p}, the probability of success of interest, is
  \deqn{actual coverage(p) = \sum_{x=0}^n {I(x,p) {n \choose x} p^x (1-p)^{(n-x)}},}
  where \eqn{I(x)} is an indicator function that determines if a confidence interval covers \eqn{p} when \eqn{X = x}
  which is from Vollset(1993).

  The binomial distribution with arguments \code{size} = \eqn{n} and
  \code{prob} = \eqn{p} has pmf

  \deqn{p(x) = {n \choose x} p^x (1-p)^{(n-x)}}{%
    p(x) = choose(n, x) p^x (1-p)^(n-x)}

  for \eqn{x = 0, 1, \ldots, n}.

  The algorithm for computing the actual coverage for a particular probability of
  success begins by calculating all possible lower and upper bounds associated
  with the confidence interval procedure specified by the \code{intervalType} argument.
  The appropriate binomial probabilities are summed to determine the actual coverage
  at \eqn{p}.
}

\examples{
  binomTestCoverage(6, 0.4)

  binomTestCoverage(n = 10, p = 0.3, alpha = 0.01, intervalType = "Wilson-Score")

}

\seealso{
  \code{\link{dbinom}}
}

\references{
  Vollset, S.E. (1993). Confidence Intervals for a Binomial Proportion. Statistics in Medicine, 12, 809-824.
}

\author{Hayeon Park (\email{hpark03@email.wm.edu}),
  Larry Leemis (\email{leemis@math.wm.edu})
  }

\keyword{actual coverage}
\keyword{confidence interval}
\keyword{binomial proportion}
\keyword{Clopper-Pearson}
\keyword{Wilson-Score}
\keyword{Wald}
\keyword{Agresti-Coull}
\keyword{Jeffreys}
\keyword{Arcsine}
\keyword{Blaker}
\concept{confidence intervals for binomial proportions}
