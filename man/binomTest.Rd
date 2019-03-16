\name{binomTest}
\alias{binomTest}

\title{Confidence Intervals for Binomial Proportions}
\description{
  Generates lower and upper confidence interval limits for a binomial proportion
  using different types of confidence intervals.
}
\usage{
 binomTest(n, x,
           alpha = 0.05,
           intervalType = "Clopper-Pearson")
}

\arguments{
  \item{n}{sample size}
  \item{x}{number of successes}
  \item{alpha}{significance level for confidence interval}
  \item{intervalType}{type of confidence interval used; either "Clopper-Pearson", "Wald", "Wilson-Score", "Jeffreys", "Agresti-Coull", "Arcsine", or "Blaker"}
}

\details{
  Generates a lower and upper confidence interval limit for a binomial proportion using
  \itemize{
    \item various types of confidence intervals,
    \item various sample sizes, and
    \item various numbers of successes.
  }

  When the \code{binomTest} function is called, it returns a two-element vector in which
  \itemize{
    \item the first element is the lower bound of the confidence interval, and
    \item the second element is the upper bound of the confidence interval.
  }

  This confidence interval is constructed by calculating lower and upper bounds associated with the confidence interval procedure specified by the \code{intervalType} argument. Lower bounds that are negative are set to 0 and upper bounds that are greater than 1 are set to 1.
}

\author{Hayeon Park (\email{hpark03@email.wm.edu}),
  Larry Leemis (\email{leemis@math.wm.edu})}

\seealso{
  \code{\link{dbinom}}
}

\examples{
  binomTest(10, 6)
  binomTest(100, 30, intervalType = "Agresti-Coull")
}

\keyword{confidence interval}
\keyword{binomial proportion}
\keyword{binomial distribution}
\keyword{Clopper-Pearson}
\keyword{Wilson-Score}
\keyword{Wald}
\keyword{Agresti-Coull}
\keyword{Jeffreys}
\keyword{Arcsine}
\keyword{Blaker}
\concept{confidence intervals for binomial proportions}
