\name{binomTestEnsemble}
\alias{binomTestEnsemble}

\title{Ensemble Confidence Intervals for Binomial Proportions}
\description{
  Generates lower and upper confidence interval limits for a binomial proportion
  using an ensemble of confidence intervals.
}
\usage{
  binomTestEnsemble(n, x,
                    alpha = 0.05,
                    CP = TRUE,
                    WS = TRUE,
                    JF = TRUE,
                    AC = TRUE,
                    AR = TRUE)
}

\arguments{
  \item{n}{sample size}
  \item{x}{number of successes}
  \item{alpha}{significance level for confidence interval}
  \item{CP}{logical; if \code{TRUE} (default), include Clopper-Pearson confidence interval procedure in the ensemble}
  \item{WS}{logical; if \code{TRUE} (default), include Wilson-Score confidence interval procedure in the ensemble}
  \item{JF}{logical; if \code{TRUE} (default), include Jeffreys confidence interval procedure in the ensemble}
  \item{AC}{logical; if \code{TRUE} (default), include Agresti-Coull confidence interval procedure in the ensemble}
  \item{AR}{logical; if \code{TRUE} (default), include Arcsine confidence interval procedure in the ensemble}
}

\details{
  Generates lower and upper confidence interval limits for a binomial proportions using
  \itemize{
    \item various sample sizes,
    \item various numbers of successes, and
    \item various combinations of confidence intervals.
  }

  When the \code{binomTestEnsemble} function is called, it returns a two-element vector in which
  \itemize{
    \item the first element is the lower bound of the Ensemble confidence interval, and
    \item the second element is the upper bound of the Ensemble confidence interval.
  }

 To construct an Ensemble confidence interval that attains an actual coverage that
     is close to the stated coverage, the five constituent confidence interval procedures
     can be combined.  Since these intervals vary in width, the lower limits and the actual
     coverage of the constituent confidence intervals at the maximum likelihood estimator
     are calculated.  Likewise, the upper limits and the actual coverage of the constituent
     confidence intervals at the maximum likelihood estimator are calculated.
     The centroids of the lower and upper constituent confidence intervals for points falling
     below and above the stated coverage are connected with a line segment.  The point of
     intersection of these line segments and the stated coverage gives the lower and upper
     bound of the Ensemble confidence interval.  Special cases to this approach are given in
     the case of (a) the actual coverages all fall above or below the stated coverage, and
     (b) the slope of the line connecting the centroids is infinite.


If only one of the logical arguments is \code{TRUE}, the code returns a simple confidence interval of that one procedure.

The Wald confidence interval is omitted because it degenerates in actual coverage for \eqn{x = 0} and \eqn{x = n}.
}

\references{
Park, H., Leemis, L. (2019), "Ensemble Confidence Intervals for Binomial Proportions", Statistics in Medicine, 38 (18), 3460-3475.
}

\author{Hayeon Park (\email{hpark031@gmail.com}),
  Larry Leemis (\email{leemis@math.wm.edu})
}

\examples{
  binomTestEnsemble(10, 3)
  binomTestEnsemble(100, 82, CP = FALSE, AR = FALSE)
  binomTestEnsemble(33, 1, CP = FALSE, JF = FALSE, AC = FALSE, AR = FALSE)
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
\keyword{Ensemble}
\concept{confidence intervals for binomial proportions}
