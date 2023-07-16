\name{binomTestCoveragePlot}
\alias{binomTestCoveragePlot}

\title{Coverage Plots for Binomial Proportions}
\description{
  Generates plots for the actual coverage of a binomial proportion
  using various types of confidence intervals. Plots the actual
  coverage for a given sample size and stated nominal coverage \eqn{1 -} \code{alpha}.
}
\usage{
  binomTestCoveragePlot(n,
                        alpha = 0.05,
                        intervalType = "Clopper-Pearson",
                        plo = 0,
                        phi = 1,
                        clo = 1 - 2 * alpha,
                        chi = 1,
                        points = 5 + floor(250 / n),
                        showTrueCoverage = TRUE,
                        gridCurves = FALSE)
}

\arguments{
  \item{n}{sample size}
  \item{alpha}{significance level for confidence interval}
  \item{intervalType}{type of confidence interval used; either "Clopper-Pearson", "Wald", "Wilson-Score", "Jeffreys", "Agresti-Coull", "Arcsine", or "Blaker"}
  \item{plo}{lower limit for percentile (horizontal axis)}
  \item{phi}{upper limit for percentile (horizontal axis)}
  \item{clo}{lower limit for coverage (vertical axis)}
  \item{chi}{upper limit for coverage (vertical axis)}
  \item{points}{number of points plotted in each segment of the plot; if default, varies with 'n' (see above)}
  \item{showTrueCoverage}{logical; if \code{TRUE} (default), a solid red line will appear at \eqn{1 -} \code{alpha}}
  \item{gridCurves}{logical; if \code{TRUE}, display acceptance curves in gray}
  }

\details{
  Generates an actual coverage plot for binomial proportions using
  \itemize{
    \item various types of confidence intervals, and
    \item various sample sizes.
  }

  When the function is called with default arguments,
  \itemize{
    \item the horizontal axis is the percentile at which the coverage is evaluated,
    \item the vertical axis is the actual coverage percentage at each percentile, that is,
    the probability that the true value at a percentile is contained in the corresponding
    confidence interval, and
    \item the solid red line is the stated coverage of \eqn{1 -} \code{alpha}.
  }
  The actual coverage for a particular value of \eqn{p}, the percentile of interest, is
  \deqn{c(p) = \sum_{x=0}^n {I(x,p) {n \choose x} p^x (1-p)^{n-x}},}
  where \eqn{I(x,p)} is an indicator function that determines whether a confidence interval covers \eqn{p} when \eqn{X = x}
  (see Vollset, 1993).

  The binomial distribution with arguments \code{size} = \eqn{n} and
  \code{prob} = \eqn{p} has probability mass function

  \deqn{p(x) = {n \choose x} p^x (1-p)^{n-x}}{%
    p(x) = choose(n, x) p^x (1-p)^(n-x)}

  for \eqn{x = 0, 1, \ldots, n}.

  The algorithm for plotting the actual coverage begins by calculating all possible lower and upper        bounds associated with the confidence interval procedure specified by the \code{intervalType} argument.
  These values are concatenated into a vector which is sorted. Negative values and values that exceed 1 are removed from this vector. These values are the breakpoints in the actual coverage function. The \code{points} argument gives the number of points plotted on each segment of the graph of the actual coverage.

  The \code{plo} and \code{phi} arguments can be used to expand or compress the plots horizontally.

  The \code{clo} and \code{chi} arguments can be used to expand or compress the plots vertically.

  By default, the \code{showTrueCoverage} argument plots a solid horizontal
  red line at the height of the stated coverage. The actual coverage is
  plotted with solid black lines for each segment of the actual coverage.

  The \code{gridCurves} argument is assigned a logical value which indicates whether the acceptance curves giving all possible actual coverage values should be displayed as gray curves.
}

\author{Hayeon Park (\email{hpark031@gmail.com}),
  Larry Leemis (\email{leemis@math.wm.edu})
  }

\references{
  Vollset, S.E. (1993). Confidence Intervals for a Binomial Proportion. Statistics in Medicine, 12, 809--824.
}

\seealso{
  \code{\link{dbinom}}
}

\examples{
  binomTestCoveragePlot(6)
  binomTestCoveragePlot(10, intervalType = "Wilson-Score", clo = 0.8)
  binomTestCoveragePlot(n = 100, intervalType = "Wald", clo = 0, chi = 1, points = 30)
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
