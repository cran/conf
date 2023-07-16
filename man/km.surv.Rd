\name{km.surv}
\alias{km.surv}

\title{Probability Mass Functions for the support of the
Kaplan-Meier product-limit estimator for various cumulative probabilities associated with \code{X}}
\description{
 Plot the probability mass functions for the support
 values of the Kaplan-Meier product-limit estimator for
 a given sample size \code{n} with a probability of observing  a failure \code{h} at various times of interest expressed as the cumulative probability \code{perc} associated with \code{X =  min(T, C)}, where \code{T} is the failure time and \code{C} is the censoring time, under a random-censoring scheme.
}
\usage{
km.surv(n, h, lambda, ev, line, graydots, gray.cex,
        gray.outline, xfrac)
}

\arguments{
  \item{n}{sample size}
  \item{h}{probability of observing a failure}
  \item{lambda}{plotting frequency of the probability mass
                functions (default is 10)}
  \item{ev}{option to plot the expected values of the support
            values (default is \code{FALSE})}
  \item{line}{option to connect the expected values with
              lines (default is \code{FALSE})}
  \item{graydots}{option to express the weight of the support
                  values using grayscale (default is
                  \code{FALSE})}
  \item{gray.cex}{option to change the size of the gray dots
                  (default is 1)}
  \item{gray.outline}{option to display outlines of the
                      gray dots (default is \code{TRUE})}
  \item{xfrac}{option to label support values on the y-axis
              as exact fractions (default is \code{TRUE})}
}

\details{
  The Kaplan-Meier product-limit estimator is used to
  estimate the survivor function for a data set of
  positive values in the presence of right censoring.
  The \code{km.surv} function plot the probability mass
  functions for the support values of the Kaplan-Meier
  product-limit estimator for a given sample size \code{n}
  with a probability of observing  a failure \code{h} at
  various times of interest expressed as the cumulative
  probability \code{perc} associated with \code{X =  min(T,
  C)}, where \code{T} is the failure time and \code{C} is the
  censoring time, under a random-censoring scheme.

  The \code{n} argument must be a positive integer denoting
  the sample size. Allowable limits are from 1 to 23.
  Larger values of \code{n} are not allowed because of CPU
  and memory limitations.

  The default method to plot the probability mass functions
  uses the area of a dot to indicate the relative probability
  of a support value. An alternative is to plot the
  probability mass functions using grayscales (by setting
  \code{graydots = TRUE}). One of the two approaches might
  work better in different scenarios.

  The expected values are calculated by removing the
  probability of \code{NA} and normalizing the rest of the
  probabilities.
  }

\value{
  The \code{km.surv} function doesn't return any value.
}

\references{
Qin, Y., Sasinowska, H., Leemis, L. (2023), "The Probability Mass
Function of the Kaplan-Meier Product-Limit Estimator",
\eqn{The American Statistician}, Volume 77, Number 1,
102-110.
}

\author{Yuxin Qin (\email{yqin08@wm.edu}),
  Heather Sasinowska (\email{hdsasinowska@wm.edu}),
  Larry Leemis (\email{leemis@math.wm.edu})}

\seealso{
  \code{\link{survfit}}
}

\examples{
km.surv(n = 4, h = 2/3, lambda = 100, ev = TRUE, line = TRUE)
km.surv(n = 5, h = 3/4, lambda = 50, graydots = TRUE, gray.cex = 0.6, gray.outline = FALSE)
km.surv(n = 7, h = 1/5, lambda = 30, graydots = TRUE, gray.cex = 0.6, xfrac = FALSE)
}

\keyword{Kaplan-Meier estimator}
\keyword{product-limit estimator}
\keyword{right-censored data}
\keyword{right censoring}
\keyword{nonparametric methods}
\concept{Kaplan-Meier product-limit estimator support values}
