\name{km.pmf}
\alias{km.pmf}

\title{
Probability Mass Function for the support of the Kaplan-Meier product-limit estimator
}
\description{
Generates the probability mass function for the support values
of the Kaplan-Meier product-limit estimator for a particular sample
size \code{n}, probability of observing a failure \code{h} at the time of interest expressed as the cumulative probability \code{perc} associated with \code{X = min(T, C)}, where \code{T} is the failure time and \code{C} is the censoring time under a random-censoring scheme.
}
\usage{
km.pmf(n, h, perc, plot, sep, xfrac, cex.lollipop)
}

\arguments{
  \item{n}{sample size}
  \item{h}{probability of observing a failure, in other words, \code{P(X = T)}}
  \item{perc}{cumulative probability associated with \code{X = min(T, C)}}
  \item{plot}{option to plot the probability mass function (default is \code{TRUE})}
  \item{sep}{option to show the breakdown of the probability for each support value (see function \code{km.outcomes} for details on the breakdown) (default is \code{TRUE})}
  \item{xfrac}{option to label support values on the x-axis as exact fractions (default is \code{TRUE})}
  \item{cex.lollipop}{size of the dots atop the spikes}
}

\details{
  The Kaplan-Meier product-limit estimator is used to
  estimate the survivor function for a data set of
  positive values in the presence of right censoring. 
  The \code{km.pmf} function generates the probability mass function for the support values
of the Kaplan-Meier product-limit estimator for a particular sample
size \code{n}, probability of observing a failure \code{h} at the time of interest expressed as the cumulative probability \code{perc} associated with \code{X = min(T, C)}, where \code{T} is the failure time and \code{C} is the censoring time under a random-censoring scheme.
  
  The \code{n} argument must be a positive integer denoting
  the sample size. Allowable limits are from 1 to 23.
  Larger values of \code{n} are not allowed because of CPU
  and memory limitations. 
  
  For larger sample size \code{n}, it is recommended to set 
  \code{sep = FALSE}, \code{xfrac = FALSE}, and 
  \code{cex.lollipop = 0.01} for a better visual effect.
}

\value{
  The \code{km.pmf} function returns a dataframe with
  2 columns. The column named \code{S} stores all the support 
  values for the Kaplan-Meier product-limit estimator 
  with sample size \code{n}, including \code{NA}. The 
  column named \code{P} stores the associated probabilities.
}

\references{
Qin, Y., Sasinowska, H., Leemis, L. (2023), "The Probability Mass
Function of the Kaplan-Meier Product-Limit Estimator",
\eqn{The American Statistician}, Volume 77, Number 1, 102-110.
}

\author{Yuxin Qin (\email{yqin08@wm.edu}),
  Heather Sasinowska (\email{hdsasinowska@wm.edu}),
  Larry Leemis (\email{leemis@math.wm.edu})
}

\seealso{
  \code{\link{survfit}}
}

\examples{
km.pmf(4, 1/3, 0.75)

km.pmf(8, 1/2, 0.75, sep = FALSE, xfrac = FALSE, cex.lollipop = 0.01)
}

\keyword{Kaplan-Meier estimator}
\keyword{product-limit estimator}
\keyword{right-censored data}
\keyword{right censoring}
\keyword{nonparametric methods}
\concept{Kaplan-Meier product-limit estimator support values}
