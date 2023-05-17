\name{km.outcomes}
\alias{km.outcomes}

\title{
Outcomes for the Kaplan-Meier product-limit estimator
}
\description{
Generates a matrix containing all possible outcomes (all possible sequences of failure times and right-censoring times) of the value
of the Kaplan-Meier product-limit estimator for a particular sample
size \code{n}.
}
\usage{
km.outcomes(n)
}

\arguments{
  \item{n}{sample size}
}

\details{
  The Kaplan-Meier product-limit estimator is used to
  estimate the survivor function for a data set of
  positive values in the presence of right censoring. 
  The \code{km.outcomes} function generates a matrix with
  all possible combinations of observed failures and
  right censored values and the resulting support values
  for the Kaplan-Meier product-limit estimator for a sample of
  size \code{n}.
  
  The \code{n} argument must be a positive integer denoting
  the sample size. Allowable limits are from 1 to 24.
  Larger values of \code{n} are not allowed because of CPU
  and memory limitations. 
  
  In order to keep the support values as exact fractions,
  the numerators and denominators are stored separately in
  the \code{a} matrix in the columns named \code{num} and
  \code{den}. The support values are stored as numeric
  values in the column named \code{S(t)}.
}

\value{
  The \code{km.outcomes} function returns a matrix with
  \ifelse{html}{\out{2<sup>n+1</sup>-1}}{\eqn{2 ^ {n + 1} - 1}} rows and \ifelse{html}{\out{n + 4}}{\eqn{n + 4}} columns. The location \code{l} indicates the position where the time of interest falls within the observed events.
  The meaning of the columns is as follows.
  \itemize{
    \item \code{l}: number of observed events (failures times or
    censoring times) between times 0 and the observation time;
    \item \code{d1, d2, ..., dn}: equals 0 if the event corresponds
          to a censored observation, equals 1 if the event
          corresponds to a failure;
    \item \code{S(t)}: numeric value of the associated support value; 
    \item \code{num}: numerator of the support value as a fraction;
    \item \code{den}: denominator of the support value as a fraction.
  }
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
km.outcomes(3)
}

\keyword{Kaplan-Meier estimator}
\keyword{product-limit estimator}
\keyword{right-censored data}
\keyword{right censoring}
\keyword{nonparametric methods}
\concept{Kaplan-Meier product-limit estimator support values}
