\name{km.support}
\alias{km.support}

\title{Support values for the Kaplan-Meier product-limit estimator}
\description{
  Calculate the support values for the Kaplan-Meier product-limit
  estimator for a particular sample size \code{n} using an induction
  algorithm.
}
\usage{
km.support(n)
}

\arguments{
  \item{n}{sample size}
}

\details{
  The Kaplan-Meier product-limit estimator is used to
  estimate the survivor function for a data set of
  positive values in the presence of right censoring.
  The \code{km.support} function calculates the support values for the
  Kaplan-Meier product-limit estimator for a sample of
  size \code{n} using an induction algorithm
  described in Qin et al. (2023).

  The \code{n} argument must be a positive integer denoting
  the sample size. Allowable limits are from 1 to 35.
  Larger values of \code{n} are not allowed because of CPU
  and memory limitations.

  The numerators and denominators are temporarily converted to
  complex numbers within the \code{km.support} function in order to
  eliminate duplicate support values using the \code{unique} function.
  }

\value{
  The \code{km.support} function returns a list with two components.
\itemize{
\item \code{num}: a vector of integers containing the numerators of the
support values
\item \code{den}: a vector of integers containing the associated
denominators of the support values
}
The support values are not returned in sorted order.
Zero and one, which are always a part of the support,
are given as 0 / 1 and 1 / 1.
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
  #  display unsorted numerators and denominators of support values for n = 4
  km.support(4)
  #  display sorted support values for n = 4 as exact fractions
  n <- 4
  s <- km.support(n)
  i <- order(s$num / s$den)
  m <- length(s$num)
  f <- ""
  for (j in i[2:(m - 1)]) f <- paste(f, s$num[j], "/", s$den[j], ", ", sep = "")
  cat(paste("The ", m, " support values for n = ", n, " are: 0, ", f, "1.\n", sep = ""))
  #  print sorted support values for n = 4 as numerics
  print(s$num[i] / s$den[i])
}

\keyword{Kaplan-Meier estimator}
\keyword{product-limit estimator}
\keyword{right-censored data}
\keyword{right censoring}
\keyword{nonparametric methods}
\concept{Kaplan-Meier product-limit estimator support values}
