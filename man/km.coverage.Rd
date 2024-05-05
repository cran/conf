\name{km.coverage}
\alias{km.coverage}

\title{Actual Coverage Plots for Several Confidence Interval Procedures for the Survivor Function
for Randomly Right-Censored Data}

\description{
This function calculates and plots the actual coverage functions associated with five confidence
interval procedures associated with the Kaplan--Meier product--limit estimator (KMPLE) for a randomly right-censored data set
with exponential failure times and exponential censoring times.

}

\usage{
 km.coverage(n, lambdaT, lambdaC,
             alpha = 0.1,
             interval = c("Greenwood"),
             show = TRUE, table = FALSE, value = 0)
}

\arguments{
  \item{n}{Sample size; must be between 1 to 15.}
  \item{lambdaT}{Lambda for the exponential failure time distribution.}
  \item{lambdaC}{Lambda for the exponential censoring time distribution.}
  \item{alpha}{Significance level for confidence interval; if equals to 0.1, for example, it indicates a 90 percent confidence interval.}
  \item{interval}{Type of confidence interval used; supported types are "Greenwood", "Exp-Greenwood", "Log-Log", "Arcsine", and "Peto". This parameter can be a vector containing one or more interval types.}
 \item{show}{A logical value indicating whether to display additional potential Probability Mass Function (PMF) lines in gray. When set to \code{TRUE}, the function will plot these gray lines, representing potential curves that the coverage lines could follow. This feature adds a visual layer to the plot, helping to illustrate the range of possible actual coverage values.}
  \item{table}{Logical value. If TRUE, a table of outcomes is printed. The table aids in understanding how the plot is generated.}
  \item{value}{A specific value of p for generating the table.}
}

\details{
The Kaplan--Meier product--limit estimator provides a point estimator for a data set of randomly
right-censored observations.
The \code{km.coverage} function computes the actual coverage function of up to five different
types of confidence interval procedures for the survivor function for various values of its arguments
when failure times are IID exponential random variables with failure rate \code{lambdaT} and
censoring times are IID exponential random variables with censoring rate \code{lambdaC}.
The \code{km.coverage} function
provides options for plotting and tabular display of the outcomes.

  Depending on the interval type, \code{km.coverage} calls internal plotting functions: \code{plot_Greenwood}, \code{plot_ExpGreenwood}, \code{plot_LogLog}, \code{plot_Arcsine}, and \code{plot_Peto}.

  When \code{km.coverage} is called, it generates an actual coverage graph (and a table) based on
  \itemize{
    \item various types of confidence intervals,
    \item various sample sizes,
    \item various failure / censoring distributions (exponential)
    \item if the "table" argument is set to TRUE, it shows the corresponding contribution table at the specified value of p.
  }

The ‘interval’ argument can be one of the character strings \code{"Greenwood"}, \code{"Exp-Greenwood"}, \code{"Log-Log"}, \code{"Arcsine"}, and \code{"Peto"}, or it can be a vector of the confidence interval procedure names, or it can be the character string \code{"all"}, which will plot all five actual coverage functions.  The order in which the actual coverage functions fall in the vector controls the order in which the actual coverage functions are plotted.  The order is significant because it is often the case that one of the actual coverage functions obscures part of another actual coverage function on the plot. For example, the first listed interval in the vector is plotted first. Therefore, it may be covered up by the plots of subsequent intervals in the vector.

This actual coverage plot is constructed by calculating lower and upper bounds associated with
the confidence interval procedure at various permutations of the failure/censoring ordering for a
particular n, the number of objects on test. We calculate the contribution to the actual coverage of
each permutation case based on the probability mass function of the KMPLE. At each value of p
corresponding to S(t), we sum up the contributions which correspond to lower bounds and upper
bounds covering p. These sums give the actual coverage curves.

}

\value{
  This function primarily generates plots illustrating the coverage function for various Kaplan-Meier product-limit confidence intervals. If the \code{table} argument is set to \code{TRUE}, it additionally prints a detailed table of outcomes corresponding to the plot. The table includes calculated coverage values for a specified value of p, providing a complete explanation of how plots are generated. For example, for a value of p along the horizontal axis, the rows of the table show all possible failure/censoring orderings, along with the associated contributions. Adding up the contributions from each row in which the interval covers p gives the point on the plot at p.

  When there is more than one interval, the \code{table} option is disabled.
}

\author{Xingyu Wang (\email{xingyu.wang@yahoo.com}),

  Larry Leemis (\email{lmleem@wm.edu})

  Heather Sasinowska (\email{hdsasinowska@wm.edu})
}

\examples{
  km.coverage(5, 0.5, 1, 0.1, "Greenwood")
  km.coverage(5, 0.5, 1, 0.1, "Greenwood", FALSE)
  km.coverage(5, 0.5, 1, 0.1, "Arcsine", FALSE, TRUE, 0.3)
  km.coverage(5, 0.5, 1, 0.1, "all")
  km.coverage(5, 0.5, 1, 0.1, c("Arcsine", "Greenwood", "Exp-Greenwood"))
}

\keyword{Kaplan Meier Estimator}
\keyword{Right Censoring}
\keyword{Greenwood’s Formula}
\keyword{Survival Analysis}
\keyword{Nonparametric Methods}
