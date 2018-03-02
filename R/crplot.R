#' Plotting Two-Dimensional Confidence Regions
#'
#' @description
#' Plots the two-dimensional confidence region for probability distribution (Weibull or
#' inverse Gaussian) parameters corresponding to a user given dataset and level of significance.
#'
#' @param dataset a 1 x n vector of dataset values.
#' @param alpha significance level; resulting plot illustrates a 100(1 - alpha)\% confidence region.
#' @param distn distribution to fit the dataset to; accepted values: \code{"weibull"}, \code{"invgauss"}.
#' @param heuristic numeric value selecting method for plotting: 0 for elliptic-oriented point distribution, and
#' 1 for smoothing boundary search heuristic.
#' @param maxdeg maximum angle tolerance between consecutive plot segments in degrees.
#' @param ellipse_n number of roughly equidistant confidence region points to plot using the
#' elliptic-oriented point distribution (must be a multiple of four because its algorithm
#' exploits symmetry in the quadrants of an ellipse).
#' @param pts displays confidence region boundary points identified if \code{TRUE}.
#' @param mlelab logical argument to include the maximum
#' likelihood estimate coordinate point (default is \code{TRUE}).
#' @param sf specifies the number of significant figures on axes labels.
#' @param mar specifies margin values for \code{par(mar = c( ))} (see \code{mar} in \code{\link{par}}).
#' @param xlab string specifying the x axis label.
#' @param ylab string specifying the y axis label.
#' @param main string specifying the plot title.
#' @param xlas numeric in {0, 1, 2, 3} specifying the style of axis labels (see \code{las} in \code{\link{par}}).
#' @param ylas numeric in {0, 1, 2, 3} specifying the style of axis labels (see \code{las} in \code{\link{par}}).
#' @param origin logical argument to include the plot origin (default is \code{FALSE}).
#' @param xlim two element vector containing horizontal axis minimum and maximum values.
#' @param ylim two element vector containing vertical axis minimum and maximum values.
#' @param tol the \code{\link{uniroot}} parameter specifying its required accuracy.
#' @param info logical argument to return plot information: MLE prints to screen; (x, y) plot point coordinates
#' and corresponding phi angles (with respect to MLE) are returned as a list.
#' @param showplot logical argument specifying if a plot is output; altering from its default of \code{TRUE} is
#' only logical assuming \code{crplot} is run for its data only (see the \code{info} argument).
#' @import stats
#' @import graphics
#' @export
#' @return if the optional argument \code{info = TRUE} is included then a list of plot coordinates and phi angles is returned
# @aliases
#' @keywords Graphical Methods, Parameter Estimation, Numerical Optimization
# @examples
#' @seealso \code{\link{uniroot}}
#' @author Christopher Weld (\email{ceweld@email.wm.edu})
#' @author Lawrence Leemis (\email{leemis@math.wm.edu})
#' @keywords confidence region, confidence intervals, statistical graphics, data visualization
#'
#' @usage
#' crplot(dataset, alpha, distn,
#'                 heuristic = 1,
#'                 maxdeg    = 5,
#'                 ellipse_n = 4,
#'                 pts       = TRUE,
#'                 mlelab    = TRUE,
#'                 sf        = c(5, 5),
#'                 mar       = c(4, 4.5, 2, 1.5),
#'                 xlab      = "",
#'                 ylab      = "",
#'                 main      = "",
#'                 xlas      = 1,
#'                 ylas      = 2,
#'                 origin    = FALSE,
#'                 xlim      = NULL,
#'                 ylim      = NULL,
#'                 tol       = .Machine$double.eps ^ 0.5,
#'                 info      = FALSE,
#'                 showplot  = TRUE )
#'
#' @details
#' This function supports two-dimensional confidence region plots for Weibull or inverse Gaussian parameters.
#' The first input argument (shape for Weibull, mean for inverse Gaussian) is given on its horizontal axis,
#' and the second (scale for Weibull, shape for inverse Gaussian) on the vertical axis.  It requires
#' \itemize{
#' \item a vector of dataset values,
#' \item the level of significance (alpha), and
#' \item a distribution (Weibull or inverse Gaussian) to fit the data to.
#' }
#' Two heuristics (and their associated combination) are available to plot confidence regions.  Along
#' with their descriptions, they are:
#' \enumerate{
#' \item \emph{Smoothing Boundary Search Heuristic (default)}.  This heuristic plots more points in areas of
#' greater curvature to ensure a smooth appearance throughout the confidence region boundary.  Its
#' \code{maxdeg} parameter specifies the maximum tolerable angle between three successive points.
#' Lower values of \code{maxdeg} result in smoother plots, and its default value of 5 degrees
#' provides adequate smoothing in most circumstances.  Values of \code{maxdeg} \eqn{\le} 3 are not
#' permitted due to their complicating implications to trigonometric numerical approximations near 0
#' and 1.
#' \item \emph{Elliptic-Oriented Point Distribution}.  This heuristic allows the user to specify
#' a number of points to plot along the confidence region boundary at roughly uniform intervals.
#' Its name is derived from the technique it uses to choose these points---an extension of the Steiner
#' generation of a non-degenerate conic section, also known as the parallelogram method---which identifies
#' points along an ellipse that are approximately equidistant.  To exploit the computational benefits of
#' ellipse symmetry over its four quadrants, \code{ellipse_n} value must be divisible by four.}
#' By default, \code{crplot} implements the smoothing boundary search heuristic.  Alternatively,
#' the user can plot using the elliptic-oriented point distribution algorithm, or a combination
#' of them both.  Combining the two techniques initializes the plot using the elliptic-oriented point
#' distribution algorithm, and then subsequently populates additional points in areas of high curvature
#' (those outside of the maximum angle tolerance parameterization) in accordance with the smoothing
#' boundary search heuristic.  This combination results when the smoothing boundary search heuristic
#' is specified in conjunction with an \code{ellipse_n} value greater than four.
#'
#' Both of the aforementioned heuristics use a radial profile log likelihood function to identify
#' points along the confidence region boundary.  It cuts the log likelihood function in a directional
#' azimuth from its MLE, and locates the associated confidence region boundary point using the
#' asymptotic results associated with the ratio test statistic \eqn{-2 [log L(\theta) - log L(\theta hat)]}
#' which converges in distribution to the chi-square distribution with two degrees of freedom (for
#' a two parameter distribution).
#'
#' @examples
#' ## plot the 95% confidence region for Weibull shape and scale parameters
#' ## corresponding to the given ballbearing dataset
#' ballbearing <- c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.48, 51.84,
#'                  51.96, 54.12, 55.56, 67.80, 68.64, 68.64, 68.88, 84.12,
#'                  93.12, 98.64, 105.12, 105.84, 127.92, 128.04, 173.40)
#' crplot(dataset = ballbearing, distn = "weibull", alpha = 0.05)
#'
#' ## repeat this plot using the elliptic-oriented point distribution heuristic
#' crplot(dataset = ballbearing, distn = "weibull", alpha = 0.05,
#'        heuristic = 0, ellipse_n = 80)
#'
#' ## combine the two heuristics, compensating any elliptic-oriented point verticies whose apparent
#' ## angles > 6 degrees with additional points, and expand the plot area to include the origin
#' crplot(dataset = ballbearing, distn = "weibull", alpha = 0.05,
#'        maxdeg = 6, ellipse_n = 80, origin = TRUE)
#'
#' ## next use the inverse Gaussian distribution and show no plot points
#' crplot(dataset = ballbearing, distn = "invgauss", alpha = 0.05,
#'        pts = FALSE)


#######################################################################################
# This R function plots confidence regions.
#
# Name             : crplot.R
# Authors          : Chris Weld & Larry Leemis
# Language         : R (part of conf package)
# Latest Revision  : Feb 2018
#######################################################################################
crplot <- function(dataset,
                   alpha,
                   distn,
                   heuristic = 1,
                   maxdeg = 5,
                   ellipse_n = 4,
                   pts = TRUE,
                   mlelab = TRUE,
                   sf = c(5, 5),
                   mar = c(4, 4.5, 2, 1.5),
                   xlab = "",
                   ylab = "",
                   main = "",
                   xlas = 1,
                   ylas = 2,
                   origin = FALSE,
                   xlim = NULL,
                   ylim = NULL,
                   tol = .Machine$double.eps ^ 0.5,
                   info = FALSE,
                   showplot = TRUE) {

  # parameter error checking ###########################################################

  if (missing(dataset)) stop ("argument 'dataset' is missing, with no default")
  if (missing(alpha)) stop ("argument 'alpha' is missing, with no default")
  if (missing(distn)) stop ("argument 'distn' is missing, with no default")

  if (is.null(dataset) || length(dataset) == 1 || !is.numeric(dataset))
    stop("'dataset' must be a numeric vector with length > 1")

  if (distn != "weibull" && distn != "invgauss")
    stop("'distn' invalid; only 'weibull' and 'invgauss' are supported")

  if (distn == "weibull" && min(dataset) <= 0)
    stop("'dataset' parameter contains infeasible Weibull distribution outcome(s) <= 0")

  if (distn == "invgauss" && min(dataset) <= 0)
    stop("'dataset' parameter contains infeasible inverse Gaussian distribution outcome(s) <= 0")

  if (is.null(alpha) || alpha <= 0 || alpha >= 1 || !is.numeric(alpha) || length(alpha) != 1)
    stop("'alpha' numeric significance level parameter required such that 0 < alpha < 1")

  if (is.null(heuristic) || !is.numeric(heuristic) || (heuristic != 0 && heuristic != 1))
    stop("'heuristic' parameter must be 0 (elliptic-oriented points) or 1 (search heuristic)")

  if (is.null(maxdeg) || !is.numeric(maxdeg) || maxdeg < 3 || maxdeg >= 90)
    stop("'maxdeg' parameter must be a numeric angle tollerance in degrees such that 3 <= maxdeg < 90")

  if ((ellipse_n == 4) && (heuristic == 0))
    stop("'ellipse_n' (number of plot points) >= 8 required for 'heuristic = 0' (elliptic-oriented points)")

  if (!is.numeric(ellipse_n) || length(ellipse_n) != 1 || (ellipse_n) %% 4 != 0 || ellipse_n <= 0 )
    stop("'ellipse_n' must be a positive numeric multiple of 4")

  if (!is.logical(pts) || length(pts) != 1)
    stop("'pts' must be a single logical parameter")

  if (!is.logical(mlelab) || length(mlelab) != 1)
    stop("'mlelab' must be a single logical parameter")

  if (length(sf) != 2 || !is.numeric(sf) || floor(sf)[1] != sf[1] || floor(sf)[2] != sf[2] || min(sf) <= 0)
    stop("'sf' must be a numeric vector with positive integer values of length two")

  if (length(mar) != 4 || !is.numeric(mar) || min(mar) < 0)
    stop("'mar' must be a vector of length four with positive numeric entries")

  if (!(xlas %in% c(1, 2, 3, 4)))
    stop("'xlas' must be a numeric value in {0, 1, 2, 3}")

  if (!(ylas %in% c(1, 2, 3, 4)))
    stop("'ylas' must be a numeric value in {0, 1, 2, 3}")

  if (!is.logical(origin) || length(origin) != 1)
    stop("'origin' must be a single logical parameter")

  if (origin == TRUE && !(is.null(xlim)) && xlim[1] > 0)
    stop("given 'origin' = TRUE, xlim---when specified---must include the origin")

  if (origin == TRUE && !(is.null(ylim)) && ylim[1] > 0)
    stop("given 'origin' = TRUE, ylim---when specified---must include the origin")

  if (!(is.null(xlim)) && (length(xlim) != 2 || !is.numeric(xlim) || xlim[1] >= xlim[2]) )
    stop("'xlim' must be a numeric vector of length two with xlim[1] < xlim[2]")

  if (!(is.null(ylim)) && (length(ylim) != 2 || !is.numeric(ylim) || ylim[1] >= ylim[2]) )
    stop("'ylim' must be a numeric vector of length two with ylim[1] < ylim[2]")

  if (is.null(tol) || !is.numeric(tol) || length(tol) != 1)
    stop("'tol' numeric parameter given is invalid (default .Machine$double.eps^0.3)")

  if (!is.logical(info) || length(info) != 1)
    stop("'info' must be a single logical parameter")

  if (!is.logical(showplot) || length(showplot) != 1)
    stop("'showplot' must be a single logical parameter")

  if (showplot == FALSE && info == FALSE)
    stop("'showplot' and 'info' are both FALSE; without either requirement crplot will not execute")

  ######################################################################################

  # sequence of code given below:
  #
  # First the following functions---all contained within crplot---are given:
  #
  # 1. steinerframe.   This function takes in a, b, and # points and returns associated (x, y) points around a
  #                    rectangle that enable graphic portrayal of an ellipse via the parallelogram method.
  #                    Matrix returned is (n x 4), with columns 1 & 2 giving x-y coordinates of points
  #                    from V1 (R side) vertex, and 3 & 4 corresponding to points from V2 (L side) vertex.
  #
  # 2. angles.         This function takes in a, b, and # points and returns the angles associated with
  #                    graphic presentation via the parallelogram method, also using the "points" fn.
  #                    This function also calls the steinerframe function.
  #
  # 3. mlesolve.       This function calculates and returns the MLE.  Weibull MLEs use Braxton Fixed Point
  #                    Algorithm with initial Menon estimate to begin iterating the algorithm (ref: pg 338
  #                    in Leemis Reliability text)
  #
  # 4. llrsolve.       Using the asymptotically chisquared likelihood ratio statistic, this function returns
  #                    a function whose positive values reside within the (1 - alpha)% confidence region,
  #                    and negative values reside outside of it.
  #
  # 5. crsolve.        This function, leveraging the llrsolve function, calls uniroot to determine the
  #                    cartesian coordinates of confidence region boundary points corresponding to the input
  #                    vector of phi angles (given with respect to the MLE).
  #
  # 6. crsmooth.       Confidence region smoothing function.  Given an initial vector of angles to
  #                    calculate confidence region points with, it iteratively adds additional angles
  #                    until all apparent angle bends are within the specified angle tolerance.
  #
  # Utilizing the above listed functions, a short block of code at the bottom of this file then creates
  # the confidence region plot.

  ######################################################################################

  # "steinerframe" function
  # - three inputs: relative x and y axis lengths and # points (a, b, & n, respectively)
  # - returns matrix of (x, y) coordinate pairs along the perimeter of the rectangle defined by a and b corresponding to
  # those necessary to identify points on an ellipse via the parallelogram method.
  # - points given assume its center is at the origin
  # - define the n input parameter as the number of segments in 2b (or equivelantly in a/2)
  steinerframe = function(a, b, n){
    space13 <- (2 * a) / n                                             # sides 1 and 3 horizontal spacing
    space24 <- (2 * b) / n                                             # sides 2 and 4 verticle spacing
    # v1 gives (x, y) to-coordinate of segments extending from (a, 0)
    v1x <- rev(seq(-a, a - space13, by = space13))                     # side 3, R to L
    v1x <- c(v1x, rep(-a, 2 * n - 2))                                  # side 2, T to B
    v1x <- c(v1x, seq(-a, a - space13, by = space13))                  # side 1, L to R
    v1y <- rep(2 * b, n)                                 # side 3, R to L
    v1y <- c(v1y, rev(seq(space24, 2 * b - space24, by = space24)))    # side 2, T to middle
    v1y <- c(v1y, rev(seq(-2 * b + space24, -space24, by = space24)))  # side 2, middle to B
    v1y <- c(v1y, rep(-2 * b, n))                        # side 1, L to R
    # v2 gives (x, y) to-coordinate of segments extending from (-a, 0)
    v2x <- rep(a, n - 1)                                               # side 4, middle to T
    v2x <- c(v2x, rev(seq(-a + space13, a, by = space13)))             # side 3, R to L
    v2x <- c(v2x, seq(-a + space13, a, by = space13))                  # side 1, L to R
    v2x <- c(v2x, rep(a, n - 1))                                       # side 4, B to middle
    v2y <- seq(space24, 2 * b - space24, by = space24)                 # side 4, middle to T
    v2y <- c(v2y, rep(2 * b, n))                                       # side 3, R to L
    v2y <- c(v2y, rep(-2 * b, n))                                      # side 1, L to R
    v2y <- c(v2y, seq(-2 * b + space24, -space24, by = space24))       # side 4, B to middle
    # assemble and return (x, y) coordinate pairs that enable generation of points along the
    # perimeter of an eclipse via the parallelogram method
    v1 <- cbind(v1x, v1y)
    v2 <- cbind(v2x, v2y)
    v <- cbind(v1, v2)
    return(v)
  }

  ######################################################################################

  # "angles" function
  # - three inputs: relative x and y axis lengths and # points (a, b, & n, respectively)
  # - invokes the parrallelogram method by leveraging the steinerframe function (above) to produce
  # roughly equally-spaced points along the circumference of an ellipse
  # - returns a vector of angles (in radians) associated with those roughly equally-spaced points on the ellipse
  angles = function(a, b, npoints) {

    # isolate and store coordinate pairs along rectangle boundary using steinerframe function
    ntop <- npoints / 4
    ntot <- 6 * ntop
    v12 <- steinerframe(a, b, ntop)
    v1 <- v12[, 1:2]
    v2 <- v12[, 3:4]
    v1fixed <- cbind(c(rep(a, length(v1) / 2)), c(rep(0, length(v1) / 2)))
    v2fixed <- -v1fixed

    # identify intercept point lying on the ellipse boudary (ex, ey) using m1 and m2 slopes from v1 and v2, respectively
    m1 = (v1fixed[ , 2] - v1[ , 2]) / (v1fixed[ , 1] - v1[ , 1])
    m2 = (v2fixed[ , 2] - v2[ , 2]) / (v2fixed[ , 1] - v2[ , 1])
    ex = -(m1 + m2) * a / (m1 - m2)
    ey = m1 * ex + m1 * a

    # determine phi values corresponding to (ex, ey) coordinate pairs
    theta = atan(ey / ex)                                             # gives (-pi / 2) < theta < (pi / 2)
    theta360 = theta[1:(length(theta) / 2)]                           # eliminate redundancies
    theta360 = c(theta360, theta360 + pi, 0.000001, pi)               # mirror on horizontal axis and augment w/ ~0 and 180 degrees
    theta360neg = sign(sign(sign(theta360) * (1 - sign(theta360))))   # vector with "-1" entry where negative theta360 entries occur
    theta360 = theta360 - 2 * pi * theta360neg                        # convert neg angle entries to pos equivelant angles
    theta360 = sort(theta360)
    n = length(theta)

    return(theta360)
  }

  ######################################################################################

  # mlesolve -------------------------------------------------------------------
  # This function calculates and returns the maximum likelihood estimator for the specified distribution
  # as well as the value of its loglikelihood function at the MLE.  Weibull MLEs are found using the
  # Braxton Fixed Point Algorithm, and inverse Gaussian MLEs are found using its closed form solution.
  # Returned is mle.list with values list("theta1.hat", "theta2.hat", "mleLLvalue").
  mlesolve = function(x, cen, epsilon = 0.000000001){
    n <- length(x)
    r <- sum(cen)

    # Weibull MLE
    if (distn == "weibull"){
      # c0 is an educated guess for the initial Menon estimate to begin iterating through the algorithm
      # reference: Appendix B (pg 338) in Leemis Reliability text
      temp1 <- sum(log(x) ^ 2)
      temp2 <- (sum(log(x)) ^ 2)
      temp3 <- 6 / ((n - 1) * pi ^ 2)
      temp4 <- (temp3 * (temp1 - temp2 / n))
      c0 <- 1 / sqrt(temp4)
      s1 <- 0
      for(i in 1:n) {
        if (cen[i] == 1) s1 = s1 + log(x[i])
      }
      repeat{
        s2 <- sum(x ^ c0)
        s3 <- sum(x ^ c0 * log(x))
        q <- (r * s2) / (r * s3 - s1 * s2)
        c0 <- (c0 + q) / 2
        if(abs(c0 - q) < epsilon){
          break
        }
      }
      kappa.hat <- c0
      lambda.hat <- 1 / (sum(x ^ kappa.hat) / r) ^ (1 / kappa.hat)
      # compute log likelihood function for censored weibull
      # reference: pg 246 in Leemis Reliability text
      temp1 <- sum(x ^ kappa.hat)
      temp2 <- 0
      for(i in 1:n){
        if(cen[i] == 1){
          temp2 <- temp2 + log(x[i])
        }
      }
      mleLL <- (r * log(kappa.hat) + (kappa.hat * r * log(lambda.hat)) + ((kappa.hat - 1) * temp2) -
                     (lambda.hat ^ kappa.hat * temp1))
      mle.list <- list("theta1.hat" = kappa.hat, "theta2.hat" = lambda.hat, "mleLLvalue" = mleLL)
    }

    # inverse Gaussian MLE
    else if (distn == "invgauss"){
      mu.hat <- sum(x) / length(x)
      tempsum <- 0
      for (i in 1:length(x)) {
        tempsum <- tempsum + (1 / x[i])
      }
      lambda.hat <-  ((1 / (length(x))) * tempsum -
                        (length(x) / sum(x))) ^ -1
      # determine loglikelihood value at MLE
      temp1 <- 0
      temp2 <- 0
      for(i in 1:n){
        if(cen[i] == 1){
          temp1 <- temp1 + log(x[i])
          temp2 <- temp2 + ((x[i] - mu.hat) ^ 2) / x[i]
        }
      }
      mleLL <- (n / 2) * log(lambda.hat) - (n / 2) * log(2 * pi) - (3 / 2) * temp1 -
        (lambda.hat / (2 * mu.hat ^ 2)) * temp2
      mle.list <- list("theta1.hat" = mu.hat, "theta2.hat" = lambda.hat, "mleLLvalue" = mleLL)
    }

    invisible(mle.list)
  }

  ######################################################################################

  # llrsolve (loglikelihood ratio equation solve supporting uniroot) ----------------------------
  # using the fact that the likelihood ratio statistic is asymptotically chi-squared(2) for the two
  # parameter estimation case (reference page 248 Leemis Reliability text):
  # 2 * [loglikelihood(theta1.hat, theta2.hat) - loglikelihood(theta1, theta2)] ~ chisquared(2)
  # this function returns the corresponding value of logliklihood(theta1, theta2) for the given
  # significance level.  It returns a function whose positive values reside within the (1 - alpha)%
  # confidence region, and negative values reside outside of it.
  # An example line of code implimenting the uniroot function by using llrsolve is:
  # g <- uniroot(llrsolve, lower = 0, upper = tempUpper, phi = phi[i], MLE = MLEHAT, x = samp,
  #             cen = cen, chi2 = qchisq(1-alpha, 2))
  llrsolve = function(d, phi, MLEHAT, mleLLvalue, x, cen, chi2){
    n <- length(x)
    r <- sum(cen)

    # Weibull distribution
    if (distn == "weibull") {
      kappa.hat <- MLEHAT[1, ]
      lambda.hat <- MLEHAT[2, ]
      #if (bartlett == FALSE) {
      temp <- mleLLvalue - (chi2 / 2)
      #}
      #else {
      #  c <- 1 / (1 + ((1 / n) * 1.73482) / 2)
      #  temp <- mleLLvalue - (chi2 / (2 * c))
      #}
      kappa <- kappa.hat + (d * cos(phi))
      lambda <- lambda.hat + (d * sin(phi))
      n <- length(x)
      r <- sum(cen)
      temp1 <- sum(x ^ kappa)
      temp2 <- 0
      for(i in 1:n){
        if(cen[i] == 1){
          temp2 <- temp2 + log(x[i])
        }
      }
      valuereturn <- (r * log(kappa) + (kappa * r * log(lambda)) + ((kappa - 1) * temp2) - (lambda ^ kappa * temp1) - temp)
    }

    # inverse Gaussian distribution
    if(distn == "invgauss"){
      mu.hat <- MLEHAT[1, ]
      lambda.hat <- MLEHAT[2, ]

      # plug mleLL and chisquare values into loglikelihood ratio eqn
      temp <- mleLLvalue - (chi2 / 2)
      mu <- mu.hat + (d * cos(phi))
      lambda <- lambda.hat + (d * sin(phi))
      temp1 <- 0
      temp2 <- 0
      for(i in 1:n){
        if(cen[i] == 1){
          temp1 <- temp1 + log(x[i])
          temp2 <- temp2 + ((x[i] - mu) ^ 2) / x[i]
        }
      }

      # return value for invgauss distn (positive when inside conf region, neg outside)
      valuereturn <- ((n / 2) * log(lambda) - (n / 2) * log(2 * pi) - (3 / 2) * temp1 -
                  (lambda / (2 * mu ^ 2)) * temp2 - temp)
    }

    invisible(valuereturn)
  }

  ######################################################################################

  # crsolve -------------------------------------------------------------------
  # leveraging the llrsolve---which returns a function containing positive values within the (1 - alpha)%
  # confidence region, and negative values outside it---this function calls uniroot to determine the cartesian
  # coordinates of points along this boundary corresponding to the input vector of phi angles (given in
  # radians from an MLE perspective with 0 corresponding to the positive horizontal axis)
  crsolve = function(samp, cen, alpha = alpha, mle.list = mle.list, phi){
    n <- length(samp)
    r <- sum(cen)
    theta1.hat <- mle.list$theta1.hat
    theta2.hat <- mle.list$theta2.hat
    mleLLvalue <- mle.list$mleLLvalue
    MLEHAT <- matrix(c(theta1.hat, theta2.hat), ncol =1)
    npoints <- length(phi)
    d_phi <- numeric(npoints)                               # will hold the radial length d for each phi
    cartesian <- matrix(0, nrow = length(phi), ncol = 2)    # This holds cartesian coordinates for each phi
    # finds radial distance and catesian coord of boundary point for a given angle
    for (j in 1:npoints){
      done <- 0
      for (umult in c(10, 100, 500)) {   # uniroot 'upper' argument will use this multiplier to set >> upper bounds in search of root
        if (done == 0) {
          if (phi[j] <= pi / 2) {                                                # phi angles of 0 to pi / 2
            tempUpper <- umult * max(theta1.hat, theta2.hat)
          }
          else if (phi[j] <= pi +  atan(theta2.hat/theta1.hat)) {                 # phi angles of pi / 2 through intercept with origin
            temp <- theta1.hat / (-cos(phi[j]))    # an upper bound for confidence region feasible points; -cos() because cos(phi[i]) < 0
            tempUpper <- temp * 0.99               # set search limit just before upper limit (y-axis)
          }
          else {
            temp <- theta2.hat/ (-sin(phi[j]))   # an upper bound for confidence region feasible points; -sin() because sin(phi[j]) < 0
            tempUpper <- temp * 0.99             # set search limit just before x-axis
          }
          if ((tempUpper > umult * theta1.hat) && (tempUpper > umult * theta2.hat)) {
            tempUpper <- umult * max(c(theta1.hat, theta2.hat))     # arbitrary upper bound for phi near pi/2; accept risk CR bound <= umult * max(mle_parameter)
          }
          g <- suppressWarnings(try(g1 <- uniroot(llrsolve, lower = 0, upper = tempUpper,
                                                  phi = phi[j], MLE = MLEHAT, mleLLvalue = mleLLvalue, x = samp,
                                                  cen = cen, chi2 = qchisq(1 - alpha, 2), tol = tol)))

          # errors in uniroot (calculating "g") are possible, and a result of poor upper bound selection
          # the next section of code troubleshoots those errors in two ways:
          # 1. it decreases the upper value by a factor of 1/2 (per iteration), and also
          # 2. it identifies a relatively small upper value (based on min(theta1.hat, theta2.hat)) and then increments it
          # The first option is valuable when the upper limit couldn't be evaluated (-infinity)
          # The second option is valuable when vastly different parameter orders of magnitude push the search too far for phi calcs ~parallel the smaller axis
          z <- 0                       # counter for "try-error" recovery attempts below
          tempUpperLo <- 0             # initialize; will incrementally push upward seeking uniroot upper param
          tempUpperHi <- tempUpper     # initialize; will incrementally push upward seeking uniroot upper param
          while ((class(g) == "try-error") && (z < 6) && (tempUpperLo < tempUpperHi)) {
            #print(paste0("-------------- problematic phi value: ", phi[j]))
            z <- z + 1
            tempUpperHi <- tempUpper * (0.5 ^ z)
            tempUpperLo <- 5 ^ (z - 1) * min(c(theta1.hat, theta2.hat))
            #print(paste0(z, " ****************************"))
            #print(paste0("Error; correction is sought with uniroot upper bound modifications... ", tempUpperHi))
            g <- suppressWarnings(try(g1 <- uniroot(llrsolve, lower = 0, upper = tempUpperHi,
                                     phi = phi[j], MLE = MLEHAT, mleLLvalue = mleLLvalue, x = samp,
                                     cen = cen, chi2 = qchisq(1 - alpha, 2), tol = tol)))
            if (class(g) == "try-error") {
              #print(paste0("...............try pushing up the min value to ", tempUpperLo))
              g <- suppressWarnings(try(g1 <- uniroot(llrsolve, lower = 0, upper = tempUpperLo,
                           phi = phi[j], MLE = MLEHAT, mleLLvalue = mleLLvalue, x = samp,
                           cen = cen, chi2 = qchisq(1 - alpha, 2), tol = tol)))
            }
            if (class(g) != "try-error") {
              print("...*error overcome* through uniroot 'upper' argument modification; algorithm resuming...")
            }
          }
          if (class(g) != "try-error") {
            done <- 1
          }
        }   # end if (done == 0)
      }     # end for (umult in ...)
      if (class(g) == "try-error") {
        print("------------------------------------------------------------------------------------------")
        print("Uniroot failure.  Challenging parameters and/or shape requires customizing uniroot bounds.")
        print("------------------------------------------------------------------------------------------")
      }
      d_phi[j] <- g$root                                                    # saves radial length per phi
      cartesian[j,] <- MLEHAT + d_phi[j] * (c(cos(phi[j]), sin(phi[j])))    # saves CR coordinates per phi
    }       # end for (j in 1:npoints)
    invisible(cartesian)     # return confidence region coordinates w/o displaying to screen
  }         # end crsolve

  ######################################################################################

  crsmooth = function(maxdeg, samp, cen, alpha, mle.list, ellipse_n)  {
    #rm(cr)                               # remove previously stored values
    maxrad <- maxdeg * (pi / 180)         # allowable angle converted to radians
    maxcurrent <- pi                      # initialize as > maxrad to enter while loop
    theta1.hat <-  mle.list$theta1.hat
    theta2.hat <- mle.list$theta2.hat
    if(ellipse_n == 4){
      phi <- seq(0, 2*pi, length.out = 5)   # initialize phi search with length.out equidistant points
      phi <- phi[2:length(phi)]             # eliminate redundant point
    }
    else {
      phi <- angles(theta1.hat, theta2.hat, ellipse_n)
    }
    transf.scale <- 1                     # initialize transformation scaling factor at 1
    count <- 0                            # administrative counter for while loop
    d01 <- 0
    d02 <- 0
    d12 <- 0
    cr <- crsolve(samp, cen, alpha = alpha, mle.list = mle.list, phi)      # confidence region points w initial phi's
    while (maxcurrent > maxrad) {                                          # in radians
      count <- count + 1
      if (count != 1){
        index_phinew <- match(phinew, phi)
        index_phiold <- match(phiold, phi)
        crnew <- crsolve(samp, cen, alpha = alpha, mle.list = mle.list, phinew)     # uniroot calcs for CR points of new phis
        crold <- cr
        cr <- matrix(rep_len(0, length(phi) * 2), length(phi), 2)          # initialize to store CR coordinate
        cr[index_phiold, ] <- crold
        cr[index_phinew, ] <- crnew
      }
      phiold <- phi                        # store phi values for next iteration
      xspan <- max(cr[,1]) - min(cr[,1])
      yspan <- max(cr[,2]) - min(cr[,2])
      transf.scale <- xspan / yspan
      cr.transf <- cr
      cr.transf[, 2] <- transf.scale * cr[, 2]

      # this sub-section plots the progression of added points, showing interim steps:
      #plot(cr, main = "in-progress build of confidence region", cex = 0.7, pch = 4, col = 'gray30',
      #     axes = FALSE, xlab = xlab, ylab = ylab)
      #lines(cr, col = 'yellow3', lty = 4)
      #segments(cr[dim(cr)[1], 1], cr[dim(cr)[1], 2], cr[1, 1], cr[1, 2], col = 'yellow3', lty = 4)
      #points(theta1.hat, theta2.hat, pch = 8, cex = 0.8, col = 'gray30')
      #axis(side = 1)
      #axis(side = 2)

      from <- cr.transf
      to1 <- rbind(cr.transf[nrow(cr.transf), ], cr.transf[1:(nrow(cr.transf) - 1), ])  # offset confidence region to previous point
      to2 <- rbind(cr.transf[2:nrow(cr.transf), ], cr.transf[1, ])                      # offset confidence region to next point
      d01 <- sqrt((cr.transf[, 1] - to1[, 1]) ^ 2 + (cr.transf[, 2] - to1[, 2]) ^ 2)    # distance to previous point
      d02 <- sqrt((cr.transf[, 1] - to2[, 1]) ^ 2 + (cr.transf[, 2] - to2[, 2]) ^ 2)    # distance to next point
      d12 <- sqrt((to1[, 1] - to2[, 1]) ^ 2 + (to1[, 2] - to2[, 2]) ^ 2)                # next point to prev point dist
      bmidpt <- (from + to1) / 2    # the before-mid-point; average of the point and its preceeding point
      amidpt <- (from + to2) / 2    # the after-mid-point; average of the point and its preceeding point

      # identify phi corresponding to each mid-point after transformed back to cr coordinates
      # angles (b: before, a: after) represent radians off of horizontal for vector from MLE to CR point (0-90 degrees)
      boffangle <- atan(abs((bmidpt[,2] / transf.scale - theta2.hat) / (bmidpt[,1] - theta1.hat)))
      aoffangle <- atan(abs((amidpt[,2] / transf.scale - theta2.hat) / (amidpt[,1] - theta1.hat)))

      # bphi assumes phi values for "before" points as they relate to the original, cr confidence region scale
      # when the elbows-algorithm below determines a point before &/or after is required, these points are inserted
      bphi <-
        ifelse((bmidpt[,1] >  theta1.hat) & (bmidpt[,2] >= (transf.scale * theta2.hat)), 1, 0) * boffangle +             # 1st quad
        ifelse((bmidpt[,1] <= theta1.hat) & (bmidpt[,2] >  (transf.scale * theta2.hat)), 1, 0) * ((pi) - boffangle) +    # 2nd quad
        ifelse((bmidpt[,1] <  theta1.hat) & (bmidpt[,2] <= (transf.scale * theta2.hat)), 1, 0) * ((pi) + boffangle) +    # 3rd quad
        ifelse((bmidpt[,1] >= theta1.hat) & (bmidpt[,2] <  (transf.scale * theta2.hat)), 1, 0) * ((2 * pi) - boffangle)  # 4th quad

      # aphi assumes phi values for "after" points as they relate to the original, cr confidence region scale
      # when the elbows-algorithm below determines a point before &/or after is required, these points are inserted
      aphi <-
        ifelse((amidpt[,1] >  theta1.hat) & (amidpt[,2] >= (transf.scale * theta2.hat)), 1, 0) * aoffangle +             # 1st quad
        ifelse((amidpt[,1] <= theta1.hat) & (amidpt[,2] >  (transf.scale * theta2.hat)), 1, 0) * ((pi) - aoffangle) +    # 2nd quad
        ifelse((amidpt[,1] <  theta1.hat) & (amidpt[,2] <= (transf.scale * theta2.hat)), 1, 0) * ((pi) + aoffangle) +    # 3rd quad
        ifelse((amidpt[,1] >= theta1.hat) & (amidpt[,2] <  (transf.scale * theta2.hat)), 1, 0) * ((2 * pi) - aoffangle)  # 4th quad

      elbowcalc <- (d01^2 + d02^2 - d12^2)/(2 * d01 * d02)              # law of cosines acos(__) parameter
      if (min(elbowcalc) < -1) {
        #print(paste0("Warning: roundoff issue ~ -1 suspected and reset to -1; acos(value) with value approx:", min(elbowcalc)))
        elbowcalc <- ifelse(elbowcalc < -1, -1, elbowcalc)
      }
      if (max(elbowcalc) > 1) {
        #print(paste0("WARNING: roundoff issue ~ 1 suspected and reset to 1; acos(value) with value approx:", max(elbowcalc)))
        elbowcalc <- ifelse(elbowcalc > 1, 1, elbowcalc)
      }
      elbows <- pi - acos(elbowcalc)
      maxcurrent <- max(elbows)
      #print(paste0("iteration of new phi's: ", count))
      #print(paste0("max current: ", maxcurrent))
      #toobig <- ifelse(elbows > maxrad, 1, 0)
      #print(paste0("number of out-of-tollerance points: ", sum(toobig)))

      good <- ifelse(elbows <= maxrad, 1, 0)
      goodminus <- c(good[length(good)], good[1:(length(good) - 1)])
      goodplus <- c(good[2:length(good)], good[1])
      goodtot <- good + 2 * goodplus + 4 *goodminus
      # for 'goodtot' above:
      # 7 = current, previous, and next points all valid                     (insert 0 points)
      # 6 = previous and next points   valid, but current point invalid      (insert 2 points)
      # 5 = current and previous point valid, next point invalid             (insert 0 points)
      # 4 = previous point             valid, current and next invalid       (insert 1 after)
      # 3 = current and next points    valid, previous point invalid         (insert 0 points)
      # 2 = next point is              valid, current and previous invalid   (insert 1 before)
      # 1 = current point              valid, previous and next are invalid  (insert 0 points)
      # 0 = current, previous, and next points all invalid                   (insert 2 points)
      phinew <- c(0)
      for (i in 1:(length(phi))) {                          # cycle through points, adding phi where necessary
        phinewi <- c(0)
        if (goodtot[i] == 0 | goodtot[i] == 6) {            # insert 2 points
          phinewi <- c(bphi[i], aphi[i])
        }
        else if (goodtot[i] == 2) {                         # insert 1 point before
          phinewi <- c(bphi[i])
        }
        else if (goodtot[i] == 4) {                         # insert 1 point after
          phinewi <- c(aphi[i])
        }

        if (sum(phinewi) != 0) {
          if (sum(phinew) == 0) { phinew <- phinewi }       # if first entry in phinew
          else { phinew <- c(phinew, phinewi) }             # o.w. augment to existing phinew
        }
      }
      crlist <- list("x" = cr[, 1], "y" = cr[, 2], "phi" = phi)
      phinew <- sort(unique(phinew))
      if (sum(phinew) != 0) {phi <- sort(unique(c(phi, phinew))) }
      if (count >= 50) {
        maxcurrent <- maxrad
        warning("max iteration tolerance hit; working confidence region plot is shown, however, maxdeg constraint is not met")
        #print("WARNING: max iteration tolerance hit; maxdeg constraint not met, however, working plot is shown.")
        #print("...this is attributable to either regions inaccessible via a radial azimuth from its MLE, or a uniroot tol argument too small.")
      }
    }
    invisible(crlist)
  }

  ############################################################################################
  ### end of functions; function calls below determine and then plot the confidence region ###
  ############################################################################################

  #ballbearing <- c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.48, 51.84, 51.96, 54.12, 55.56,
  #  67.80, 68.64, 68.64, 68.88, 84.12, 93.12, 98.64, 105.12, 105.84, 127.92, 128.04, 173.40)
  cen <- c(rep(1, length(dataset)))
  mydata <- dataset
  mle.list <- mlesolve(mydata, cen = cen)
  theta1.hat <- mle.list$theta1.hat
  theta2.hat <- mle.list$theta2.hat
  mleLLvalue <- mle.list$mleLLvalue

  # identify confidence region boundary coordinates
  if (heuristic == 1) {
    crlist <- crsmooth(maxdeg = maxdeg, samp = mydata, cen = cen, alpha = alpha,
                       mle.list = mle.list, ellipse_n = ellipse_n)
  }
  else if (heuristic == 0) {
    phi <- angles(a = theta1.hat, b = theta2.hat, npoints = ellipse_n)
    cr <- crsolve(samp = dataset, cen = cen, alpha = alpha, mle.list = mle.list, phi = phi)
    crlist <- list("x" = cr[, 1], "y" = cr[, 2], "phi" = phi)
  }

  # assemble conf region points and "close" the gap between first & last points graphically
  cr <- cbind(crlist$x, crlist$y)
  cr <- rbind(cr, cr[1, ])
  phi <- crlist$phi

  if (showplot == TRUE) {
  # construct the confidence region plot
  # override default plot characteristics if user specified
  # i.e. mar, xlim, ylim, main, xlab, ylab, etc.

  # margin adjustments
  par(xpd = F)
  if (!is.null(mar)) {
    par(mar = mar)
  }

  # axis tick-marks based on mlelab and origin
  xlabs <- c(min(cr[, 1]), max(cr[, 1]))
  ylabs <- c(min(cr[, 2]), max(cr[, 2]))
  if (mlelab == TRUE) {
    xlabs <- c(xlabs, theta1.hat)
    ylabs <- c(ylabs, theta2.hat)
  }
  if (origin == TRUE) {
    xmin <- 0
    ymin <- 0
    xlabs <- c(0, xlabs)
    ylabs <- c(0, ylabs)
  }
  else {
    xmin <- min(cr[,1])
    ymin <- min(cr[,2])
  }

  # axis limits and tick-marks if xlim and/or ylim entries:
  if(is.null(xlim)) {
    xlim <- c(xmin, max(cr[,1]))
  }
  else {
    xlabs <- c(xlim, xlabs)
    xlabs <- xlabs[xlabs >= xlim[1] & xlabs <= xlim[2]]
  }
  if(is.null(ylim)) {
    ylim <- c(ymin, max(cr[,2]))
  }
  else {
    ylabs <- c(ylim, ylabs)
    ylabs <- ylabs[ylabs >= ylim[1] & ylabs <= ylim[2]]
  }

  # plot
    plot(cr, xlab = xlab, ylab = ylab, main = main, ylim = ylim, xlim = xlim,
         axes = FALSE, type = 'l')
    my.atx <- unique(round(xlabs, sf[1]))
    my.aty <- unique(round(ylabs, sf[2]))
    axis(side = 1, at = my.atx, las = xlas)
    axis(side = 2, at = my.aty, las = ylas)
    if (mlelab == TRUE) {
      points(theta1.hat, theta2.hat, pch = 3)
    }
    if (pts == TRUE) points(unique(cr), lwd = 0.65)
    # enable plotting beyond margins for any post-processing user add-ons
    par(xpd = TRUE)
  }

  #print("dataset was: ")
  #print(dataset)
  # returned values (if requested), otherwise only output to screen # boundary points.
  if (length(dataset) <= 5) {
    warning("small sample size is ill-suited to invoke the asymptotic properties assumed by the confidence region plot")
    #print("WARNING: small sample sizes can yield irregular confidence region shapes that are unatainable using crplot.")
    #print("WARNING: small sample sizes violate the asymptotic properties assumption used to produce the confidence region.")
  }
  if (info == TRUE) {
    print(paste0("Confidence region plot complete; made using ", length(phi)," boundary points."))
    if (distn == "weibull") {
      print(paste0("MLE value is: (kappa.hat = ", theta1.hat, ", lambda.hat = ", theta2.hat,")"))
      crlist <- list("kappa" = crlist$x, "lambda" = crlist$y, "phi" = crlist$phi)
    }
    else if (distn == "invgauss") {
      print(paste0("MLE value is: (mu.hat = ", theta1.hat, ", lambda.hat = ", theta2.hat,")"))
      crlist <- list("mu" = crlist$x, "lambda" = crlist$y, "phi" = crlist$phi)
    }
  return(crlist)
  }
  else if (info == FALSE) {
    return(paste0("Confidence region plot complete; made using ", length(phi)," boundary points."))
  }
}
