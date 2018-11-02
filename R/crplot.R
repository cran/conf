#' Plotting Two-Dimensional Confidence Regions
#'
#' @description
#' Plots the two-dimensional confidence region for probability distribution parameters (supported distribution
#' suffixes: cauchy, gamma, invgauss, lnorm, llogis, logis, norm, unif, weibull) corresponding to a user given
#' dataset and level of significance.  See the CRAN website https://CRAN.R-project.org/package=conf for a link
#' to a \code{crplot} vignette.
#'
#' @param dataset a 1 x n vector of dataset values.
#' @param alpha significance level; resulting plot illustrates a 100(1 - alpha)\% confidence region.
#' @param distn distribution to fit the dataset to; accepted values: \code{'cauchy'}, \code{'gamma'}, \code{'invgauss'},
#' \code{'logis'}, \code{'llogis'}, \code{'lnorm'}, \code{'norm'}, \code{'unif'}, \code{'weibull'}.
#' @param cen a vector of binary values specifying right-censored values as 0, and 1 (default) otherwise
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
#' @param xyswap logical argument to switch the axes that the distribution parameter are shown.
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
#' @param repair logical argument to repair regions inaccessible using a radial angle from its MLE due to multiple
#' roots at select \eqn{\phi} angles.
#' @param jumpshift \% of vertical or horizontal "gap" (near uncharted region) to angle jump-center location towards
#' @param jumpuphill significance level increase to alpha ("uphill" on confidence region) to locate the jump-center
#' @param showplot logical argument specifying if a plot is output; altering from its default of \code{TRUE} is
#' only logical assuming \code{crplot} is run for its data only (see the \code{info} argument).
#' @import stats
#' @import graphics
#' @importFrom fitdistrplus mledist
#' @importFrom STAR llogisMLE gammaMLE
#' @export
#' @return if the optional argument \code{info = TRUE} is included then a list of plot coordinates and phi angles is returned
#' @concept confidence region plot
#' @keywords confidence region, confidence intervals, statistical graphics, data visualization, graphical methods,
#' parameter estimation, numerical optimization
#' @references Jaeger, A. (2016), "Computation of Two- and Three-Dimensional Confidence Regions with the Likelihood Ratio",
#' The American Statistician, 49, 48--53.
#' @references Weld, C., Loh, A., Leemis, L. (in press), "Plotting Likelihood-Ratio Based Confidence Regions for
#' Two-Parameter Univariate Probability Models, The American Statistician.
#' @seealso \code{\link{coversim}}, \code{\link{uniroot}}
#' @author Christopher Weld (\email{ceweld@email.wm.edu})
#' @author Lawrence Leemis (\email{leemis@math.wm.edu})
#'
#' @usage
#' crplot(dataset, alpha, distn,
#'                 cen       = rep(1, length(dataset)),
#'                 heuristic = 1,
#'                 maxdeg    = 5,
#'                 ellipse_n = 4,
#'                 pts       = TRUE,
#'                 mlelab    = TRUE,
#'                 sf        = c(5, 5),
#'                 mar       = c(4, 5, 2, 1.5),
#'                 xyswap    = FALSE,
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
#'                 repair    = TRUE,
#'                 jumpshift = 0.5,
#'                 jumpuphill = min(alpha, 0.01),
#'                 showplot  = TRUE )
#'
#' @details
#' This function plots confidence regions for a variety of two-parameter distributions.  It requires
#' \itemize{
#' \item a vector of dataset values,
#' \item the level of significance (alpha), and
#' \item a distribution to fit the data to.
#' }
#' Plots display according to probability density function parameterization given later in this section.
#' Two heuristics (and their associated combination) are available to plot confidence regions.  Along
#' with their descriptions, they are:
#' \enumerate{
#' \item \emph{Smoothing Boundary Search Heuristic (default)}.  This heuristic plots more points in areas of
#' greater curvature to ensure a smooth appearance throughout the confidence region boundary.  Its
#' \code{maxdeg} parameter specifies the maximum tolerable angle between three successive points.
#' Lower values of \code{maxdeg} result in smoother plots, and its default value of 5 degrees
#' provides adequate smoothing in most circumstances.  Values of \code{maxdeg} \eqn{\le} 3 are not
#' recommended due to their complicating implications to trigonometric numerical approximations near 0
#' and 1; their use may result in plot errors.
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
#' The default axes convention in use by \code{crplot} are
#'
#' \tabular{lcc}{
#' \tab Horizontal \tab Vertical\cr
#' Distribution  \tab  Axis  \tab Axis\cr
#' Caucy \tab \eqn{a} \tab \eqn{s}\cr
#' gamma \tab \eqn{\theta} \tab \eqn{\kappa}\cr
#' inverse Gaussian \tab \eqn{\mu} \tab \eqn{\lambda}\cr
#' log logistic \tab \eqn{\lambda} \tab \eqn{\kappa}\cr
#' log normal \tab \eqn{\mu} \tab \eqn{\sigma}\cr
#' logistic \tab \eqn{\mu} \tab \eqn{\sigma}\cr
#' normal \tab \eqn{\mu} \tab \eqn{\sigma}\cr
#' uniform \tab \eqn{a} \tab \eqn{b}\cr
#' Weibull \tab \eqn{\kappa} \tab \eqn{\lambda}
#' }
#'
#' where each respective distribution is defined below.
#'
#' \itemize{
#' \item The Cauchy distribution
#' for the real-numbered location parameter \eqn{a}, scale parameter \eqn{s}, and \eqn{x} is a real number,
#' has the probability density function
#' \deqn{1 / (s \pi (1 + ((x - a) / s) ^ 2)).}
#'
#' \item The gamma distribution
#' for shape parameter \eqn{\kappa > 0}, scale parameter \eqn{\theta > 0}, and \eqn{x > 0},
#' has the probability density function
#' \deqn{1 / (Gamma(\kappa) \theta ^ \kappa) x ^ {(\kappa - 1)} exp(-x / \theta).}
#'
#' \item The inverse Gaussian distribution
#' for mean \eqn{\mu > 0}, shape parameter \eqn{\lambda > 0}, and \eqn{x > 0},
#' has the probability density function
#' \deqn{\sqrt (\lambda / (2 \pi x ^ 3)) exp( - \lambda (x - \mu) ^ 2 / (2 \mu ^ 2 x)).}
#'
#' \item The log logistic distribution
#' for scale parameter \eqn{\lambda > 0}, shape parameter \eqn{\kappa > 0}, and \eqn{x \ge 0},
#' has a probability density function
#' \deqn{(\kappa \lambda) (x \lambda) ^ {(\kappa - 1)} / (1 + (\lambda x) ^ \kappa) ^ 2.}
#'
#' \item The log normal distribution
#' for the real-numbered mean \eqn{\mu} of the logarithm, standard deviation \eqn{\sigma > 0}
#' of the logarithm, and \eqn{x > 0},
#' has the probability density function
#' \deqn{1 / (x \sigma \sqrt(2 \pi)) exp(-(\log x - \mu) ^ 2 / (2 \sigma ^ 2)).}
#'
#' \item The logistic distribution
#' for the real-numbered location parameter \eqn{\mu}, scale parameter \eqn{\sigma}, and \eqn{x} is a real number,
#' has the probability density function
#' \deqn{(1 / \sigma) exp((x - \mu) / \sigma) (1 + exp((x - \mu) / \sigma)) ^ -2}
#'
#' \item The normal distribution
#' for the real-numbered mean \eqn{\mu}, standard deviation \eqn{\sigma > 0}, and \eqn{x} is a real number,
#' has the probability density function
#' \deqn{1 / \sqrt (2 \pi \sigma ^ 2) exp(-(x - \mu) ^ 2 / (2 \sigma ^ 2)).}
#'
#' \item The uniform distribution for real-valued parameters \eqn{a} and \eqn{b} where \eqn{a < b}
#' and \eqn{a \le x \le b},
#' has the probability density function
#' \deqn{1 / (b - a).}
#'
#' \item The Weibull distribution
#' for scale parameter \eqn{\lambda > 0}, shape parameter \eqn{\kappa > 0}, and \eqn{x > 0},
#' has the probability density function
#' \deqn{\kappa (\lambda ^ \kappa) x ^ {(\kappa - 1)} exp(-(\lambda x) ^ \kappa).}
#' }
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
# Latest Revision  : May 2018
#######################################################################################
crplot <- function(dataset,
                   alpha,
                   distn,
                   cen = rep(1, length(dataset)),
                   heuristic = 1,
                   maxdeg = 5,
                   ellipse_n = 4,
                   pts = TRUE,
                   mlelab = TRUE,
                   sf = c(5, 5),
                   mar = c(4, 5, 2, 1.5),
                   xyswap = FALSE,
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
                   repair = TRUE,
                   jumpshift = 0.5,
                   jumpuphill = min(alpha, 0.01),
                   showplot = TRUE) {

  # parameter error checking ###########################################################

  if (missing(dataset)) stop ("argument 'dataset' is missing, with no default")
  if (missing(alpha)) stop ("argument 'alpha' is missing, with no default")
  if (missing(distn)) stop ("argument 'distn' is missing, with no default")

  if (is.null(dataset) || length(dataset) == 1 || !is.numeric(dataset))
    stop("'dataset' must be a numeric vector with length > 1")

  if (!is.element(distn, c("weibull", "invgauss", "norm", "lnorm", "logis", "llogis", "gamma", "unif", "cauchy")))
    stop("'distn' invalid; supported are: 'cauchy', 'gamma', 'invgauss', 'lnorm', 'logis', 'llogis', 'norm', 'unif', 'weibull'")

  if (distn == "weibull" && min(dataset) <= 0)
    stop("'dataset' parameter contains infeasible Weibull distribution outcome(s) <= 0")

  if (distn == "invgauss" && min(dataset) <= 0)
    stop("'dataset' parameter contains infeasible inverse Gaussian distribution outcome(s) <= 0")

  if (distn == "llogis" && min(dataset) <= 0)
    stop("'dataset' parameter contains infeasible loglogistic distribution outcome(s) <= 0")

  if (distn == "gamma" && min(dataset) <= 0)
    stop("'dataset' parameter contains infeasible gamma distribution outcome(s) <= 0")

  if (is.null(alpha) || alpha <= 0 || alpha >= 1 || !is.numeric(alpha) || length(alpha) != 1)
    stop("'alpha' numeric significance level parameter required such that 0 < alpha < 1")

  if (!is.numeric(cen) || !all(cen %in% 0:1))
    stop("'cen' must be a vector of binary (0 or 1) values with length(dataset) entries")

  if ((distn == "unif") && (sum(cen) != length(dataset))) {
    if (max(as.numeric(dataset[which(cen == 0)] %in% max(dataset))) == 1) {
      stop("undefined 'unif' confidence region when max(dataset) corresponds to a (cen = 0) censored value")
    }
  }

  if (is.null(heuristic) || !is.numeric(heuristic) || (heuristic != 0 && heuristic != 1))
    stop("'heuristic' parameter must be 0 (elliptic-oriented points) or 1 (search heuristic)")

  if (is.null(maxdeg) || !is.numeric(maxdeg) || maxdeg < 1 || maxdeg >= 90)
    stop("'maxdeg' parameter must be a numeric angle tollerance in degrees such that 1 <= maxdeg < 90")

  if ((ellipse_n == 4) && (heuristic == 0))
    stop("'ellipse_n' (number of plot points) >= 8 required for 'heuristic = 0' (elliptic-oriented points)")

  if (!is.numeric(ellipse_n) || length(ellipse_n) != 1 || (ellipse_n) %% 4 != 0 || ellipse_n <= 0 )
    stop("'ellipse_n' must be a positive numeric multiple of 4")

  if (!is.logical(pts) || length(pts) != 1)
    stop("'pts' must be a single logical parameter")

  if (!is.logical(mlelab) || length(mlelab) != 1)
    stop("'mlelab' must be a single logical parameter")

  if (length(sf) != 2 || !is.numeric(sf) || floor(sf)[1] != sf[1] || floor(sf)[2] != sf[2])
    stop("'sf' must be a vector of integers with length two")

  if (length(mar) != 4 || !is.numeric(mar) || min(mar) < 0)
    stop("'mar' must be a vector of length four with positive numeric entries")

  if (!is.logical(xyswap) || length(xyswap) != 1)
    stop("'xyswap' must be a single logical parameter")

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

  if (!is.logical(repair) || length(repair) != 1)
    stop("'repair' must be a single logical parameter")

  if (jumpshift <= 0 || jumpshift >= 1 || !is.numeric(jumpshift) || length(jumpshift) != 1)
    stop("'jumpshift' must be a single numeric value 0 < jumpshift < 1")

  if (jumpuphill <= 0 || !is.numeric(jumpuphill) || length(jumpuphill) != 1)
    stop("'jumpuphill' must be a single numeric value such that 0 < (alpha + jumpuphill) < 1")

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
  #                    in Leemis Reliability text).
  #
  # 4. asesolve.       This function calculates and returns the asymptotic standard error for the specified
  #                    distribution.  Returned is ase.list with values list("theta1.ase", "theta2.ase").
  #
  # 5. llrsolve.       Using the asymptotically chisquared likelihood ratio statistic, this function returns
  #                    a function whose positive values reside within the (1 - alpha)% confidence region,
  #                    and negative values reside outside of it.
  #
  # 6. crsolve.        This function, leveraging the llrsolve function, calls uniroot to determine the
  #                    cartesian coordinates of confidence region boundary points corresponding to the input
  #                    vector of phi angles (given with respect to the MLE).
  #
  # 7. crsmooth.       Confidence region smoothing function.  Given an initial vector of angles to
  #                    calculate confidence region points with, it iteratively adds additional angles
  #                    until all apparent angle bends are within the specified angle tolerance.  Also
  #                    contains a recursive routine to identify and fix "deadspace" inaccessible by a
  #                    radial angle from the MLE by identifying and smoothing the CR from an alternate
  #                    off-MLE jump-center location within the CR.
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
  # Gamma and log logistic use the STAR package to identify their MLEs.
  # Returned is mle.list with values list("theta1.hat", "theta2.hat", "mleLLvalue").
  mlesolve = function(x, cen, epsilon = 0.000000001){
    n <- length(x)
    r <- sum(cen)

    # Cauchy MLE
    if (distn == "cauchy"){
      # with censoring:
      left <- replace(x, which(cen == -1), rep(NA, length(which(cen == -1))))
      right <- replace(x, which(cen == 0), rep(NA, length(which(cen == 0))))
      xx <- data.frame(left = left, right = right)
      xxmle <- fitdistrplus::mledist(xx, "cauchy", silent = TRUE)
      mleLL <- xxmle$loglik
      a.hat <- as.numeric(xxmle$estimate[1][1])
      alpha.hat <- as.numeric(xxmle$estimate[2][1])
      mle.list <- list("theta1.hat" = a.hat, "theta2.hat" = alpha.hat, "mleLLvalue" = mleLL)
    }

    # Weibull MLE
    if (distn == "weibull"){
      # c0 is an educated guess for the initial Menon estimate to begin iterating through the algorithm
      # reference: Appendix B (pg 338) in Leemis Reliability text
      temp1 <- sum(log(x) ^ 2)
      temp2 <- (sum(log(x)) ^ 2)
      temp3 <- 6 / ((n - 1) * pi ^ 2)
      temp4 <- (temp3 * (temp1 - temp2 / n))
      c0 <- 1 / sqrt(temp4)                     # initial estimate (ref: p 338 Reliability)

      # Fixed point algorithm (FPA) supporting function for the Weibull MLE
      # given by Qiao and Tsokos in their article Estimation of the three parameter
      # Weibull probability distribution in Mathematics and Computers in Simulation,
      # 1995, pgs 173--185
      # (two-parameter summary is given on pg 174--175)
      s1 <- sum(log(x) * cen)
      repeat{                                   # ref: bottom p 246 Reliability
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
      # using STAR package for Weibull MLE calcluations generates warning messages:
      #kappa.hat <- weibullMLE(x)$estimate[1]
      #lambda.hat <- 1 / weibullMLE(x)$estimate[2]

      # compute log likelihood function for censored weibull
      # reference: pg 246 in Leemis Reliability text
      temp1 <- sum(x ^ kappa.hat)
      temp2 <- sum(cen * log(x))
      mleLL <- (r * log(kappa.hat) + (kappa.hat * r * log(lambda.hat)) + ((kappa.hat - 1) * temp2) -
                     (lambda.hat ^ kappa.hat * temp1))
      mle.list <- list("theta1.hat" = kappa.hat, "theta2.hat" = lambda.hat, "mleLLvalue" = mleLL)
    }

    # inverse Gaussian MLE
    else if (distn == "invgauss"){
      # original formulation (without censored values)
      #mu.hat <- sum(x) / length(x)
      #tempsum <- sum(1 / x)
      #lambda.hat <-  ((1 / (length(x))) * tempsum -
      #                  (length(x) / sum(x))) ^ -1
      # original loglikelihood value at MLE (without censored values)
      #temp1 <- sum(log(x))
      #temp2 <- sum((x - mu.hat) ^ 2 / x)
      #mleLL <- (n / 2) * log(lambda.hat) - (n / 2) * log(2 * pi) - (3 / 2) * temp1 -
      #  (lambda.hat / (2 * mu.hat ^ 2)) * temp2

      # using STAR package in order to incorporate censored values:
      # note: conf shape parameter (lambda) = (1 / sigma2) from STAR package
      temp1 <- STAR::invgaussMLE(x, si = as.double(cen))     # STAR requires cen package as double
      mu.hat <- temp1$estimate[['mu']]
      lambda.hat <- 1 / temp1$estimate[['sigma2']]
      mleLL <- temp1$logLik                      # log likelihood value at its maximum
      mle.list <- list("theta1.hat" = mu.hat, "theta2.hat" = lambda.hat, "mleLLvalue" = mleLL)
    }

    # normal MLE
    else if (distn == "norm"){
      # without censoring:
      #n <- length(x)
      #mu.hat <- (1 / n) * sum(x)
      #sigma2.hat <- (1 / n) * sum((x - mu.hat)^2)  # this is sigma-squared, but will return sigma
      #mleLL <- (-n / 2)*log(2*pi) - (n / 2)*log(sigma2.hat) - (1 / (2 * sigma2.hat)) * sum((x - mu.hat)^2)
      #mle.list <- list("theta1.hat" = mu.hat, "theta2.hat" = sqrt(sigma2.hat), "mleLLvalue" = mleLL)
      # with censoring:
      left <- replace(x, which(cen == -1), rep(NA, length(which(cen == -1))))
      right <- replace(x, which(cen == 0), rep(NA, length(which(cen == 0))))
      xx <- data.frame(left = left, right = right)
      xxmle <- fitdistrplus::mledist(xx, "norm", silent = TRUE)
      mleLL <- xxmle$loglik
      mu.hat <- xxmle$estimate[1][1]
      sigma.hat <- xxmle$estimate[2][1]
      mle.list <- list("theta1.hat" = mu.hat, "theta2.hat" = sigma.hat, "mleLLvalue" = mleLL)
    }

    # lognormal MLE
    else if (distn == "lnorm"){
      # without censoring:
      #n <- length(x)
      #mu.hat <- (1 / n) * sum(log(x))
      #sigma2.hat <- (1 / n) * sum((log(x) - mu.hat)^2)  # this is sigma-squared, but will return sigma
      #mleLL <- (-n / 2)*log(2*pi) - (n / 2)*log(sigma2.hat) - sum(log(x)) - (1 / (2 * sigma2.hat)) * sum((log(x) - mu.hat)^2)
      #mle.list <- list("theta1.hat" = mu.hat, "theta2.hat" = sqrt(sigma2.hat), "mleLLvalue" = mleLL)
      # with censoring:
      left <- replace(x, which(cen == -1), rep(NA, length(which(cen == -1))))
      right <- replace(x, which(cen == 0), rep(NA, length(which(cen == 0))))
      xx <- data.frame(left = left, right = right)
      xxmle <- fitdistrplus::mledist(xx, "lnorm", silent = TRUE)
      mleLL <- xxmle$loglik
      mu.hat <- xxmle$estimate[1][1]
      sigma.hat <- xxmle$estimate[2][1]
      mle.list <- list("theta1.hat" = mu.hat, "theta2.hat" = sigma.hat, "mleLLvalue" = mleLL)
    }

    # logistic MLE
    else if (distn == "logis"){
      # with censoring:
      left <- replace(x, which(cen == -1), rep(NA, length(which(cen == -1))))
      right <- replace(x, which(cen == 0), rep(NA, length(which(cen == 0))))
      xx <- data.frame(left = left, right = right)
      xxmle <- fitdistrplus::mledist(xx, "logis", silent = TRUE)
      mleLL <- xxmle$loglik
      mu.hat <- xxmle$estimate[1][1]        # location
      sigma.hat <- xxmle$estimate[2][1]     # scale
      mle.list <- list("theta1.hat" = mu.hat, "theta2.hat" = sigma.hat, "mleLLvalue" = mleLL)
    }

    # loglogistic MLE
    # admin note: STAR expresses the log logistic distribution with a much different parameterization
    # although a counterintuitive equality, its "location" parameter = 1 / log(scale) via our parameterization
    # and its "scale" parameter = (1 / shape) via our parameterization
    else if (distn == "llogis"){
      temp1 <- STAR::llogisMLE(x, si = as.double(cen))
      lambda.hat <- 1 / exp(temp1$estimate[['location']])
      kappa.hat <- 1 / temp1$estimate[['scale']]
      mleLL <- temp1$logLik                       # log likelihood value at its maximum
      mle.list <- list("theta1.hat" = lambda.hat, "theta2.hat" = kappa.hat, "mleLLvalue" = mleLL)
    }

    # gamma MLE
    else if (distn == "gamma"){
      temp1 <- STAR::gammaMLE(x, si = as.double(cen))
      kappa.hat <- temp1$estimate[['shape']]
      theta.hat <- temp1$estimate[['scale']]
      mleLL <- temp1$logLik                      # log likelihood value at its maximum
      mle.list <- list("theta1.hat" = theta.hat, "theta2.hat" = kappa.hat, "mleLLvalue" = mleLL)
    }

    # uniform MLE
    else if (distn == "unif") {
      a.hat <- min(x)
      b.hat <- max(x)
      #mleLL <- -n * log(b.hat - a.hat)
      mleLL <- sum(log(b.hat - x)[cen == 0]) - n * log(b.hat - a.hat)  # with censoring
      mle.list <- list("theta1.hat" = a.hat, "theta2.hat" = b.hat, "mleLLvalue" = mleLL)
    }

    invisible(mle.list)
  }


  ######################################################################################

  # asesolve -------------------------------------------------------------------
  # This function calculates and returns the asymptotic standard error for the specified distribution
  # Returned is ase.list with values list("theta1.ase", "theta2.ase").
  # Distns without std error easily accessible use the relative MLE sizes to estimate the aspect ratio
  asesolve = function(x, cen, theta1.hat, theta2.hat) {
    #print("entering asesolve")
    n <- length(x)
    r <- sum(cen)
    if (distn == "cauchy") {
      left <- replace(x, which(cen == -1), rep(NA, length(which(cen == -1))))
      right <- replace(x, which(cen == 0), rep(NA, length(which(cen == 0))))
      xx <- data.frame(left = left, right = right)
      hes <- fitdistrplus::mledist(xx, "cauchy", silent = TRUE)$hessian      # using plotdistrplus
      OI <- solve(hes)
      se <- sqrt(OI)
      theta1ase <- se[1]    # a location parameter
      theta2ase <- se[4]    # alpha scale parameter
    }
    else if (distn == "gamma") {
      theta1ase <- as.numeric(STAR::gammaMLE(x, si = as.double(cen))$se[2])    # scale parameter theta
      theta2ase <- as.numeric(STAR::gammaMLE(x, si = as.double(cen))$se[1])    # shape parameter kappa
    }
    else if (distn == "invgauss") {
      theta1ase <- as.numeric(STAR::invgaussMLE(x, si = as.double(cen))$se[1])       # mu
      inv_theta2ase <- as.numeric(STAR::invgaussMLE(x, si = as.double(cen))$se[2])   # (1 / theta) se
      theta2ase <- abs(1 / ((1 / theta2.hat) + (inv_theta2ase / 2)) -     # theta se approximation
                         1 / ((1 / theta2.hat) - (inv_theta2ase / 2)))
    }
    else if (distn == "norm") {
      left <- replace(x, which(cen == -1), rep(NA, length(which(cen == -1))))
      right <- replace(x, which(cen == 0), rep(NA, length(which(cen == 0))))
      xx <- data.frame(left = left, right = right)
      hes <- fitdistrplus::mledist(xx, "norm", silent = TRUE)$hessian       # using plotdistrplus
      OI <- solve(hes)
      se <- sqrt(OI)
      theta1ase <- se[1]    # mu
      theta2ase <- se[4]    # sigma
    }
    else if (distn == "logis") {
      left <- replace(x, which(cen == -1), rep(NA, length(which(cen == -1))))
      right <- replace(x, which(cen == 0), rep(NA, length(which(cen == 0))))
      xx <- data.frame(left = left, right = right)
      hes <- fitdistrplus::mledist(xx, "logis", silent = TRUE)$hessian       # using plotdistrplus
      OI <- solve(hes)
      se <- sqrt(OI)
      theta1ase <- se[1]    # mu (location)
      theta2ase <- se[4]    # sigma (scale)
    }
    else if (distn == "llogis") {
      alt_param <- STAR::llogisMLE(x, si = as.double(cen))
      alt_theta1ase <- as.numeric(alt_param$se[1])     # (1 / log(lambda)) se
      inv_theta2ase <- as.numeric(alt_param$se[2])     # (1 / kappa) se
      high1 <- as.numeric(alt_param$estimate[1] + alt_theta1ase * 0.5)
      low1 <- as.numeric(alt_param$estimate[1] - alt_theta1ase * 0.5)
      theta1ase <- abs((1 / exp(high1)) - (1 / exp(low1)))              # approximation IAW param conversion
      theta2ase <- abs(1 / ((1 / theta2.hat) + (inv_theta2ase / 2)) -   # approximation IAW param conversion
                         1 / ((1 / theta2.hat) - (inv_theta2ase / 2)))
    }
    else if (distn == "lnorm") {
      left <- replace(x, which(cen == -1), rep(NA, length(which(cen == -1))))
      right <- replace(x, which(cen == 0), rep(NA, length(which(cen == 0))))
      xx <- data.frame(left = left, right = right)
      hes <- fitdistrplus::mledist(xx, "lnorm", silent = TRUE)$hessian      # using plotdistrplus
      OI <- solve(hes)
      se <- sqrt(OI)
      theta1ase <- se[1]    # mu
      theta2ase <- se[4]    # sigma
    }
    else if (distn == "weibull") {
      khat <- theta1.hat
      lhat <- theta2.hat
      # 2nd partial derivatives at MLE (ref: pg 247 Reliability text):
      p2kap <- r / (khat ^ 2) + sum((lhat * x) ^ khat * (log(lhat * x))^2)
      p2lam <- khat * r / (lhat ^ 2) + khat *(khat - 1) * lhat ^ (khat - 2) * sum(x ^ khat)
      p2lamkap <- -length(x) / lhat + (lhat ^ (khat - 1)) * (khat * sum(x ^ khat * log(x)) +
                                                               (1 + khat * log(lhat)) * sum(x ^ khat))
      # asymptotic standard error
      hes <- matrix(c(p2kap, p2lamkap, p2lamkap, p2lam), nrow = 2)
      OI <- solve(hes)
      se_kappa <- sqrt(OI[1])
      se_lambda <- sqrt(OI[4])
      theta1ase <- se_kappa
      theta2ase <- se_lambda
    }
    ase.list <- list("theta1.ase" = theta1ase, "theta2.ase" = theta2ase)
    invisible(ase.list)
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
    n <- length(x)                     # sample size (censored and uncensored)
    r <- sum(cen)                      # number of uncensored values
    temp <- mleLLvalue - (chi2 / 2)    # for later use with respect to loglikelihood ratio eqn

    # Cauchy distribution
    if (distn == "cauchy") {
      a.hat <- MLEHAT[1, ]
      alpha.hat <- MLEHAT[2, ]
      a <- a.hat + (d * cos(phi))
      alpha <- alpha.hat + (d * sin(phi))
      # this option isn't cooperative (possibly because fitditrplus functions?)
      #llfn <- sum(cen * log(fitdistrplus::dcauchy(x, scale = alpha, location = a))) +
      #  sum(as.numeric(cen==0) * log((1 - fitdistrplus::pnorm(x, scale = alpha, location = a))))
      # written out log likelihood function
      #llfn <- sum(cen * (log(2 * alpha * (alpha ^ 2 - 2 * x * a + a ^ 2 + x ^ 2) ^ (-1) * (pi + 2 * atan((a - x) / alpha)) ^ (-1)))) -
      #               sum(log(2 * pi) - log(pi + 2 * atan((a - x) / alpha)))
      llfn <- sum(cen * (log(2 * alpha))) - sum(cen * log(alpha ^ 2 - 2 * x * a + a ^ 2 + x ^ 2)) - sum(cen * log(pi + 2 * atan((a - x) / alpha))) -
        sum(log(2 * pi) - log(pi + 2 * atan((a - x) / alpha)))
      valuereturn <- llfn - temp   # ID d such that 0 = (log likelihood fn) -  (mleLLvalue - (chi2 / 2))
    }

    # Weibull distribution
    if (distn == "weibull") {
      kappa.hat <- MLEHAT[1, ]
      lambda.hat <- MLEHAT[2, ]
      kappa <- kappa.hat + (d * cos(phi))
      lambda <- lambda.hat + (d * sin(phi))
      temp1 <- sum(x ^ kappa)
      temp2 <- sum(cen * log(x))
      llfn <- r * log(kappa) + (kappa * r * log(lambda)) + ((kappa - 1) * temp2) - (lambda ^ kappa * temp1)
      valuereturn <- llfn - temp   # ID d such that 0 = (log likelihood fn) -  (mleLLvalue - (chi2 / 2))
    }

    # inverse Gaussian distribution
    if(distn == "invgauss"){
      mu.hat <- MLEHAT[1, ]
      lambda.hat <- MLEHAT[2, ]
      mu <- mu.hat + (d * cos(phi))
      lambda <- lambda.hat + (d * sin(phi))
      if (n == r) {     # no censored values; use this formulation because it is more robust
        temp1 <- sum(log(x))
        temp2 <- sum((x - mu) ^ 2 / x)
        # return value for invgauss distn (positive when inside conf region, neg outside)
        # original formulation (without censoring):
        llfn <- (n / 2) * log(lambda) - (n / 2) * log(2 * pi) - (3 / 2) * temp1 - (lambda / (2 * mu ^ 2)) * temp2
      }
      else {            # with censored values
        # note: conf shape parameter (lambda) = (1 / sigma2) from STAR package
        llfn <- sum(cen * log(STAR::dinvgauss(x, mu = mu, sigma2 = 1 / lambda))) +
          sum(as.numeric(cen==0) * log((1 - STAR::pinvgauss(x, mu = mu, sigma2 = 1 / lambda))))
      }
      valuereturn <- llfn - temp   # ID d such that 0 = (log likelihood fn) -  (mleLLvalue - (chi2 / 2))
    }

    # normal distribution
    if (distn == "norm") {
      mu.hat <- MLEHAT[1, ]
      sigma.hat <- MLEHAT[2, ]
      mu <- mu.hat + (d * cos(phi))
      sigma <- sigma.hat + (d * sin(phi))
      if (n == r) {     # no censored values
        temp1 <- sum((x - mu) ^ 2)
        llfn <- (-n / 2) * log(2 * pi) - (n / 2) * log(sigma^2) - (1 / (2*sigma^2)) * temp1
      }
      else {            # no censored values:
        llfn <- sum(cen * log(dnorm(x, mean = mu, sd = sigma))) +
                    sum(as.numeric(cen==0) * log((1 - pnorm(x, mean = mu, sd = sigma))))
      }
      valuereturn <- llfn - temp   # ID d such that 0 = (log likelihood fn) -  (mleLLvalue - (chi2 / 2))
    }

    # lognormal distribution
    if (distn == "lnorm") {
      mu.hat <- MLEHAT[1, ]
      sigma.hat <- MLEHAT[2, ]
      mu <- mu.hat + (d * cos(phi))
      sigma <- sigma.hat + (d * sin(phi))
      temp1 <- sum(log(x))
      temp2 <- sum((log(x) - mu) ^ 2)
      #if (n == r) {     # no censored values
      #  llfn <- (-n / 2) * log(2 * pi) - (n / 2) * log(sigma^2) - temp1 - (1 / (2*sigma^2)) * temp2
      #}
      #else {            # with censored values:
      llfn <- sum(cen * log(dlnorm(x, meanlog = mu, sdlog = sigma))) +
        sum(as.numeric(cen==0) * log((1 - plnorm(x, meanlog = mu, sdlog = sigma))))
      valuereturn <- llfn - temp   # ID d such that 0 = (log likelihood fn) -  (mleLLvalue - (chi2 / 2))
      #}
    }

    # logistic distribution
    if (distn == "logis") {
      mu.hat <- MLEHAT[1, ]
      sigma.hat <- MLEHAT[2, ]
      mu <- mu.hat + (d * cos(phi))
      sigma <- sigma.hat + (d * sin(phi))
      llfn <- sum(cen * log(dlogis(x, location = mu, scale = sigma))) +
        sum(as.numeric(cen==0) * log((1 - plogis(x, location = mu, scale = sigma))))
      valuereturn <- llfn - temp   # ID d such that 0 = (log likelihood fn) -  (mleLLvalue - (chi2 / 2))
    }

    # loglogistic distribution
    if (distn == "llogis") {
      lambda.hat <- MLEHAT[1, ]
      kappa.hat <- MLEHAT[2, ]
      lambda <- lambda.hat + (d * cos(phi))
      kappa <- kappa.hat + (d * sin(phi))
      if (n == r) {     # no censored values; use this formulation because it is more robust
        temp1 <- sum(log(x))
        temp2 <- sum(log(1 + (x / (1 / lambda)) ^ kappa))
        llfn <- n * log(kappa) + n * kappa * log(lambda) + (kappa - 1) * temp1 - 2 * temp2
      }
      else {            # with censored values
        temp1 <- (kappa - 1) * sum(cen * log(lambda * x))
        temp2 <- log(1 + (lambda * x) ^ kappa)
        llfn <- r * (log(lambda) + log(kappa)) + temp1 - sum(cen * temp2) - sum(temp2)
      }
      valuereturn <- llfn - temp   # ID d such that 0 = (log likelihood fn) -  (mleLLvalue - (chi2 / 2))
    }

    # gamma distribution
    if (distn == "gamma") {
      theta.hat <- MLEHAT[1, ]
      kappa.hat <- MLEHAT[2, ]
      theta <- theta.hat + (d * cos(phi))
      kappa <- kappa.hat + (d * sin(phi))
      temp1 <- sum(log(x))
      temp2 <- sum(x / theta)
      if (n == r) {     # no censored values; use this formulation because it is more robust
        llfn <- -n * kappa * log(theta) + (kappa - 1) * temp1 - temp2 - n * log(gamma(kappa))
      }
      else {            # with censored values
        llfn <- sum(cen * log(dgamma(x, scale = theta, shape = kappa))) +
          sum(as.numeric(cen==0) * log((1 - pgamma(x, scale = theta, shape = kappa))))
        #llfn <- sum(cen * log(dgamma(x, scale = theta, shape = kappa) /
        #              (1 - pgamma(x, scale = theta, shape = kappa)))) +
        #              sum(log(1 - pgamma(x, scale = theta, shape = kappa)))
      }
      valuereturn <- llfn - temp   # ID d such that 0 = (log likelihood fn) -  (mleLLvalue - (chi2 / 2))
    }

    # uniform distribution
    if (distn == "unif") {
      a.hat <- MLEHAT[1, ]
      b.hat <- MLEHAT[2, ]
      a <- a.hat + (d * cos(phi))
      b <- b.hat + (d * sin(phi))
      #llfn <- -n * log(b - a)      # original formulation (without censoring)
      llfn <- sum(log(b - x)[cen == 0]) - n * log(b - a)  # with censoring
      valuereturn <- llfn - temp   # ID d such that 0 = (log likelihood fn) -  (mleLLvalue - (chi2 / 2))
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
          # future improvement note: location param distns such as norm & cauchy need variance account for tempUpper estimate; low var & high magnitude location will be problematic
          # if: arbitrary upper bound for 1st quad relative to MLE, or 1st & 2nd quad for distns with (-inf, inf) x-domain
          # else if: distiguish upper bound along y-axis assuming x-domain isn't (-inf, inf)
          # else: upper bound along x-axis
          xinfdomain <- c("cauchy", "norm", "lnorm", "logis", "unif")
          if ((phi[j] <= pi / 2) || ((distn %in% xinfdomain) && (phi[j] <= pi))) {              # phi angles of 0 to pi / 2
            tempUpper <- umult * max(theta1.hat, theta2.hat)
          }
          else if ((phi[j] <= pi +  atan(theta2.hat/theta1.hat))  && !(distn %in% xinfdomain)) {    # phi angles of pi / 2 through intercept with origin
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
          #print(theta2.hat)
          #print(paste0("phi is: ", phi[j]))
          #print(paste0("tempUpper is: ", tempUpper))
          g <- suppressWarnings(try(g1 <- uniroot(llrsolve, lower = 0, upper = abs(tempUpper),
                                                  phi = phi[j], MLE = MLEHAT, mleLLvalue = mleLLvalue, x = samp,
                                                  cen = cen, chi2 = qchisq(1 - alpha, 2), tol = tol), silent = TRUE))
          # errors in uniroot (calculating "g") are possible, and a result of poor upper bound selection
          # the next section of code troubleshoots those errors in two ways:
          # 1. it decreases the upper value by a factor of 1/2 (per iteration), and also
          # 2. it identifies a relatively small upper value (based on min(theta1.hat, theta2.hat)) and then increments it
          # The first option is valuable when the upper limit couldn't be evaluated (-infinity)
          # The second option is valuable when vastly different parameter orders of magnitude push the search too far for phi calcs ~parallel the smaller axis
          z <- 0                       # counter for "try-error" recovery attempts below
          tempUpperLo <- 0             # initialize; will incrementally push upward seeking uniroot upper param
          tempUpperHi <- tempUpper     # initialize; will incrementally push upward seeking uniroot upper param
          #print(phi[j])
          #print(j)
          #if(class(g) == "list") {print(g$root)}
          while ((class(g) == "try-error") && (z < 6) && (tempUpperLo < tempUpperHi)) {
            z <- z + 1
            tempUpperHi <- tempUpper * (0.5 ^ z)
            tempUpperLo <- 5 ^ (z - 1) * min(c(abs(theta1.hat), abs(theta2.hat)))
            #print(paste0("-------------- problematic phi value: ", phi[j]))
            #print(paste0(z, " ****************************"))
            #print(paste0("Error; correction sought with uniroot upper bound modifications... ", tempUpperLo, " and ", tempUpperHi)))
            g <- suppressWarnings(try(g1 <- uniroot(llrsolve, lower = 0, upper = abs(tempUpperHi),
                                     phi = phi[j], MLE = MLEHAT, mleLLvalue = mleLLvalue, x = samp,
                                     cen = cen, chi2 = qchisq(1 - alpha, 2), tol = tol), silent = TRUE))
            if (class(g) == "try-error") {
              #print(paste0("...............try pushing up the min value to ", tempUpperLo))
              g <- suppressWarnings(try(g1 <- uniroot(llrsolve, lower = 0, upper = abs(tempUpperLo),
                           phi = phi[j], MLE = MLEHAT, mleLLvalue = mleLLvalue, x = samp,
                           cen = cen, chi2 = qchisq(1 - alpha, 2), tol = tol), silent = TRUE))
            }
            if (class(g) != "try-error") {
              #print("...*error overcome* through uniroot 'upper' argument modification; algorithm resuming...")
            }
          }
          if (class(g) != "try-error") {
            done <- 1
          }
        }   # end if (done == 0)
      }     # end for (umult in ...)

      if (class(g) == "try-error") {
        print("-------------------------------------------------------------------------------------------------")
        print("R uniroot failure searching for confidence region boundary---challenging parameters and/or shape.")
        print("Unable to produce a confidence region for the given sample and/or parameterization.")
        #print("Uniroot failure.  Challenging parameters and/or shape requires customizing uniroot bounds.")
        print("-------------------------------------------------------------------------------------------------")
      }
      d_phi[j] <- g$root                                                    # saves radial length per phi
      cartesian[j,] <- MLEHAT + d_phi[j] * (c(cos(phi[j]), sin(phi[j])))    # saves CR coordinates per phi
    }       # end for (j in 1:npoints)
    invisible(cartesian)     # return confidence region coordinates w/o displaying to screen
  }         # end crsolve

  ######################################################################################

  crsmooth = function(maxdeg, samp, cen, alpha, mle.list, ellipse_n,
                      xrepair, phinewstore, repairinfo, jumpxy, repairpass,
                      repaircrlist, repairq)  {
    #rm(cr)                               # remove previously stored values
    maxrad <- maxdeg * (pi / 180)         # allowable angle converted to radians
    maxcurrent <- pi                      # initialize as > maxrad to enter while loop
    theta1.hat <-  mle.list$theta1.hat
    theta2.hat <- mle.list$theta2.hat

    # admin note:
    # design to keep CR points when beginning a jump-center heuristic proves slower than without keeping them,
    # therefore those sections of code are bypassed below (numberrepairs < 5) line.  This is attributed at least
    # in-part due to its greater number of points to require angle analysis, etc.  The code supporting it is, however, kept
    # in the event further inquiry into that options is sought in the future.  A better options would be to present
    # an alternative that only looks to verify maxdeg between the points defined by its border region rather than
    # the complete 2pi range around its jump-center.  That enhancement is left for a rainy day...

    # initiallize CR parameters for first-time use
    numberrepairs <- sum(which(is.numeric(repairinfo[1, ] != 0)))
    #if (is.null(repairinfo)) {                         # always keep points when repairing (NOT recommended; takes longer)
    #if (is.null(repairinfo) || numberrepairs < 2) {    # only keep points for 1 repair quadrant
    #if (is.null(repairinfo) || numberrepairs < 5)  {   # never keep points (recommended)
    if(ellipse_n == 4){
      phi <- seq(0, 2*pi, length.out = 5)   # initialize phi search with length.out equidistant points
      phi <- phi[2:length(phi)]             # eliminate redundant point
    }
    else{
      #if !(distn %in% c("list_distns_here_WITHOUT_asesolve_ability")) {    # use if distn arises where ASE not possible
      ase <- asesolve(x = dataset, cen, theta1.hat, theta2.hat)
      phi <- angles(a = ase$theta1.ase, b = ase$theta2.ase, npoints = ellipse_n)
      #}
      #else {                                                                  # estimate aspect ratio from MLE
      #  phi <- angles(a = theta1.hat, b = theta2.hat, npoints = ellipse_n)
      #}
    }
    #}

    # (repair alternatives below did not save time and/or were subject to lack of graph detail & are therefore commented-out)
    # (however, kept for future reference)
    #if (!is.null(repairinfo)) {               # an inaccessible region repair pass; establish parameters to reflect current CR plot
    #  repairLindex <- which(repaircrlist$x == repairinfo[2, repairq])  # left border
    #  repairRindex <- which(repaircrlist$x == repairinfo[5, repairq])  # right border
    #  checkrange <- min(c(repairLindex, repairRindex)):max(c(repairLindex, repairRindex))
    #  maxdiff <- which(abs(diff(repaircrlist$x[checkrange])) == max(abs(diff(repaircrlist$x[checkrange]))))
    #  repairborder1index <- min(c(repairLindex, repairRindex)) + maxdiff - 1  # 1st phi angle index
    #  # identify phi value from jump-center to the two points bordering the inaccessible region
    #  repairphi <- atan((repaircrlist$y - as.numeric(repairinfo[9, repairq])) /
    #              (repaircrlist$x - as.numeric(repairinfo[8, repairq])))
    #  # break into cases where jump-center x-value < point x-value, and > point x-value to ensure 0 < phi < 2pi appropriately
    #  # if jump-center x-value is < border x-value correct negative phi values to become 2pi - phi
    #  # if jump-center x-value is > border x-value correct value to pi + phi
    #  indexRjumpc <- which(repaircrlist$x > as.numeric(repairinfo[8, repairq]))         # index of points to right of jump-center
    #  indexLjumpc <- which(repaircrlist$x < as.numeric(repairinfo[8, repairq]))         # index of points to left of jump-center
    #  repairphi[indexRjumpc] <- (sign(repairphi[indexRjumpc]) - 1) * (-pi) + repairphi[indexRjumpc]
    #  repairphi[indexLjumpc] <- pi + repairphi[indexLjumpc]
    #  #phi <- repairphi
    #  #cr <- crsolve(samp, cen, alpha = alpha, mle.list = mle.list, phi)      # confidence region points w initial phi's
    #
    #  # ID phi values from jump-center to borders of inaccessible region being assessed so that can terminate
    #  # while loop when phinew additions are no longer made in that region (ignoring issues elsewhere)
    #  repairborderphis <- repairphi[repairborder1index:(repairborder1index + 1)]
    #  #phi <- sort(c(phi, repairborderphis))
    #
    #  # note: the approach below does not work because it re-orders sequence each iteration according to phi values, which makes
    #  # the multi-point roots placed in in-correct sequence ("hopping" back and forth across CR)
    #  #cr <- matrix(c(repaircrlist$x, repaircrlist$y), nrow = length(repaircrlist$x), ncol = 2)
    #}
    cr <- crsolve(samp, cen, alpha = alpha, mle.list = mle.list, phi)      # confidence region points w initial phi's

    transf.scale <- 1                     # initialize transformation scaling factor at 1
    count <- 0                            # administrative counter for while loop
    d01 <- 0
    d02 <- 0
    d12 <- 0
    while (maxcurrent > maxrad) {                                          # in radians
      count <- count + 1
      if (count != 1){
        index_phinew <- match(phinew, phi)
        index_phiold <- match(phiold, phi)
        crnew <- crsolve(samp, cen, alpha = alpha, mle.list = mle.list, phi = phinew)     # uniroot calcs for CR points of new phis
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
      # (below commented-out section is problematic for implimentation but kept here for possible future attempts)
      # Its intent was to cut off search when phinew points outside the region of interest were all that remain
      # ...doing so, however, caused pre-mature termination in some circumstance (i.e. crplot(c(1.9, 2), 0.05, "invgauss"))
      #if (!is.null(repairinfo) && (count > 10)) {
      #  phinew <- phinew[intersect(which(phinew >= min(repairborderphis[1], repairborderphis[2])),
      #                             which(phinew <= max(repairborderphis[1], repairborderphis[2])))]
      #  if (length(phinew) == 0) {
      #    count <- 50  # force exit of while loop because phinew additions no longer apply to region between borderphis
      #  }
      #}
      if (sum(phinew) != 0) {phi <- sort(unique(c(phi, phinew))) }
      if (count >= 50) {                                    # indicates there is an inaccessible region

        ############################################################################################
        # May 2018 addition enables charting previously inaccessible regions by identifying
        # jump-center locations away from the MLE but within the CR where a similar smoothing
        # algorithm is done using the radial log likelihood approach.  Recursively calls crsmooth.
        ############################################################################################
        if ((repairpass != TRUE) && (repair == TRUE)) {

        #jumpshift <- 0.5                   # % along available shift range for phi to locate jump point
        #jumpuphill <- min(alpha, 0.01)     # percentage uphill from CR boarder

        xrepair <- xrepair + 1
        # identify all areas needing repair and their points to facilitate reconstruction
        # this is done during the first pass (xrepair = 1) only, since these variables change as areas are repaired
        if (xrepair == 1) {
          #print("Creating and storing repair info")
          phinewstore <- phinew
          leftphi <- leftx <- lefty <- rep(0, 4)
          rightphi <- rightx <- righty <- rep(0, 4)
          jumpphi <- gaptype <- phi1 <- phi2 <- rep(0, 4)
          jumpxy <- matrix(rep(0, 8), ncol = 2)

          # ensure a user-entered jumpuphill value is not infeasible
          if ((alpha + jumpuphill) > 1) {
            warning("infeasible jumpuphill value ((alpha + jumpuphill) > 1); restoring default of min(alpha, 0.01)")
            jumpuphill <- min(alpha, 0.01)
            #jumpuphill <- 0.9 * alpha + 0.1     # (alpha + 0.1(1 - alpha))
          }

          #####################################################
          # Quadrant I (with respect to MLE) repairs
          if (any(phinewstore < (pi / 2))) {
            #print("Creating jump-point to access and patch unreachable area (in Quad I relative to MLE)")
            print("The 1st quadrant (with respect to the MLE) currently does not offer repairs (no known cases with issues here; yours may be the first!)")
            q <- 1      # quad I
          }

          #####################################################
          # Quadrant II (with respect to MLE) identification
          if (length(phinewstore[(phinewstore > pi / 2) & (phinewstore < pi)]) > 0) {
            #print("Creating jump-point to access and patch unreachable area (in Quad II relative to MLE)")
            q <- 2      # quad II

            # identify the boarders of the unreachable area
            # offset two points above and below the inaccessible region to ensure it is partitioned with one reference point on
            # each side (o.w. rare roundoff issues may result with both co-located)
            quadindex <- intersect(which(phinewstore > (pi / 2)), which(phinewstore < pi))    # unique to this quadrant
            leftphi[q] <- max(phinewstore[quadindex])                                            # unique to this quadrant
            leftx[q] <- cr[min(which(crlist$phi > leftphi[q])) + 1, 1]
            lefty[q] <- cr[min(which(crlist$phi > leftphi[q])) + 1, 2]
            rightphi[q] <- min(phinewstore[quadindex])                                           # unique to this quadrant
            rightx[q] <- cr[max(which(crlist$phi < rightphi[q])) - 1, 1]
            righty[q] <- cr[max(which(crlist$phi < rightphi[q])) - 1, 2]
            #points(leftx[q], lefty[q], pch = 16, col = "green")
            #points(rightx[q], righty[q], pch = 16, col = "blue")
            # create a "new" center (instead of MLE) to access previously unreachable CR points as follows:
            # assess shift range available in x-direction within CR boundary,
            # determine MLE phi corresponding to a % along that available shift range
            # offset slightly inside CR via "uphill" alpha adjustment to locate jump x-y coordinates
            shiftmaxindex <- min(intersect(which(cr[, 1] < rightx[q]), which(cr[, 2] < righty[q])))         # unique to this quadrant
            shift <- jumpshift * (rightx[q] - cr[shiftmaxindex, 1])                                      # unique to this quadrant
            jumpphi[q] <- pi / 2 + atan((theta1.hat - (rightx[q] - shift)) / (righty[q] - theta2.hat))      # unique to this quadrant
            jumpxy[q,] <- crsolve(samp, cen, alpha = alpha + jumpuphill, mle.list = mle.list, phi = jumpphi[q])
            #points(jumpxy[q, 1], jumpxy[q, 2], pch = 16, col = "red")
            phi1[q] <- pi / 2 + atan((jumpxy[q,1] - leftx[q]) / (lefty[q] - jumpxy[q,2]))      # region of interest: (0, phi1)
            phi2[q] <- 2 * pi - atan((jumpxy[q,2] - righty[q]) / (rightx[q] - jumpxy[q,1]))    # region of interest: (phi2, 2pi]
          }

          #####################################################
          # Quadrant III (with respect to MLE) identification
          if (length(phinewstore[(phinewstore > pi) & (phinewstore < 3 * pi / 2)]) > 0) {
            #print("Creating jump-point to access and patch unreachable area (in Quad III relative to MLE)")
            q <- 3      # quad III

            # identify the borders of the unreachable area
            quadindex <- intersect(which(phinewstore > pi), which(phinewstore < 3 * pi / 2))  # unique to this quadrant
            leftphi[q] <- max(phinewstore[quadindex])                                            # unique to this quadrant
            leftx[q] <- cr[min(which(crlist$phi > leftphi[q])) + 1, 1]
            lefty[q] <- cr[min(which(crlist$phi > leftphi[q])) + 1, 2]
            rightphi[q] <- min(phinewstore[quadindex])                                           # unique to this quadrant
            rightx[q] <- cr[max(which(crlist$phi < rightphi[q])) - 1, 1]
            righty[q] <- cr[max(which(crlist$phi < rightphi[q])) - 1, 2]
            #points(leftx[q], lefty[q], pch = 16, col = "green")
            #points(rightx[q], righty[q], pch = 16, col = "blue")
            # create a "new" center (instead of MLE) to access previously unreachable CR points as follows:
            # assess shift range available in y-direction within CR boundary,
            # determine MLE phi corresponding to a % along that available shift range
            # offset slightly inside CR via "uphill" alpha adjustment to locate jump x-y coordinates
            # if statement first identifies if uncharted region is above or below border of unreachable are
            if (cr[max(which(crlist$phi < rightphi[q])) - 2, 2] > cr[min(which(crlist$phi > leftphi[q])) + 2, 2]) {
              gaptype[q] <- "|"   # vertical gap segment that jumpphi will cross
              shiftmaxindex <- min(intersect(which(cr[, 1] > rightx[q]), which(cr[, 2] < righty[q])))      # unique to this quadrant
              shift <- jumpshift * (righty[q] - cr[shiftmaxindex, 2])                                   # unique to this quadrant
              jumpphi[q] <- pi + atan((theta2.hat - righty[q] + shift) / (theta1.hat - rightx[q]))         # unique to this quadrant
            }
            else if (cr[max(which(crlist$phi < rightphi[q])) - 2, 2] < cr[min(which(crlist$phi > leftphi[q])) + 2, 2]) {
              gaptype[q] <- "-"
              shiftmaxindex <- max(intersect(which(cr[, 1] < leftx[q]), which(cr[, 2] > lefty[q])))      # unique to this quadrant
              shift <- jumpshift * (leftx[q] - cr[shiftmaxindex, 1])                                   # unique to this quadrant
              jumpphi[q] <- pi + atan((theta2.hat - lefty[q]) / (theta1.hat - leftx[q] + shift))         # unique to this quadrant
            }
            jumpxy[q,] <- crsolve(samp, cen, alpha = alpha + jumpuphill, mle.list = mle.list, phi = jumpphi[q])
            #points(cr[shiftmaxindex, 1], cr[shiftmaxindex, 2], pch = 16, col = "brown")
            #points(jumpxy[q, 1], jumpxy[q, 2], pch = 16, col = "red")
            phi1[q] <- pi + atan((jumpxy[q, 2] - lefty[q]) / (jumpxy[q, 1] - leftx[q]))      # leftmost phi relative to jump-center
            phi2[q] <- atan((righty[q] - jumpxy[q, 2]) / (rightx[q] - jumpxy[q, 1]))         # rightmost phi relative to jump-center
          }

          #####################################################
          # Quadrant IV (with respect to MLE) repairs
          if (any(phinewstore > (3 * pi / 2))) {
            #print("Creating jump-point to access and patch unreachable area (in Quad IV relative to MLE)")
            q <- 4               # quad IV

            # identify the boarders of the unreachable area
            quadindex <- which(phinewstore > (3 * pi / 2))                               # unique to this quadrant
            leftphi[q] <- max(phinew[quadindex])                                            # unique to this quadrant; phi angle for left point (w.r.t. MLE prespective)
            leftx[q] <- cr[min(which(crlist$phi > leftphi[q])) + 1, 1]
            lefty[q] <- cr[min(which(crlist$phi > leftphi[q])) + 1, 2]
            rightphi[q] <- min(phinew[quadindex])                                           # unique to this quadrant
            rightx[q] <- cr[max(which(crlist$phi < rightphi[q])) - 1, 1]
            righty[q] <- cr[max(which(crlist$phi < rightphi[q])) - 1, 2]
            #points(leftx[q], lefty[q], pch = 16, col = "green")
            #points(rightx[q], righty[q], pch = 16, col = "blue")

            # create a "new" center (instead of MLE) to access previously unreachable CR points as follows:
            # assess shift range available in y-direction within CR boundary,
            # determine MLE phi corresponding to a % along that available shift range
            # offset slightly inside CR via "uphill" alpha adjustment to locate jump x-y coordinates
            shiftmaxindex <- max(intersect(which(cr[, 1] < leftx[q]), which(cr[, 2] < lefty[q])))      # unique to this quadrant
            #points(cr[shiftmaxindex, 1], cr[shiftmaxindex, 2], pch = 16, col = "purple")
            shift <- jumpshift * (lefty[q] - cr[shiftmaxindex, 2])
            jumpphi[q] <- 2 * pi - atan((theta2.hat - lefty[q] + shift) / (leftx[q] - theta1.hat))      # unique to this quadrant
            jumpxy[q,] <- crsolve(samp, cen, alpha = alpha + jumpuphill, mle.list = mle.list, phi = jumpphi[q])
            #points(jumpxy[q, 1], jumpxy[q, 2], pch = 16, col = "orange")
            phi1[q] <- pi / 2 + atan((jumpxy[q,1] - leftx[q]) / (lefty[q] - jumpxy[q,2]))      # region of interest: (0, phi1)
            phi2[q] <- 2 * pi - atan((jumpxy[q,2] - righty[q]) / (rightx[q] - jumpxy[q,1]))    # region of interest: (phi2, 2pi]
          }

        # repair info must be kept un-impacted by upcoming recursive crsmooth calls; data stored accordingly here
        done <- rep(FALSE, 4)
        jumpx <- jumpxy[, 1]
        jumpy <- jumpxy[, 2]
        repairinfo <- rbind(leftphi, leftx, lefty, rightphi, rightx, righty, jumpphi, jumpx, jumpy, phi1, phi2, done, gaptype)
        #print(repairinfo)
        #print(jumpxy)
        #stop()
        warning("alternate-centerpoint(s) used to repair plot regions inaccessible via a radial angle from its MLE")
        #message("alternate-centerpoint(s) used to repair plot regions inaccessible via a radial angle from its MLE")
        }                                    # end if xrepair == 1
        else {                               # all iterations following 1st pass must re-assume repair parameters
          # each column in repairinfo represents a quadrant (with respect to the MLE)
          # "left" items are indicative of the left point from the perspective of the MLE (higher phi value); right is lower phi value
          #print("Loading repair info")
          leftphi <- repairinfo[1,]        # phi value to left point
          leftx <- repairinfo[2,]       # x value of left point
          lefty <- repairinfo[3,]       # y value of left point
          rightphi <- repairinfo[4,]
          rightx <- repairinfo[5,]
          righty <- repairinfo[6,]
          jumpphi <- repairinfo[7,]     # angle from MLE to locate jump-center
          jumpx <- repairinfo[8,]       # x value of jump center
          jumpy <- repairinfo[9,]       # y value of jump center
          phi1 <- repairinfo[10,]        # angle from jump-center to inaccessible region border point 1 (smaller phi angle)
          phi2 <- repairinfo[11,]        # angle from jump-center to inaccessible region border point 2 (larger phi angle)
          done <- repairinfo[12,]       # will record after quadrant repairs are done
          gaptype <- repairinfo[13,]    # identifies if "-" or "|" type repair needed (gap segment that jumpphi angle crosses)
        }

        ###############################################################################
        # half-way through repairs; relevant points and angles recorded above.
        # next, recursively call crsmooth on jump-centers to populate missing regions.
        ###############################################################################

        #####################################################
        # Quadrant I (with respect to MLE) repairs
        q <- 1
        if ((length(phinewstore[phinewstore < pi / 2] > 0)) && (done[q] != TRUE)) {
          #print("Repairing (Quad I relative to MLE)")
          message("quad I (relative to the MLE) repair not present because no confirmed cases yet (this might be the first!)")
          repairinfo[12, q] <- TRUE      # annotate this quad is done
        }

        #####################################################
        # Quadrant II (with respect to MLE) repairs
        q <- 2        # quad II
        if ((length(phinewstore[(phinewstore > pi / 2) & (phinewstore < pi)]) > 0) && (done[q] != TRUE)) {
          #message("...repairing (Quad II relative to MLE)")
          repairinfo[12, q] <- TRUE      # annotate this quad is done
          gaptype[q] <- "-"              # assumed b/c no "|" cases currently known, requires modification similar to Quad III otherwise

          # store the "new" centerpoint as theta1.hat & thata2.hat (albeit not an MLE; stays consistent with existing code)
          # maintain the MLE log-likelihood value as mleLLvalues because it is still needed in llrsolve
          # run smoothing search algorithm from jump point, and trim results to region of interest (per phi1, phi2 values)
          # note: since multiple roots available per phi, also check that y-values are not outside region of interest
          jumpinfo <- list("theta1.hat" = jumpxy[q, 1], "theta2.hat" = jumpxy[q, 2], "mleLLvalue" = mle.list$mleLLvalue)
          jpoints <- crsmooth(maxdeg = maxdeg, samp = mydata, cen = cen, alpha = alpha, mle.list = jumpinfo, ellipse_n = ellipse_n,
                              xrepair = xrepair, phinewstore = phinewstore, repairinfo = repairinfo, jumpxy = jumpxy, repairpass = TRUE,
                              repaircrlist = crlist, repairq = q)

          # identify boundary points produced using jump-center that are relevant to keep (occur in the previously inaccessible region)
          keepindex <- which((jpoints$phi < phi1[q] | jpoints$phi > phi2[q]) & jpoints$y > righty[q])
          jkeep <- list("x" = jpoints$x[keepindex],
                        "y" = jpoints$y[keepindex],
                        "phi" = jpoints$phi[keepindex] )
                        #"y" = jpoints$y[(which(jpoints$phi < phi1[q] | jpoints$phi > phi2[q] & jpoints$y > righty[q]))],
                        #"phi" = jpoints$phi[(which(jpoints$phi < phi1[q] | jpoints$phi > phi2[q] & jpoints$y > righty[q]))])
          #points(jkeep$x, jkeep$y, pch = 16, col = "blue", cex = 0.5)

          # identify phi angle corresponding to MLE to remain consistent
          # to get phi value in range [0, 2pi) need to adjust off tan which computes from -pi/2 to pi/2:
          # if (jpoints$y - theta2.hat) is positive, add pi/2
          # if (jpoints$y - theta2.hat) is negative, add 3pi/2
          phiactual <- (pi - (pi / 2) * sign(jkeep$y - theta2.hat)) +
            atan(-(jkeep$x - theta1.hat) / (jkeep$y - theta2.hat))

          # insert additional CR boundary points into existing list
          # sequence angles being kept in the proper order
          go <- c(which(jkeep$phi > phi1[q]), which(jkeep$phi <= phi1[q]))
          # ensure new points integrate into the combined confidence region plot at right location
          # these steps are unique to this quadrant
          options <- which(crlist$y > jumpxy[q, 2])                         # candidates for insertion after are above the jump-center
          insertafter <- which(crlist$phi == min(crlist$phi[options])) - 1  # "new" points fall before the lowest phi value among those options
          crlist$x <- append(crlist$x, jkeep$x[go], after = insertafter)
          crlist$y <- append(crlist$y, jkeep$y[go], after = insertafter)
          crlist$phi <- append(crlist$phi, phiactual[go], after = insertafter)

        }

        #####################################################
        # Quadrant III (with respect to MLE) repairs
        q <- 3
        if ((length(phinewstore[(phinewstore > pi) & (phinewstore < 3 * pi / 2)]) > 0) && (done[q] != TRUE)) {
          #message("...repairing (Quad III relative to MLE)")
          repairinfo[12, q] <- TRUE      # annotate this quad is done

          # store the "new" centerpoint as theta1.hat & thata2.hat (albeit not an MLE; stays consistent with existing code)
          # maintain the MLE log-likelihood value as mleLLvalues because it is still needed in llrsolve
          # run smoothing search algorithm from jump point, and trim results to region of interest (per phi1, phi2 values)
          # note: since multiple roots available per phi, also check that y-values are not outside region of interest
          jumpinfo <- list("theta1.hat" = jumpxy[q, 1], "theta2.hat" = jumpxy[q, 2], "mleLLvalue" = mle.list$mleLLvalue)
          jpoints <- crsmooth(maxdeg = maxdeg, samp = mydata, cen = cen, alpha = alpha, mle.list = jumpinfo, ellipse_n = ellipse_n,
                              xrepair = xrepair, phinewstore = phinewstore, repairinfo = repairinfo, jumpxy = jumpxy, repairpass = TRUE,
                              repaircrlist = crlist, repairq = q)

          # identify phi angle corresponding to MLE to remain consistent
          # to get phi value in range [0, 2pi) need to adjust off tan which computes from -pi/2 to pi/2:
          # if (jpoints$y - theta2.hat) is positive, add pi/2
          # if (jpoints$y - theta2.hat) is negative, add 3pi/2
          phiactual <- (pi - (pi / 2) * sign(jpoints$y - theta2.hat)) +
            atan(-(jpoints$x - theta1.hat) / (jpoints$y - theta2.hat))

          # identify boundary points produced using jump-center that are relevant to keep (occur in the previously inaccessible region)
          # conditional on what side of the inaccessible region border the uncharted region lies
          if (gaptype[q] == "|") {
            keepindex <- which((jpoints$phi < phi1[q]) & (jpoints$x < rightx[q]) & (phiactual < rightphi[q]))
          }
          else if (gaptype[q] == "-") {
            keepindex <- which((jpoints$phi > phi1[q]) & (jpoints$x < leftx[q]) & (phiactual > rightphi[q]))
          }
          jkeep <- list("x" = jpoints$x[keepindex],
                        "y" = jpoints$y[keepindex],
                        "phi" = jpoints$phi[keepindex] )
          #points(jkeep$x, jkeep$y, pch = 16, col = "yellow", cex = 0.5)

          # update phiactual angles w.r.t. MLE to represent only jump-center points kept (jkeep variable)
          phiactual <- (pi - (pi / 2) * sign(jkeep$y - theta2.hat)) +
            atan(-(jkeep$x - theta1.hat) / (jkeep$y - theta2.hat))

          # insert additional CR boundary points into existing list
          # sequence angles being kept in the proper order
          go <- c(which(jkeep$phi > phi1[q]), which(jkeep$phi <= phi1[q]))
          # ensure new points integrate into the combined confidence region plot at right location
          # these steps are unique to this quadrant
          if (gaptype[q] == "|") {
            options <- which(crlist$x < jumpxy[q, 1])                         # candidates for insertion after are left of the jump-center
            insertafter <- which(crlist$phi == min(crlist$phi[options])) - 1  # "new" points fall before the lowest phi value among those options
          }
          else if (gaptype[q] == "-") {
            options <- which(crlist$y < jumpxy[q, 2])                         # candidates for insertion after are below the jump-center
            insertafter <- which(crlist$phi == max(crlist$phi[options]))      # "new" points fall after the highest phi value among those options
          }
          crlist$x <- append(crlist$x, jkeep$x[go], after = insertafter)
          crlist$y <- append(crlist$y, jkeep$y[go], after = insertafter)
          crlist$phi <- append(crlist$phi, phiactual[go], after = insertafter)

        }

        #####################################################
        # Quadrant IV (with respect to MLE) repairs
        q <- 4
        if (any(phinewstore > (3 * pi / 2)) && (done[q] != TRUE)) {
          #message("...repairing (Quad IV relative to MLE)")
          repairinfo[12, q] <- TRUE      # annotate this quad is done
          gaptype[q] <- "|"              # assumed b/c no "|" cases currently known, requires modification similar to Quad III otherwise

          # store the "new" centerpoint as theta1.hat & thata2.hat (albeit not an MLE; stays consistent with existing code)
          # maintain the MLE log-likelihood value as mleLLvalues because it is still needed in llrsolve
          # run smoothing search algorithm from jump point, and trim results to region of interest (per phi1, phi2 values)
          # note: since multiple roots available per phi, also check that y-values are not outside region of interest
          jumpinfo <- list("theta1.hat" = jumpxy[q, 1], "theta2.hat" = jumpxy[q, 2], "mleLLvalue" = mle.list$mleLLvalue)
          jpoints <- crsmooth(maxdeg = maxdeg, samp = mydata, cen = cen, alpha = alpha, mle.list = jumpinfo, ellipse_n = ellipse_n,
                              xrepair = xrepair, phinewstore = phinewstore, repairinfo = repairinfo, jumpxy = jumpxy, repairpass = TRUE,
                              repaircrlist = crlist, repairq = q)

          # identify boundary points produced using jump-center that are relevant to keep (occur in the previously inaccessible region)
          keepindex <- which((jpoints$phi < phi1[q] | jpoints$phi > phi2[q]) & jpoints$x > leftx[q])
          jkeep <- list("x" = jpoints$x[keepindex],
                        "y" = jpoints$y[keepindex],
                        "phi" = jpoints$phi[keepindex])
          #points(jkeep$x, jkeep$y, pch = 16, col = "yellow", cex = 0.5)

          # identify phi angle corresponding to MLE to remain consistent
          # to get phi value in range [0, 2pi) need to adjust off tan which computes from -pi/2 to pi/2:
          # if (jpoints$y - theta2.hat) is positive, add pi/2
          # if (jpoints$y - theta2.hat) is negative, add 3pi/2
          phiactual <- (pi - (pi / 2) * sign(jkeep$y - theta2.hat)) +
            atan(-(jkeep$x - theta1.hat) / (jkeep$y - theta2.hat))

          # insert additional CR boundary points into existing list
          # sequence angles being kept in the proper order
          go <- c(which(jkeep$phi > phi1[q]), which(jkeep$phi <= phi1[q]))
          # ensure new points integrate into the combined confidence region plot at right location
          # these steps are unique to this quadrant
          options <- which(crlist$x > jumpxy[q, 1])                         # candidates for insertion after are right of jump-center
          insertafter <- which(crlist$phi == max(crlist$phi[options]))      # "new" points fall after the highest phi value among those options
          crlist$x <- append(crlist$x, jkeep$x[go], after = insertafter)
          crlist$y <- append(crlist$y, jkeep$y[go], after = insertafter)
          crlist$phi <- append(crlist$phi, phiactual[go], after = insertafter)

        }
        #print(jumpxy)   # print jump-center locations corresponding to each respective quadrant (quadrants are relative to MLE)
      }                  # end repairpass
      maxcurrent <- maxrad
      if (repair == FALSE) {
        warning("max iteration tolerance hit; working confidence region plot is shown, however, maxdeg constraint is not met")
        #print("WARNING: max iteration tolerance hit; maxdeg constraint not met, however, working plot is shown.")
        #print("...this is attributable to either regions inaccessible via a radial azimuth from its MLE, or a uniroot tol argument too small.")
      }
    }                # end (if (count > 50) && repairpass != TRUE)
    }                # end (while (maxcurrent > maxrad))

    invisible(crlist)
  }


  ############################################################################################
  ### end of functions; function calls below determine and then plot the confidence region ###
  ############################################################################################

  #ballbearing <- c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.48, 51.84, 51.96, 54.12, 55.56,
  #  67.80, 68.64, 68.64, 68.88, 84.12, 93.12, 98.64, 105.12, 105.84, 127.92, 128.04, 173.40)
  #cen <- c(rep(1, length(dataset)))  # delete after censering successfully incorporated
  mydata <- dataset
  mle.list <- mlesolve(mydata, cen = cen)
  theta1.hat <- mle.list$theta1.hat
  theta2.hat <- mle.list$theta2.hat
  mleLLvalue <- mle.list$mleLLvalue

  # identify confidence region boundary coordinates
  if (distn == "unif") {
      phi <- c(pi / 2, pi)
      cr <- crsolve(samp = dataset, cen = cen, alpha = alpha, mle.list = mle.list, phi = phi)
      crlist <- list("x" = append(cr[, 1], theta1.hat), "y" = append(cr[, 2], theta2.hat), "phi" = append(phi, 0))
      message("uniform distribution confidence region is built using the three points triangulating that region")
  }
  else if (heuristic == 1) {
    crlist <- crsmooth(maxdeg = maxdeg, samp = mydata, cen = cen, alpha = alpha, mle.list = mle.list, ellipse_n = ellipse_n,
                       xrepair = 0, phinewstore = NULL, repairinfo = NULL, jumpxy = NULL,
                       repairpass = FALSE, repaircrlist = NULL, repairq = NULL)
  }
  else if (heuristic == 0) {
    #if !(distn %in% c("list_distns_here_WITHOUT_asesolve_ability")) {    # use if distn arises where ASE not possible
    ase <- asesolve(x = dataset, cen, theta1.hat, theta2.hat)
    phi <- angles(a = ase$theta1.ase, b = ase$theta2.ase, npoints = ellipse_n)
    #}
    #else {                                                                  # estimate aspect ratio from MLE
    #  phi <- angles(a = theta1.hat, b = theta2.hat, npoints = ellipse_n)
    #}
    cr <- crsolve(samp = dataset, cen = cen, alpha = alpha, mle.list = mle.list, phi = phi)
    crlist <- list("x" = cr[, 1], "y" = cr[, 2], "phi" = phi)
  }

  # assemble conf region points and "close" the gap between first & last points graphically
  cr <- cbind(crlist$x, crlist$y)
  cr <- rbind(cr, cr[1, ])

  if (xyswap == TRUE) {
    # identify phi angle corresponding to MLE to remain consistent
    # "reflect" the phi angle across the x = y diagonal should hold true for axes swap
    # assess and index the regions, then perform the transformation:
    # region 1: [0 to 3pi/4]     then phiswap <- pi/2 - phi
    # region 2: (3pi/4, 7pi/4]   then phiswap <- 7pi/2 - phi
    # region 3: (7pi/4, 2pi)     then phiswap <- 5pi/2 - phi
    case1index <- which(crlist$phi <= 3 * pi / 4)   # 1st region
    case2index <- intersect(which(crlist$phi > 3 * pi / 4), which(crlist$phi <= 7 * pi / 4))   # 2nd region
    case3index <- which(crlist$phi > 7 * pi / 4)    # 3rd region
    cswap <- rep(0, length(crlist$phi))
    cswap[case1index] <- pi / 2
    cswap[case2index] <- 10 * pi / 4
    cswap[case3index] <- 5 * pi / 2
    phiswap <- cswap - crlist$phi
    crlist$phi <- phiswap
  }
  phi <- crlist$phi


  if (showplot == TRUE) {
  # construct the confidence region plot
  # override default plot characteristics if user specified
  # i.e. mar, xlim, ylim, main, xlab, ylab, etc.

  # margin adjustments
  par(xpd = FALSE)
  if (!is.null(mar)) {
    par(mar = mar)
  }

  # label axes appropriately in the absence of a user specified label
  disttype <- c("gamma", "invgauss", "llogis", "lnorm", "norm", "unif", "weibull", "cauchy", "logis")
  xaxislabel <- c(expression(theta), expression(mu), expression(lambda), expression(mu),
                  expression(mu), expression(a), expression(kappa), expression(a), expression(mu))
  yaxislabel <- c(expression(kappa), expression(lambda), expression(kappa), expression(sigma),
                  expression(sigma), expression(b), expression(lambda), expression(s), expression(sigma))
  if (!is.expression(xlab)) {
    if (xlab == "") {
      if (xyswap == FALSE) {
        xlab = xaxislabel[which(disttype == distn)]
      }
      else if (xyswap == TRUE) {
        xlab = yaxislabel[which(disttype == distn)]
      }
    }
  }
  if (!is.expression(ylab)) {
    if (ylab == "") {
      if (xyswap == FALSE) {
        ylab = yaxislabel[which(disttype == distn)]
      }
      else if (xyswap == TRUE) {
        ylab = xaxislabel[which(disttype == distn)]
      }
    }
  }

  # axis tick-marks based on mlelab and origin
  xlabs <- c(min(cr[, 1]), max(cr[, 1]))
  ylabs <- c(min(cr[, 2]), max(cr[, 2]))
  if (mlelab == TRUE) {
    xlabs <- c(xlabs, theta1.hat)
    ylabs <- c(ylabs, theta2.hat)
  }
  if (origin == TRUE) {
    if (min(cr[,1]) > 0) {
      xmin <- 0
    }
    else {
      xmin <- min(cr[,1])
    }
    if (min(cr[,2]) > 0) {
      ymin <- 0
    }
    else {
      ymin <- min(cr[,2])
    }
    xlabs <- sort(c(xmin, xlabs))
    ylabs <- sort(c(ymin, ylabs))
  }
  else {
    xmin <- min(cr[,1])
    ymin <- min(cr[,2])
  }

  # axis limits and tick-marks if xlim and/or ylim entries:
  if(is.null(xlim)) {
    if (xyswap == FALSE) {
      xlim <- c(xmin, max(cr[,1]))
    }
    else if (xyswap == TRUE) {
      xlim <- c(ymin, max(cr[,2]))
    }
  }
  else {
    xlabs <- c(xlim, xlabs)
    xlabs <- xlabs[xlabs >= xlim[1] & xlabs <= xlim[2]]
  }
  if(is.null(ylim)) {
    if (xyswap == FALSE) {
      ylim <- c(ymin, max(cr[,2]))
    }
    else if (xyswap == TRUE) {
      ylim <- c(xmin, max(cr[,1]))
    }
  }
  else {
    ylabs <- c(ylim, ylabs)
    ylabs <- ylabs[ylabs >= ylim[1] & ylabs <= ylim[2]]
  }

  # plot
    if (xyswap == FALSE) {
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
    }
    else if (xyswap == TRUE) {
      crswap <- cbind(cr[,2], cr[,1])
      plot(crswap, xlab = xlab, ylab = ylab, main = main, ylim = ylim, xlim = xlim,
           axes = FALSE, type = 'l')
      my.atx <- unique(round(ylabs, sf[1]))
      my.aty <- unique(round(xlabs, sf[2]))
      axis(side = 1, at = my.atx, las = xlas)
      axis(side = 2, at = my.aty, las = ylas)
      if (mlelab == TRUE) {
        points(theta2.hat, theta1.hat, pch = 3)
      }
      if (pts == TRUE) points(unique(crswap), lwd = 0.65)
    }
    # enable plotting beyond margins for any post-processing user add-ons
    # ...or disable
    #par(xpd = TRUE)
    par(xpd = FALSE)
  }

  #print("dataset was: ")
  #print(dataset)
  # returned values (if requested), otherwise only output to screen # boundary points.
  if ((length(dataset) <= 5) && (distn != "unif")) {
    warning("small sample size is ill-suited to invoke the asymptotic properties assumed by the confidence region plot")
    #print("WARNING: small sample sizes can yield irregular confidence region shapes that are unatainable using crplot.")
    #print("WARNING: small sample sizes violate the asymptotic properties assumption used to produce the confidence region.")
  }
  if (info == TRUE) {
    print(paste0("Confidence region plot complete; made using ", length(phi)," boundary points."))
    if (distn == "weibull") {
      #print(paste0("MLE value is: (kappa.hat = ", theta1.hat, ", lambda.hat = ", theta2.hat,")"))
      crlist <- list("kappa" = crlist$x, "lambda" = crlist$y, "phi" = crlist$phi,
                     "kappahat" = mle.list$theta1.hat, "lambdahat" = mle.list$theta2.hat)
    }
    else if (distn == "invgauss") {
      #print(paste0("MLE value is: (mu.hat = ", theta1.hat, ", lambda.hat = ", theta2.hat,")"))
      crlist <- list("mu" = crlist$x, "lambda" = crlist$y, "phi" = crlist$phi,
                     "muhat" = mle.list$theta1.hat, "lambdahat" = mle.list$theta2.hat)
    }
    else if (distn == "norm") {
      #print(paste0("MLE value is: (mu.hat = ", theta1.hat, ", sigma.hat = ", theta2.hat,")"))
      crlist <- list("mu" = crlist$x, "sigma" = crlist$y, "phi" = crlist$phi,
                     "muhat" = mle.list$theta1.hat, "sigmahat" = mle.list$theta2.hat)
    }
    else if (distn == "lnorm") {
      #print(paste0("MLE value is: (mu.hat = ", theta1.hat, ", sigma.hat = ", theta2.hat,")"))
      crlist <- list("mu" = crlist$x, "sigma" = crlist$y, "phi" = crlist$phi,
                     "muhat" = mle.list$theta1.hat, "sigmahat" = mle.list$theta2.hat)
    }
    else if (distn == "logis") {
      #print(paste0("MLE value is: (mu.hat = ", theta1.hat, ", sigma.hat = ", theta2.hat,")"))
      crlist <- list("mu" = crlist$x, "sigma" = crlist$y, "phi" = crlist$phi,
                     "muhat" = mle.list$theta1.hat, "sigmahat" = mle.list$theta2.hat)
    }
    else if (distn == "llogis") {
      #print(paste0("MLE value is: (lambda.hat = ", theta1.hat, ", kappa.hat = ", theta2.hat,")"))
      crlist <- list("lambda" = crlist$x, "kappa" = crlist$y, "phi" = crlist$phi,
                     "lambdahat" = mle.list$theta1.hat, "kappahat" = mle.list$theta2.hat)
    }
    else if (distn == "gamma") {
      #print(paste0("MLE value is: (theta.hat = ", theta1.hat, ", kappa.hat = ", theta2.hat,")"))
      crlist <- list("theta" = crlist$x, "kappa" = crlist$y, "phi" = crlist$phi,
                     "thetahat" = mle.list$theta1.hat, "kappahat" = mle.list$theta2.hat)
    }
    else if (distn == "unif") {
      #print(paste0("MLE value is: (a.hat = ", theta1.hat, ", b.hat = ", theta2.hat,")"))
      crlist <- list("a" = crlist$x, "b" = crlist$y, "phi" = crlist$phi,
                     "ahat" = mle.list$theta1.hat, "bhat" = mle.list$theta2.hat)
    }
    else if (distn == "cauchy") {
      #print(paste0("MLE value is: (a.hat = ", theta1.hat, ", s.hat = ", theta2.hat,")"))
      crlist <- list("a" = crlist$x, "s" = crlist$y, "phi" = crlist$phi,
                     "ahat" = mle.list$theta1.hat, "shat" = mle.list$theta2.hat)
    }
    return(crlist)
  }
  else if (info == FALSE) {
    return(paste0("Confidence region plot complete; made using ", length(phi)," boundary points."))
  }
}
