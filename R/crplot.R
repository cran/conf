#' Plotting Two-Dimensional Confidence Regions
#'
#' @description
#' Plotting a two-dimensional confidence region for probability distribution parameters (supported distribution
#' suffixes: cauchy, gamma, invgauss, lnorm, llogis, logis, norm, unif, weibull) corresponding to a user given
#' complete or right-censored dataset and level of significance.  See the CRAN website
#' https://CRAN.R-project.org/package=conf for a link to two \code{crplot} vignettes.
#'
#' @param dataset a vector of n data values.
#' @param distn distribution to fit the dataset to; accepted values: \code{'cauchy'}, \code{'gamma'}, \code{'invgauss'},
#' \code{'logis'}, \code{'llogis'}, \code{'lnorm'}, \code{'norm'}, \code{'unif'}, \code{'weibull'}.
#' @param alpha significance level; resulting plot illustrates a 100(1 - \code{alpha})\% confidence region.
#' @param cen a vector of binary values specifying if the corresponding data values are right-censored (0), or
#' observed (1, default); its length must match length(dataset).
#' @param heuristic numeric value selecting method for plotting: 0 for elliptic-oriented point distribution, and
#' 1 for smoothing boundary search heuristic.
#' @param maxdeg maximum angle tolerance between consecutive plot segments in degrees.
#' @param ellipse_n number of roughly equidistant confidence region points to plot using the
#' elliptic-oriented point distribution (must be a multiple of four because its algorithm
#' exploits symmetry in the quadrants of an ellipse).
#' @param pts displays confidence region boundary points identified if \code{TRUE}.
#' @param mlelab logical argument to include the maximum
#' likelihood estimate coordinate point (default is \code{TRUE}).
#' @param sf significant figures in axes labels specified using sf = c(x, y), where x and y represent the optional
#' digits argument in the R function \code{\link{round}} as it pertains to the horizontal and vertical labels.
#' @param mar specifies margin values for \code{par(mar = c( ))} (see \code{mar} in \code{\link{par}}).
#' @param xyswap logical argument to switch the axes that the distribution parameter are shown.
#' @param xlab string specifying the x axis label.
#' @param ylab string specifying the y axis label.
#' @param main string specifying the plot title.
#' @param xlas numeric value of 0, 1, 2, or 3 specifying the style of axis labels (see \code{las} in \code{\link{par}}).
#' @param ylas numeric value of 0, 1, 2, or 3 specifying the style of axis labels (see \code{las} in \code{\link{par}}).
#' @param origin logical argument to include the plot origin (default is \code{FALSE}).
#' @param xlim two-element vector containing horizontal axis minimum and maximum values.
#' @param ylim two-element vector containing vertical axis minimum and maximum values.
#' @param tol the \code{\link{uniroot}} parameter specifying its required accuracy.
#' @param info logical argument to return plot information: MLE is returned as a list; (x, y) plot point coordinates
#' and corresponding phi angles (with respect to MLE) are returned as a list.
#' @param maxcount integer value specifying the number of smoothing search iterations before terminating with \code{maxdeg} not met.
#' @param repair logical argument to repair regions inaccessible using a radial angle from its MLE due to multiple
#' roots at select \eqn{\phi} angles.
#' @param jumpshift see vignette "conf Advanced Options" for details; location (as a fractional value between 0 and 1) along the
#' vertical or horizontal "gap" (near an uncharted region) to locate a jump-center toward; can be either a scalar value (uniformly
#' applied to all jump-centers) or vector of length four (with unique values for its respective quadrants, relative to the MLE).
#' @param jumpuphill see vignette "conf Advanced Options" for details; significance level increase to \code{alpha} for the jump-center
#' (corresponds to an "uphill" location on its likelihood function); can be either a scalar value (uniformly applied to all jump-centers)
#' or vector of length four (with unique values for its respective quadrants, relative to the MLE).
#' @param jumpinfo logical argument to return plot info (see \code{info} argument) and jump-center info; returned within `repair`
#' attribute are \code{jumpuphill} value, \code{jumpshift} value, "|" or "-" gap type, jump-center(s) coordinates, and coordinates
#' of points left & right of the inaccessible region.
#' @param showjump logical argument specifying if jump-center repair reference points appear on the confidence region plot.
#' @param showplot logical argument specifying if a plot is output; altering from its default of \code{TRUE} is
#' only logical assuming \code{crplot} is run for its data only (see the \code{info} argument).
#' @param animate logical argument specifying if an animated plot build will display; the annimation sequence is given in successive plots.
#' @param delay numeric value of delay (in seconds) between successive plots when \code{animate = TRUE}.
#' @param exact logical argument specifying if alpha value is adjusted to compensate for negative coverage bias to achieve
#' (1 - alpha) coverage probability using previously recorded Monte Carlo simulation results; available for limited values of
#' alpha (roughly <= 0.2--0.3), n (typically n = 4, 5, ..., 50) and distributions (distn suffixes: weibull, llogis, norm).
#' @param silent logical argument specifying if console output should be suppressed.
#' @import stats
#' @import graphics
#' @importFrom fitdistrplus mledist
#' @export
#' @return If the optional argument \code{info = TRUE} is included then a list is returned with:
#' \itemize{
#' \item parm1*: a vector containing the associated confidence region boundary values for parameter 1
#' \item parm2*: a vector containing the associated confidence region boundary values for parameter 2
#' \item phi: a vector containing the angles used
#' \item parm1hat*: the MLE for parameter 1
#' \item parm2hat*: the MLE for parameter 2
#' }
#' *Note: "param1" and "param2" are placeholders that will be replaced with the appropriate parameter names
#' based on the probability distribution.
#' @concept confidence region plot graphics visualization coverage parameter estimation
#' @keywords distribution models univar
#' @references A. Jaeger (2016), "Computation of Two- and Three-Dimensional Confidence Regions with the Likelihood Ratio",
#' The American Statistician, 49, 48--53.
#' @references C. Weld, A. Loh, L. Leemis (2020), "Plotting Two-Dimensional Confidence Regions",
#' The American Statistician, Volume 72, Number 2, 156--168.
#' @seealso \code{\link{coversim}}, \code{\link{uniroot}}
#' @author Christopher Weld (\email{ceweld241@gmail.com})
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
#'                 sf        = NULL,
#'                 mar       = c(4, 4.5, 2, 1.5),
#'                 xyswap    = FALSE,
#'                 xlab      = "",
#'                 ylab      = "",
#'                 main      = "",
#'                 xlas      = 0,
#'                 ylas      = 0,
#'                 origin    = FALSE,
#'                 xlim      = NULL,
#'                 ylim      = NULL,
#'                 tol       = .Machine$double.eps ^ 1,
#'                 info      = FALSE,
#'                 maxcount  = 30,
#'                 repair    = TRUE,
#'                 jumpshift = 0.5,
#'                 jumpuphill = min(alpha, 0.01),
#'                 jumpinfo  = FALSE,
#'                 showjump  = FALSE,
#'                 showplot  = TRUE,
#'                 animate   = FALSE,
#'                 delay     = 0.5,
#'                 exact     = FALSE,
#'                 silent    = FALSE )
#'
#' @details
#' This function plots a confidence region for a variety of two-parameter distributions.  It requires:
#' \itemize{
#' \item a vector of dataset values,
#' \item the level of significance (alpha), and
#' \item a population distribution to fit the data to.
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
#' asymptotic results associated with the ratio test statistic \eqn{-2 [\log L(\theta) - \log L(\hat{\theta})]}
#' which converges in distribution to the chi-square distribution with two degrees of freedom (for
#' a two parameter distribution).
#'
#' The default axes convention in use by \code{crplot} are
#'
#' \tabular{lcc}{
#' \tab Horizontal \tab Vertical\cr
#' Distribution  \tab  Axis  \tab Axis\cr
#' Cauchy \tab \eqn{a} \tab \eqn{s}\cr
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
#' \deqn{(1 / \sigma) exp((x - \mu) / \sigma) (1 + exp((x - \mu) / \sigma)) ^ {-2}}
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
                   sf = NULL,
                   mar = c(4, 4.5, 2, 1.5),
                   xyswap = FALSE,
                   xlab = "",
                   ylab = "",
                   main = "",
                   xlas = 0,
                   ylas = 0,
                   origin = FALSE,
                   xlim = NULL,
                   ylim = NULL,
                   tol = .Machine$double.eps ^ 1,
                   info = FALSE,
                   maxcount = 30,
                   repair = TRUE,
                   jumpshift = 0.5,
                   jumpuphill = min(alpha, 0.01),
                   jumpinfo = FALSE,
                   showjump = FALSE,
                   showplot = TRUE,
                   animate = FALSE,
                   delay = 0.5,
                   exact = FALSE,
                   silent = FALSE) {

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

  if (!is.numeric(cen) || !all(cen %in% 0:1) || (length(dataset) != length(cen)))
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

  if (!is.null(sf) && (length(sf) != 2 || !is.numeric(sf) || floor(sf)[1] != sf[1] || floor(sf)[2] != sf[2]))
    stop("'sf' must be a vector of integers with length two")

  if (length(mar) != 4 || !is.numeric(mar) || min(mar) < 0)
    stop("'mar' must be a vector of length four with positive numeric entries")

  if (!is.logical(xyswap) || length(xyswap) != 1)
    stop("'xyswap' must be a single logical parameter")

  if (!(xlas %in% c(0, 1, 2, 3)))
    stop("'xlas' must be a numeric value in {0, 1, 2, 3}")

  if (!(ylas %in% c(0, 1, 2, 3)))
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
    stop("'tol' numeric parameter given is invalid (default .Machine$double.eps^1)")

  if (!is.logical(info) || length(info) != 1)
    stop("'info' must be a single logical parameter")

  if (!is.logical(jumpinfo) || length(jumpinfo) != 1)
    stop("'info' must be a single logical parameter")

  if (jumpinfo && !repair)
    warning("'jumpinfo' is not applicable when 'repair' is FALSE")

  if (!is.logical(repair) || length(repair) != 1)
    stop("'repair' must be a single logical parameter")

  if (min(jumpshift) <= 0 || max(jumpshift) >= 1 || !is.numeric(jumpshift) || ((length(jumpshift) != 1) && (length(jumpshift) != 4)))
    stop("'jumpshift' must be numeric value(s) between 0 and 1 with length 1 (same everywhere) or 4 (unique to each quadrant, relative to MLE)")
  #stop("'jumpshift' must be a single numeric value 0 < jumpshift < 1, ")

  if (jumpuphill <= 0 || !is.numeric(jumpuphill) || ((length(jumpuphill) != 1) && (length(jumpuphill) != 4)))
    stop("'jumpuphill' must numeric value(s) such that 0 < (alpha + jumpuphill) < 1 with length 1 (same everywhere) or 4 (unique to each quadrant, relative to MLE)")
  #stop("'jumpuphill' must be a single numeric value such that 0 < (alpha + jumpuphill) < 1")

  if (!is.logical(showjump) || length(showjump) != 1)
    stop("'showjump' must be a single logical parameter")

  if (showjump && !showplot)
    warning("'showjump' is TRUE but will not appear since 'showplot' is FALSE")

  if (!is.logical(showplot) || length(showplot) != 1)
    stop("'showplot' must be a single logical parameter")

  if (!showplot && !info && !jumpinfo && !animate)
    warning("'showplot', 'info', 'jumpinfo', and 'animate' are all FALSE; without these, crplot returns nothing to its user")

  if (!is.logical(animate) || length(animate) != 1)
    stop("'animate' must be a single logical parameter")

  if (animate && (heuristic == 0))
    warning("'animate' is not applicable when 'heuristic = TRUE' (not an iterative build process)")

  if (!is.numeric(delay) || length(delay) != 1 || delay < 0 )
    stop("'delay' must be a non-negative numeric scalar value")

  if (!animate && (delay != 0.5))
    warning("'delay' is not applicable unless 'animate = TRUE'")

  if (!is.logical(exact) || length(exact) != 1)
    stop("'exact' must be a single logical parameter")

  if (exact && !(distn %in% c("weibull", "llogis", "norm")))
    warning("'exact' not available for this distn; proceeding without it")

  if (!is.logical(silent) || length(silent) != 1)
    stop("'silent' must be a single logical parameter")

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
  # 8. nomadjust.      Adjusts nominal value (1 - alpha) to compensate for negative bias inherent in small
  #                    sample sizes.  Uses tables assembled via 20 million Monte Carlo simulation replications
  #                    of the conf coversim() function.  Available for a limited range of alpha values
  #                    (roughly <= 0.25) and distributions (distn suffixes: weibull, llogis, norm).
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
  # Gamma and log logistic use the STAR package code from Christophe Pouzat to identify their MLEs.
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
      # temp1 <- STAR::invgaussMLE(x, si = as.double(cen))     # STAR requires cen package as double
      temp1 <- invgaussMLE(x, si = as.double(cen))     # invgaussMLE code previously from STAR package
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
      # temp1 <- STAR::llogisMLE(x, si = as.double(cen))
      temp1 <- llogisMLE(x, si = as.double(cen))   # llogisMLE code previously from STAR package
      lambda.hat <- 1 / exp(temp1$estimate[['location']])
      kappa.hat <- 1 / temp1$estimate[['scale']]
      mleLL <- temp1$logLik                       # log likelihood value at its maximum
      mle.list <- list("theta1.hat" = lambda.hat, "theta2.hat" = kappa.hat, "mleLLvalue" = mleLL)
    }

    # gamma MLE
    else if (distn == "gamma"){
      # temp1 <- STAR::gammaMLE(x, si = as.double(cen))
      temp1 <- gammaMLE(x, si = as.double(cen))  # gammaMLE code previously from STAR package
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
  # Some ASEs are approximations based on the ASE of other distributions (after an appropriate
  # conversion to the parameterization in use by crplot).
  # for future reference only (this approach not yet taken):
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
      se <- suppressMessages(sqrt(OI))
      theta1ase <- se[1]    # a location parameter
      theta2ase <- se[4]    # alpha scale parameter
    }
    else if (distn == "gamma") {
      # theta1ase <- as.numeric(STAR::gammaMLE(x, si = as.double(cen))$se[2])    # scale parameter theta
      # theta2ase <- as.numeric(STAR::gammaMLE(x, si = as.double(cen))$se[1])    # shape parameter kappa
      # gammaMLE code previously from STAR package
      theta1ase <- as.numeric(gammaMLE(x, si = as.double(cen))$se[2])    # scale parameter theta
      theta2ase <- as.numeric(gammaMLE(x, si = as.double(cen))$se[1])    # shape parameter kappa
    }
    else if (distn == "invgauss") {
      # theta1ase <- as.numeric(STAR::invgaussMLE(x, si = as.double(cen))$se[1])       # mu
      # inv_theta2ase <- as.numeric(STAR::invgaussMLE(x, si = as.double(cen))$se[2])   # (1 / theta) se
      # invgaussMLE code previously from STAR package
      theta1ase <- as.numeric(invgaussMLE(x, si = as.double(cen))$se[1])       # mu
      inv_theta2ase <- as.numeric(invgaussMLE(x, si = as.double(cen))$se[2])   # (1 / theta) se
      theta2ase <- abs(1 / ((1 / theta2.hat) + (inv_theta2ase / 2)) -     # theta se approximation
                         1 / ((1 / theta2.hat) - (inv_theta2ase / 2)))
    }
    else if (distn == "norm") {
      left <- replace(x, which(cen == -1), rep(NA, length(which(cen == -1))))
      right <- replace(x, which(cen == 0), rep(NA, length(which(cen == 0))))
      xx <- data.frame(left = left, right = right)
      hes <- fitdistrplus::mledist(xx, "norm", silent = TRUE)$hessian       # using plotdistrplus
      OI <- solve(hes)
      se <- suppressMessages(sqrt(OI))
      theta1ase <- se[1]    # mu
      theta2ase <- se[4]    # sigma
    }
    else if (distn == "logis") {
      left <- replace(x, which(cen == -1), rep(NA, length(which(cen == -1))))
      right <- replace(x, which(cen == 0), rep(NA, length(which(cen == 0))))
      xx <- data.frame(left = left, right = right)
      hes <- fitdistrplus::mledist(xx, "logis", silent = TRUE)$hessian       # using plotdistrplus
      OI <- solve(hes)
      se <- suppressMessages(sqrt(OI))
      theta1ase <- se[1]    # mu (location)
      theta2ase <- se[4]    # sigma (scale)
    }
    else if (distn == "llogis") {
      # alt_param <- STAR::llogisMLE(x, si = as.double(cen))
      alt_param <- llogisMLE(x, si = as.double(cen))        # llogisMLE code previously from STAR package
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
      se <- suppressMessages(sqrt(OI))
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
      se_kappa <- suppressMessages(try(sqrt(OI[1]), silent = TRUE))
      se_lambda <- suppressMessages(try(sqrt(OI[4]), silent = TRUE))
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
        # llfn <- sum(cen * log(STAR::dinvgauss(x, mu = mu, sigma2 = 1 / lambda))) +
        #   sum(as.numeric(cen==0) * log((1 - STAR::pinvgauss(x, mu = mu, sigma2 = 1 / lambda))))
        # dinvgauss and pinvgauss code originally from the STAR package
        llfn <- sum(cen * log(dinvgauss(x, mu = mu, sigma2 = 1 / lambda))) +
          sum(as.numeric(cen==0) * log((1 - pinvgauss(x, mu = mu, sigma2 = 1 / lambda))))
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
      else {            # censored values:
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
      #for (umult in c(10, 100, 500)) {   # uniroot 'upper' argument will use this multiplier to set upper bounds in search of root (>> in 1st quad; using theta.hat in others)
      for (umult in c(3, 10, 100, 500)) {   # uniroot 'upper' argument will use this multiplier to set upper bounds in search of root (>> in 1st quad; using theta.hat in others)
        #print(paste0("umult value: ", umult, " --------------"))
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
            tempUpper <- temp * 0.99               # set search limit just before upper limit (y-axis); 0.99 seems most robust here (trial & error assessment)
          }
          else {
            temp <- theta2.hat/ (-sin(phi[j]))   # an upper bound for confidence region feasible points; -sin() because sin(phi[j]) < 0
            tempUpper <- temp * 0.99999           # set search limit just before x-axis; 0.9999 seems most robust here (trial & error assessment)
          }
          if ((tempUpper > umult * theta1.hat) && (tempUpper > umult * theta2.hat)) {
            tempUpper <- umult * max(c(theta1.hat, theta2.hat))     # arbitrary upper bound for phi near pi/2; accept risk CR bound <= umult * max(mle_parameter)
          }
          #print(paste0("phi is: ", phi[j]))
          #print(paste0("tempUpper is: ", tempUpper))
          g <- suppressWarnings(try(g1 <- uniroot(llrsolve, lower = 0, upper = abs(tempUpper),
                                                  phi = phi[j], MLE = MLEHAT, mleLLvalue = mleLLvalue, x = samp, #extendInt = "downX",
                                                  cen = cen, chi2 = qchisq(1 - alpha, 2), tol = tol),
                                    silent = TRUE))
          # if (class(g) == "list") {
          #   if (g$root == 0) {
          #     #print(g)
          #     print("THIS IS A PROBLEM.....   g is zero!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  at...")
          #     #print(phi[j])
          #     lowv <- 0
          #     while (g$root == lowv) {
          #       #if (lowv == 0) {
          #       #  lowv <- theta2.hat / 10
          #       #}
          #       #else {
          #         lowv <- lowv + (theta2.hat - lowv) * 0.1
          #       #}
          #       print("low and upper:")
          #       print(lowv)
          #       print(c(tempUpper, theta2.hat))
          #       if ((phi[j] == 3*pi / 2) ||  (phi[j] == pi / 2)) {
          #         print("vertical")
          #         g <- suppressWarnings(try(g1 <- uniroot(llrsolve, lower = lowv, upper = abs(tempUpper),
          #                                                 phi = phi[j], MLE = MLEHAT, mleLLvalue = mleLLvalue, x = samp,
          #                                                 cen = cen, chi2 = qchisq(1 - alpha, 2), tol = tol), silent = TRUE))
          #       }
          #     }
          #     if ((phi[j] == 0) || (phi[j] == pi) || (phi[j] == 2 * pi)) {
          #       print("horizontal")
          #     }
          #   }
          #}

          # errors in uniroot (calculating "g") are possible, and a result of poor upper bound selection
          # the next section of code troubleshoots those errors in two ways:
          # 1. it decreases the upper value by a factor of multipliers (see shiftHi below), and also
          # 2. it identifies a relatively small upper value (based on min(theta1.hat, theta2.hat)) and then increments it
          # The first option is valuable when the upper limit couldn't be evaluated (-infinity)
          # The second option is valuable when vastly different parameter orders of magnitude push the search too far for phi calcs ~parallel the smaller axis
          z <- 0                       # counter for "try-error" recovery attempts below
          tempUpperLo <- 0             # initialize; will incrementally push upward seeking uniroot upper param
          tempUpperHi <- tempUpper     # initialize; will incrementally push upward seeking uniroot upper param
          #if(class(g) == "list") {print(g$root)}

          # the more shiftHi multipliers, the more chances uniroot has to "fix" any numeric errors, however, this comes at a computational cost
          shiftHi <- seq(.9, .5, by = -.1)^(c(1:5))  # "upper" multipliers of c(0.900 0.640 0.343 0.1296 0.03125) when try-error
          # other schemes / options tried include:
          #shiftHi <- seq(0.95, 0.05, by = -0.1)
          #shiftHi <- c(0.99, seq(.9, .5, by = -.1)^(c(1:5)))  # "upper" multipliers of c(0.99 0.900 0.640 0.343 0.1296 0.03125) when try-error
          #shiftHi <- 0.75 ^ (c(1:6))
          #shiftHi <- seq(0.99, 0.01, by = -0.01)
          # explaination of baseUpperLo calculation given in three print statements following its calculation:
          baseUpperLo <- (tempUpper * min(shiftHi) - min(c(abs(theta1.hat), abs(theta2.hat)))) ^ (1 / (length(shiftHi) + 1))
          #print(paste0("lowest tempUpperHi value for umult ", umult, " will be: ", (tempUpper * min(shiftHi))))
          #print(paste0("baseUpperLo: ", baseUpperLo, " and length(shift(Hi)) = ", length(shiftHi)))
          #print(paste0("...so highest tempUpperLo will be: ", baseUpperLo ^ length(shiftHi)))
          while ((class(g) == "try-error") && (z < length(shiftHi))) {
            #while ((class(g) == "try-error") && (z < length(shiftHi)) && (tempUpperLo < tempUpperHi)) {     # tempUpperLo and tempUpperHi don't cross b/c baseUpperLo calc
            z <- z + 1
            #print(c(z, shiftHi[z]))
            tempUpperHi <- tempUpper * shiftHi[z]        # was 0.5 ^ z in v1.4.0 instead of shiftHi[z]; 0.75 ^ z also good option
            #tempUpperLo <- 5 ^ (z - 1) * min(c(abs(theta1.hat), abs(theta2.hat)))
            tempUpperLo <- baseUpperLo ^ (z - 1) * min(c(abs(theta1.hat), abs(theta2.hat)))
            #print(paste0("-------------- problematic phi value: ", phi[j]))
            #print(c(tempUpperLo, tempUpperHi))
            #points(theta1.hat, theta2.hat, cex = 3)        # highlight MLE / jump-center
            #points(theta1.hat - tempUpper * sin(3*pi/2 - phi[j]), theta2.hat - tempUpper * cos(3*pi/2 - phi[j]), cex = 2)                      # highlight problematic angle
            #points(theta1.hat - tempUpperHi * sin(3*pi/2 - phi[j]), theta2.hat - tempUpperHi * cos(3*pi/2 - phi[j]), cex = 2, col = "blue")    # highlight problematic angle
            #points(theta1.hat - tempUpperLo * sin(3*pi/2 - phi[j]), theta2.hat - tempUpperLo * cos(3*pi/2 - phi[j]), cex = 2, col = "orange")  # highlight problematic angle
            #lines(c(theta1.hat, theta1.hat - tempUpper * sin(3*pi/2 - phi[j])), c(theta2.hat, theta2.hat - tempUpper * cos(3*pi/2 - phi[j])), col = "red")
            #print(paste0(z, " ****************************"))
            #print(paste0("Error; correction sought with uniroot upper bound modifications... ", tempUpperLo, " and ", tempUpperHi)))
            g <- suppressWarnings(try(g1 <- uniroot(llrsolve, lower = 0, upper = abs(tempUpperHi),
                                                    phi = phi[j], MLE = MLEHAT, mleLLvalue = mleLLvalue, x = samp, #extendInt = "downX",
                                                    cen = cen, chi2 = qchisq(1 - alpha, 2), tol = tol),
                                      silent = TRUE))
            # if (class(g) == "try-error") {   # replace due to CRAN checks
            if (inherits(g, "try-error")) {
              #print(paste0("...............try pushing up the min value to ", tempUpperLo))
              g <- suppressWarnings(try(g1 <- uniroot(llrsolve, lower = 0, upper = abs(tempUpperLo),
                                                      phi = phi[j], MLE = MLEHAT, mleLLvalue = mleLLvalue, x = samp, #extendInt = "downX",
                                                      cen = cen, chi2 = qchisq(1 - alpha, 2), tol = tol),
                                        silent = TRUE))
              # if (class(g) != "try-error") {    # replace due to CRAN checks
              if (!inherits(g, "try-error")) {
                #print(paste0("tempUpperLo fix value of: ", tempUpperLo, " with umult: ", umult))
              }
            }
            else {
              #print(paste0("...*error overcome* through uniroot 'upper' argument modification by tempUpperHi factor ", shiftHi[z],
              #             " with umult: ", umult, "; algorithm resuming..."))
            }
          }
          # if (class(g) != "try-error") {   # replace due to CRAN checks
          if (!inherits(g, "try-error")) {
            done <- 1
          }
        }   # end if (done == 0)
      }     # end for (umult in ...)

      # if (class(g) == "try-error") {   # replace due to CRAN checks
      if (inherits(g, "try-error")) {
        print("-------------------------------------------------------------------------------------------------")
        print("R uniroot failure searching for confidence region boundary---challenging parameters and/or shape.")
        print("Unable to produce a confidence region for the given sample and/or parameterization.              ")
        #print("Uniroot failure.  Challenging parameters and/or shape requires customizing uniroot bounds.")
        print("-------------------------------------------------------------------------------------------------")
        stop()
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
    # repairpass = FALSE if 1st call to crsmooth (radial angles taken from MLE), and TRUE o.w. (jump-center repairs)
    # repairinfo stores jump-center & other relevant info when repairpass = TRUE
    # repairq identifies the current "quadrant" (relevant to MLE) being repaired when repairpass = TRUE
    # repaircrlist stores crlist x, y, and phi values when repairpass = TRUE

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

    # tilt <- - pi / 4                                  # "tilt" cardinal coord directions to avoid numeric difficulties at ~0, pi/2,...
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
      ase <- try(asesolve(x = dataset, cen, theta1.hat, theta2.hat), silent = TRUE)
      #}

      # circumstance where ASE is unatainable include computational singular observed info matrix (ase returned as a
      # "character" string with error message) or NaN result (both filtered with "if" statement below)
      if (is.list(ase) && !is.nan(ase$theta1.ase) && !is.nan(ase$theta2.ase)) {
        phi <- angles(a = ase$theta1.ase, b = ase$theta2.ase, npoints = ellipse_n)
      }
      else  {        # if ASE unavailable, estimate aspect ratio from MLE
        #message("ASE calc unsuccessful; using MLE to estimate aspect ratio")
        phi <- angles(a = theta1.hat, b = theta2.hat, npoints = ellipse_n)
      }
      #print("HERE IS PHI")
      #print(phi)

    }
    #}

    ## (repair alternatives below did not save time and/or were subject to lack of graph detail & are therefore commented-out)
    ## (however, kept for future reference)
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
    cr <- crsolve(samp, cen = cen, alpha = alpha, mle.list = mle.list, phi)      # confidence region points w initial phi's

    ##### check if any uniroot solutions returned ~0 (i.e. identical location to MLE), and discard if necessary
    # note: a numeric difficulty sometimes arises when phi angle is in {0, pi/2, pi, 3pi/2, 2pi} causing this error
    invalid_pt <- intersect(which(cr[, 1] == theta1.hat), which(cr[, 2] == theta2.hat))
    #print(invalid_pt)
    if (length(invalid_pt) != 0) {
      count <- 0
      adj_ellipse_n <- ellipse_n
      # if invalid CR point (== MLE),
      while (((adj_ellipse_n - length(invalid_pt)) < 6) && (count < 6)) {
        count <- count + 1
        adj_ellipse_n <- 4 * (2 ^ count)
        #print(paste0("trying: ", adj_ellipse_n))
        ase <- try(asesolve(x = dataset, cen, theta1.hat, theta2.hat), silent = TRUE)
        # circumstance where ASE is unatainable include computational singular observed info matrix (ase returned as a
        # "character" string with error message) or NaN result (both filtered with "if" statement below)
        if (is.list(ase) && !is.nan(ase$theta1.ase) && !is.nan(ase$theta2.ase)) {
          phi <- angles(a = ase$theta1.ase, b = ase$theta2.ase, npoints = adj_ellipse_n)
        }
        else  {        # if ASE unavailable, estimate aspect ratio from MLE
          #message("ASE calc unsuccessful; using MLE to estimate aspect ratio")
          phi <- angles(a = theta1.hat, b = theta2.hat, npoints = adj_ellipse_n)
        }
        cr <- crsolve(samp, cen = cen, alpha = alpha, mle.list = mle.list, phi = phi)

        invalid_pt <- intersect(which(cr[, 1] == theta1.hat), which(cr[, 2] == theta2.hat))
        #print(paste0("now invalid: ", length(invalid_pt)))
      }
      cr <- cr[-c(invalid_pt), ]
      phi <- phi[-c(invalid_pt)]
      #print("REMOVED INDEX #(s):")
      #print(invalid_pt)
      #print(phi)
    }

    transf.scale <- 1                     # initialize transformation scaling factor at 1
    count <- 0                            # administrative counter for while loop
    d01 <- 0
    d02 <- 0
    d12 <- 0
    while (maxcurrent > maxrad) {                                          # in radians
      #print(c(maxrad, maxcurrent))
      count <- count + 1
      if (count != 1){
        index_phinew <- match(phinew, phi)
        index_phiold <- match(phiold, phi)
        #print("check1")
        crnew <- crsolve(samp, cen = cen, alpha = alpha, mle.list = mle.list, phi = phinew)     # uniroot calcs for CR points of new phis
        #print("check2")
        #points(crnew, col = "yellow", pch = 16, cex = 0.5)
        #Sys.sleep(0.02)
        crold <- cr
        cr <- matrix(rep_len(0, length(phi) * 2), length(phi), 2)          # initialize to store CR coordinate
        cr[index_phiold, ] <- crold
        cr[index_phinew, ] <- crnew

      }
      phiold <- phi                        # store phi values for next iteration
      if (!repairpass) {
        xspan <- max(cr[,1]) - min(cr[,1])
        yspan <- max(cr[,2]) - min(cr[,2])
      }
      else if (repairpass) {
        xspan <- max(c(cr[,1], repaircrlist$x)) - min(c(cr[,1], repaircrlist$x))
        yspan <- max(c(cr[,2], repaircrlist$y)) - min(c(cr[,2], repaircrlist$y))
        #if (length(cr[,1]) == 4) {       # print first first pass of repair points (in cardinal directions from jump-center)
        #  print(cr)
        #}
      }
      transf.scale <- xspan / yspan
      cr.transf <- cr
      cr.transf[, 2] <- transf.scale * cr[, 2]

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

      elbowcalc <- (d01^2 + d02^2 - d12^2)/(2 * d01 * d02)    # law of cosines acos(__) parameter
      same <- sort(union(which(d01 == 0), which(d02 == 0)))   # ID any consecutive points that are virtually identical (dist to prev/next = 0)
      elbowcalc[same] <- rep(-1, length(same))                # for same points, make resulting elbow = 0 (no bend for same points); avoids NaN error
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
      if (count != 1) {philast <- phinew}                   # store previous phi additions for plotting
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

      ######### subsection to create .gif of a plot build
      # setwd("/Users/insert_path_here")
      #
      # # # after all png files are in place, to create the .gif run ImageMagick with:
      # # library(magick)
      # my_command <- 'convert *.png -delay 100 -loop 0 animation3.gif'
      # system(my_command)
      # # # then go to: https://ezgif.com/speed to slow down the speed
      #
      # # rename function creates file name with leading zeros
      # # makes it easier to process them sequentially
      # rename <- function(x, y){
      #   if (x < 10) {
      #     return(name <- paste('000', x,'plot.png',sep=''))
      #   }
      #   if (x < 100 && x >= 10) {
      #     return(name <- paste('00', x,'plot.png', sep=''))
      #   }
      #   if (x >= 100) {
      #     return(name <- paste('0', x,'plot.png', sep=''))
      #   }
      # }
      #
      # # Name & Save .png file for this snapshot of the build
      # # note: to capture jump-center repairs, it is necessary to record multiple times and increase count (i.e. to count + 30) to capture successive regions using different file names
      # name <- rename(count)   # name the next plot in the plot-build sequence
      # # if (repairpass && (repairq == 4)) {stop()}  # with multiple jump-center repairs, may need to use something like this to capture specific quadrants
      # png(name)               # saves the plot as a .png file in the working directory
      #
      ######### (end of png file save setup; NOTE: also need to un-commend next dev.off() below)

      # this sub-section plots the progression of added points, showing interim steps:
      if (animate) {
        par(xpd = FALSE)
        par(mar = c(3, 3, 1.5, 1.5))
        if (repairpass) {
          plot(c(repaircrlist$x, cr[,1]), c(repaircrlist$y, cr[,2]),
               #main = "in-progress build of confidence region\n(jump-center repairs)",
               axes = FALSE, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
        }
        else {
          plot(cr, #main = "in-progress build of confidence region",
               axes = FALSE, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
        }
        lines(cr, lty = 4, col = "black", lwd = 0.9)
        big <- max(theta1.hat, theta2.hat) * 100
        if (count != 1) {
          segments(rep(theta1.hat, length(philast)), rep(theta2.hat, length(philast)),
                   theta1.hat + big * cos(philast), theta2.hat + big * sin(philast), col = "yellow3")
        }
        else if (ellipse_n != 4) {
          segments(rep(theta1.hat, length(phi)), rep(theta2.hat, length(phi)),
                   theta1.hat + big * cos(phi), theta2.hat + big * sin(phi), col = "yellow3")
        }
        else {
          segments(rep(theta1.hat, 4), rep(theta2.hat, 4),
                   c(theta1.hat + big, theta1.hat, theta1.hat - big, theta1.hat),
                   c(theta2.hat, theta2.hat + big, theta2.hat, theta2.hat - big), col = "yellow3")
        }
        points(cr, pch = 21, bg = c("red", "blue", "red", "blue", "red", "blue", "red", "green")[goodtot + 1],
               col = c("firebrick4", "blue4", "firebrick4", "blue4", "firebrick4", "blue4", "firebrick4", "chartreuse4")[goodtot + 1])
        segments(cr[dim(cr)[1], 1], cr[dim(cr)[1], 2], cr[1, 1], cr[1, 2], col = 'black', lty = 4) # connect first and last points
        points(theta1.hat, theta2.hat, pch = 8, cex = 0.8, col = 'gray30')
        legend("topright", legend = c("complete", "near-complete", "incomplete"), lty = c(4, 4, 4), pch = c(1, 21, 21, 21),
               col = c("chartreuse4", "blue4", "firebrick4"), pt.bg = c("green", "blue", "red"), cex = 0.7,
               title = paste0(length(cr[,1]), " boundary points, iteration: ", count), bty = "n")
        #        legend("topright", legend = c("complete", "near-complete", "incomplete"), lty = 4, pch = c(21, 21, 21),
        #               col = c("chartreuse4", "blue4", "firebrick4"), pt.bg = c("green", "blue", "red"), cex = 0.8)
        axis(side = 1)
        axis(side = 2)
        Sys.sleep(delay)      # pause to view plot
      }                     # end if (animate)

      # dev.off()           # uncomment when saving .png for .gif annimation build

      crlist <- list("x" = cr[, 1], "y" = cr[, 2], "phi" = phi)
      phinew <- sort(unique(phinew))
      # (below commented-out section is problematic for implimentation but kept here for possible future attempts)
      # Its intent was to cut off search when phinew points outside the region of interest were all that remain
      # ...doing so, however, caused pre-mature termination in some circumstance (i.e. crplot(c(1.9, 2), "invgauss", 0.05))
      #if (!is.null(repairinfo) && (count > 10)) {
      #  phinew <- phinew[intersect(which(phinew >= min(repairborderphis[1], repairborderphis[2])),
      #                             which(phinew <= max(repairborderphis[1], repairborderphis[2])))]
      #  if (length(phinew) == 0) {
      #    count <- maxcount  # force exit of while loop because phinew additions no longer apply to region between borderphis
      #  }
      #}
      if (sum(phinew) != 0) {phi <- sort(unique(c(phi, phinew))) }
      if (count >= maxcount) {                              # indicates there is an inaccessible region; 20--50 typically is sufficient

        ############################################################################################
        # May 2018 addition enables charting previously inaccessible regions by identifying
        # jump-center locations away from the MLE but within the CR where a similar smoothing
        # algorithm is done using the radial log likelihood approach.  Recursively calls crsmooth.
        ############################################################################################
        if ((!repairpass) && (repair)) {        # enter if repairs requested; repairs are done within this if statement

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

            # ensure a user-entered jumpuphill value is feasible
            for (z in 1:4) {
              if ((alpha + jumpuphill[z]) > 1) {
                warning(paste0("quad", z, " infeasible jumpuphill value ((alpha + jumpuphill) > 1); restoring default of min(alpha, 0.01)"))
                jumpuphill[z] <- min(alpha, 0.01)
                #jumpuphill <- 0.9 * alpha + 0.1     # (alpha + 0.1(1 - alpha))
              }
            }

            #####################################################
            # Quadrant I (with respect to MLE) repairs
            if (any(phinewstore < (pi / 2))) {
              #print("Creating jump-point to access and patch unreachable area (in Quad I relative to MLE)")
              q <- 1      # quad I

              # identify the borders of the unreachable area
              quadindex <- which(phinewstore <= pi / 2)                                             # unique to this quadrant
              leftphi[q] <- max(phinewstore[quadindex])                                            # unique to this quadrant
              leftx[q] <- cr[min(which(crlist$phi > leftphi[q])) + 1, 1]
              lefty[q] <- cr[min(which(crlist$phi > leftphi[q])) + 1, 2]
              rightphi[q] <- min(phinewstore[quadindex])                                           # unique to this quadrant
              rightx[q] <- cr[max(which(crlist$phi < rightphi[q])) - 1, 1]
              righty[q] <- cr[max(which(crlist$phi < rightphi[q])) - 1, 2]
              #points(leftx[q], lefty[q], pch = 16, col = "green")                # show inaccessible region boundary
              #points(rightx[q], righty[q], pch = 16, col = "blue")               # show inaccessible region boundary

              # create a "new" center (instead of MLE) to access previously unreachable CR points as follows:
              # assess shift range available in y-direction within CR boundary,
              # determine MLE phi corresponding to a % along that available shift range
              # offset slightly inside CR via "uphill" alpha adjustment to locate jump x-y coordinates

              # identify if uncharted region is above or below border of unreachable are
              # (admin note: no vertical gap, "|" case, has been seen / verified to this point)
              if (cr[max(which(crlist$phi < rightphi[q])) - 2, 2] < cr[min(which(crlist$phi > leftphi[q])) + 2, 2]) {
                gaptype[q] <- "|"   # vertical gap segment that jumpphi will cross
                shiftmaxindex <- min(intersect(which(cr[, 1] < rightx[q]), which(cr[, 2] > righty[q])))     # unique to this quadrant
                shift <- jumpshift[q] * (cr[shiftmaxindex, 2] - righty[q])                                     # unique to this quadrant
                jumpphi[q] <- atan((righty[q] - theta2.hat + shift) / (rightx[q] - theta1.hat))             # unique to this quadrant
              }
              else if (cr[max(which(crlist$phi < rightphi[q])) - 2, 2] > cr[min(which(crlist$phi > leftphi[q])) + 2, 2]) {
                gaptype[q] <- "-"
                shiftmaxindex <- max(intersect(intersect(which(cr[, 1] > leftx[q]), which(cr[, 2] < lefty[q])),   # unique to this quadrant
                                               which(cr[, 2] > theta2.hat)))                          # cannot be a quad 4 point
                shift <- jumpshift[q] * (cr[shiftmaxindex, 1] - leftx[q])                                # unique to this quadrant
                jumpphi[q] <- atan((lefty[q] - theta2.hat) / (leftx[q] - theta1.hat + shift))         # unique to this quadrant
              }
              jumpxy[q,] <- crsolve(samp, cen, alpha = alpha + jumpuphill[q], mle.list = mle.list, phi = jumpphi[q])   # ID jump-center location
              #points(cr[shiftmaxindex, 1], cr[shiftmaxindex, 2], pch = 16, col = "brown")   # show jump-center "target" on CR boundary (before jumpuphill)
              #points(theta1.hat + leftx[q] - theta1.hat + shift, theta2.hat + lefty[q] - theta2.hat, col = "red", pch = 16)
              #points(jumpxy[q, 1], jumpxy[q, 2], pch = 16, col = "red")                       # plot jump-center location
              phi1[q] <- pi/2 + atan((jumpxy[q, 1] - rightx[q]) / (righty[q] - jumpxy[q, 2]))   # angle from jump-center to inaccessible region border point 1 (smaller phi angle)
              phi2[q] <- pi + atan((jumpxy[q, 2] - lefty[q]) / (jumpxy[q, 1] - leftx[q]))     # angle from jump-center to inaccessible region border point 2 (larger phi angle)

            }

            #####################################################
            # Quadrant II (with respect to MLE) identification
            if (length(phinewstore[(phinewstore > pi / 2) & (phinewstore <= pi)]) > 0) {
              #print("Creating jump-point to access and patch unreachable area (in Quad II relative to MLE)")
              q <- 2      # quad II

              # identify the boarders of the unreachable area
              # offset two points above and below the inaccessible region to ensure it is partitioned with one reference point on
              # each side (o.w. rare roundoff issues may result with both co-located)
              quadindex <- intersect(which(phinewstore > (pi / 2)), which(phinewstore <= pi))    # unique to this quadrant
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

              # identify if uncharted region is above or below border of unreachable are
              # (admin note: no vertical gap, "|" case, has been seen / verified to this point)
              if (cr[max(which(crlist$phi < rightphi[q])) - 2, 2] > cr[min(which(crlist$phi > leftphi[q])) + 2, 2]) {
                gaptype[q] <- "|"   # vertical gap segment that jumpphi will cross
                shiftmaxindex <- max(intersect(which(cr[, 1] > leftx[q]), which(cr[, 2] > lefty[q])))         # unique to this quadrant
                shift <- jumpshift[q] * (cr[shiftmaxindex, 2] - lefty[q])                                      # unique to this quadrant
                jumpphi[q] <- pi / 2 + atan((theta1.hat - leftx[q]) / ((lefty[q] + shift) - theta2.hat))      # unique to this quadrant
              }
              else if (cr[max(which(crlist$phi < rightphi[q])) - 2, 2] < cr[min(which(crlist$phi > leftphi[q])) + 2, 2]) {
                gaptype[q] <- "-"
                shiftmaxindex <- min(intersect(which(cr[, 1] < rightx[q]), which(cr[, 2] < righty[q])))         # unique to this quadrant
                shift <- jumpshift[q] * (rightx[q] - cr[shiftmaxindex, 1])                                      # unique to this quadrant
                jumpphi[q] <- pi / 2 + atan((theta1.hat - (rightx[q] - shift)) / (righty[q] - theta2.hat))      # unique to this quadrant
              }

              jumpxy[q,] <- crsolve(samp, cen, alpha = alpha + jumpuphill[q], mle.list = mle.list, phi = jumpphi[q])
              #points(jumpxy[q, 1], jumpxy[q, 2], pch = 16, col = "red")
              phi1[q] <- pi / 2 + atan((jumpxy[q,1] - leftx[q]) / (lefty[q] - jumpxy[q,2]))      # region of interest: (0, phi1)
              phi2[q] <- 2 * pi - atan((jumpxy[q,2] - righty[q]) / (rightx[q] - jumpxy[q,1]))    # region of interest: (phi2, 2pi]
            }

            #####################################################
            # Quadrant III (with respect to MLE) identification
            if (length(phinewstore[(phinewstore > pi) & (phinewstore <= 3 * pi / 2)]) > 0) {
              #print("Creating jump-point to access and patch unreachable area (in Quad III relative to MLE)")
              q <- 3      # quad III

              # identify the borders of the unreachable area
              quadindex <- intersect(which(phinewstore > pi), which(phinewstore <= 3 * pi / 2))  # unique to this quadrant
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

              # identify if uncharted region is above or below border of unreachable are
              if (cr[max(which(crlist$phi < rightphi[q])) - 2, 2] > cr[min(which(crlist$phi > leftphi[q])) + 2, 2]) {
                gaptype[q] <- "|"   # vertical gap segment that jumpphi will cross
                shiftmaxindex <- min(intersect(which(cr[, 1] > rightx[q]), which(cr[, 2] < righty[q])))      # unique to this quadrant
                shift <- jumpshift[q] * (righty[q] - cr[shiftmaxindex, 2])                                   # unique to this quadrant
                jumpphi[q] <- pi + atan((theta2.hat - righty[q] + shift) / (theta1.hat - rightx[q]))         # unique to this quadrant
              }
              else if (cr[max(which(crlist$phi < rightphi[q])) - 2, 2] < cr[min(which(crlist$phi > leftphi[q])) + 2, 2]) {
                gaptype[q] <- "-"
                shiftmaxindex <- max(intersect(which(cr[, 1] < leftx[q]), which(cr[, 2] > lefty[q])))      # unique to this quadrant
                shift <- jumpshift[q] * (leftx[q] - cr[shiftmaxindex, 1])                                   # unique to this quadrant
                jumpphi[q] <- pi + atan((theta2.hat - lefty[q]) / (theta1.hat - leftx[q] + shift))         # unique to this quadrant
              }
              jumpxy[q,] <- crsolve(samp, cen, alpha = alpha + jumpuphill[q], mle.list = mle.list, phi = jumpphi[q])
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

              # identify if uncharted region is above or below border of unreachable are
              if (cr[max(which(crlist$phi < rightphi[q])) - 2, 2] < cr[min(which(crlist$phi > leftphi[q])) + 2, 2]) {
                gaptype[q] <- "|"   # vertical gap segment that jumpphi will cross
                shiftmaxindex <- max(intersect(which(cr[, 1] < leftx[q]), which(cr[, 2] < lefty[q])))       # unique to this quadrant
                shift <- jumpshift[q] * (lefty[q] - cr[shiftmaxindex, 2])                                   # unique to this quadrant
                jumpphi[q] <- 2 * pi - atan((theta2.hat - lefty[q] + shift) / (leftx[q] - theta1.hat))      # unique to this quadrant
              }
              else if (cr[max(which(crlist$phi < rightphi[q])) - 2, 2] > cr[min(which(crlist$phi > leftphi[q])) + 2, 2]) {
                gaptype[q] <- "-"
                upR <- intersect(which(cr[, 1] > rightx[q]), which(cr[, 2] > righty[q]))   # points to upper-right quad of 'right' inaccess region border point
                shiftmaxindex <- min(intersect(upR, which(cr[, 2] < theta2.hat)))                           # unique to this quadrant
                shift <- jumpshift[q] * (cr[shiftmaxindex, 1] - rightx[q])                                  # unique to this quadrant
                jumpphi[q] <- 2 * pi - atan((theta2.hat - righty[q]) / ((rightx[q] + shift) - theta1.hat))  # unique to this quadrant
              }

              jumpxy[q,] <- crsolve(samp, cen, alpha = alpha + jumpuphill[q], mle.list = mle.list, phi = jumpphi[q])
              #points(cr[shiftmaxindex, 1], cr[shiftmaxindex, 2], pch = 16, col = "purple")
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
            if (!silent) {
              warning("alternate-centerpoint(s) used to repair plot regions inaccessible via a radial angle from its MLE")
              #message("alternate-centerpoint(s) used to repair plot regions inaccessible via a radial angle from its MLE")
            }
          }                                    # end if xrepair == 1
          else {                               # all iterations following 1st pass must re-assume repair parameters
            # each column in repairinfo represents a quadrant (with respect to the MLE)
            # "left" items are indicative of the left point from the perspective of the MLE (higher phi value); right is lower phi value
            #print("Loading repair info")
            leftphi <- repairinfo[1,]     # phi value to left point
            leftx <- repairinfo[2,]       # x value of "left" point (with respect to MLE)
            lefty <- repairinfo[3,]       # y value of "left" point (with respect to MLE)
            rightphi <- repairinfo[4,]
            rightx <- repairinfo[5,]
            righty <- repairinfo[6,]
            jumpphi <- repairinfo[7,]     # angle from MLE to locate jump-center
            jumpx <- repairinfo[8,]       # x value of jump center
            jumpy <- repairinfo[9,]       # y value of jump center
            phi1 <- repairinfo[10,]       # angle from jump-center to inaccessible region border point 1 (smaller phi angle)
            phi2 <- repairinfo[11,]       # angle from jump-center to inaccessible region border point 2 (larger phi angle)
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
          if ((length(phinewstore[phinewstore <= pi / 2] > 0)) && (done[q] != TRUE)) {
            #message("...repairing (Quad I relative to MLE)")
            repairinfo[12, q] <- TRUE      # annotate this quad is done

            # store the "new" centerpoint as theta1.hat & thata2.hat (albeit not an MLE; stays consistent with existing code)
            # maintain the MLE log-likelihood value as mleLLvalues because it is still needed in llrsolve
            # run smoothing search algorithm from jump point, and trim results to region of interest (per phi1, phi2 values)
            # note: since multiple roots available per phi, also check that y-values are not outside region of interest
            jinfo <- list("theta1.hat" = jumpxy[q, 1], "theta2.hat" = jumpxy[q, 2], "mleLLvalue" = mle.list$mleLLvalue)
            jpoints <- crsmooth(maxdeg = maxdeg, samp = mydata, cen = cen, alpha = alpha, mle.list = jinfo, ellipse_n = ellipse_n,
                                xrepair = xrepair, phinewstore = phinewstore, repairinfo = repairinfo, jumpxy = jumpxy, repairpass = TRUE,
                                repaircrlist = crlist, repairq = q)
            #          points(jpoints$x, jpoints$y, col = "purple", pch = 16, cex = 0.5)        # uncomment in order to plot jump-center repair-points

            # identify phi angle corresponding to MLE to remain consistent
            # to get phi value in range [0, 2pi) need to adjust off tan which computes from -pi/2 to pi/2:
            # if (jpoints$y - theta2.hat) is positive, add pi/2
            # if (jpoints$y - theta2.hat) is negative, add 3pi/2
            phiactual <- (pi - (pi / 2) * sign(jpoints$y - theta2.hat)) +
              atan(-(jpoints$x - theta1.hat) / (jpoints$y - theta2.hat))

            # algorithm for identifying insertafter index and keepindex adopted Nov 2018 (more robust) is:
            # - identify the midpoint of the edge bordering the inaccessible region,
            # - relative to this midpoint, identify CR points in its "furthest" quadrant (away from the MLE), and
            # - take the max or min of those points as appropriate for this quadrant (unique to current q value)
            # to identify one of the inaccessible region border points (& +/- 1 if appropriate to "jump" to other side)
            # - capture inaccessible region repair points based on which jpoints$phi angles occurs between the two boundary points
            midpt <- c(mean(c(leftx[q], rightx[q])), mean(c(lefty[q], righty[q])))
            #lowleft <- intersect(which(crlist$x < midpt[1]), which(crlist$y < midpt[2]))
            hiright <- intersect(which(crlist$x > midpt[1]), which(crlist$y > midpt[2]))
            if (gaptype[q] == "|") {
              insertafter <- min(hiright) - 1
            }
            else if (gaptype[q] == "-") {
              insertafter <- max(hiright)
            }
            # identify boundary points produced using jump-center that are relevant to keep (occur in the previously inaccessible region)
            # conditional on what side of the inaccessible region border the uncharted region lies
            angle2bp <- rep(0, 2)     # initiatlize
            for (bpoint in 1:2) {
              if (bpoint == 1) {
                bp <- c(rightx[q], righty[q])
              }
              else if (bpoint == 2) {
                bp <- c(leftx[q], lefty[q])
              }
              if ((jumpx[q] < bp[1]) && (jumpy[q] < bp[2])) {        # bp is quad 1 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpy[q] - bp[2]) / (jumpx[q] - bp[1])))
              }
              else if ((jumpx[q] > bp[1]) && (jumpy[q] < bp[2])) {   # bp is quad 2 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpx[q] - bp[1]) / (jumpy[q] - bp[2]))) + pi / 2
              }
              else if ((jumpx[q] > bp[1]) && (jumpy[q] > bp[2])) {   # bp is quad 3 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpy[q] - bp[2]) / (jumpx[q] - bp[1]))) + pi
              }
              else if ((jumpx[q] < bp[1]) && (jumpy[q] > bp[2])) {   # bp is quad 4 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpx[q] - bp[1]) / (jumpy[q] - bp[2]))) + 3 * pi / 2
              }
            }
            # print(format(angle2bp, digits = 9))
            # capture points in inaccessible region, whose jpoints$phi angle occurs between the two boundary points (angle2bp[1], angle2bp[2])
            if (angle2bp[1] < angle2bp[2]) {
              #keepindex <- intersect(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2]))
              keepindex <- intersect(
                intersect(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2])),  # w/in range-fan of jump-center bounded by infeasible region border points
                union(which(jpoints$x > theta1.hat), which(jpoints$y > theta2.hat)))            # ID points NOT in diagonally opposite quadrant (for q1, q4/1/2 repair pts are valid)
            }
            else if (angle2bp[1] > angle2bp[2]) {    # spans 0 degree (or 2pi radians) region
              #keepindex <- c(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2]))
              keepindex <- intersect(
                c(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2])),          # w/in range-fan of jump-center bounded by infeasible region border points
                union(which(jpoints$x > theta1.hat), which(jpoints$y > theta2.hat)))            # ID points NOT in diagonally opposite quadrant (for q1, q4/1/2 repair pts are valid)
            }

            # (previous technique:)
            #if (gaptype[q] == "|") {
            #  keepindex <- which(((jpoints$phi < phi1[q]) || (jpoints$phi > phi2[q])) & (jpoints$x > rightx[q]) & (phiactual < rightphi[q]))   # unique to this quad
            #}
            #else if (gaptype[q] == "-") {
            #  keepindex <- which((jpoints$phi > phi1[q]) & (jpoints$phi < phi2[q]) & (phiactual > leftphi[q]))      # unique to this quad
            #}
            #if (gaptype[q] == "|") {
            #   keepindex <- which(((jpoints$phi < phi1[q]) || (jpoints$phi > phi2[q])) & (jpoints$x > rightx[q]) & (phiactual < rightphi[q]))   # unique to this quad
            #}
            #else if (gaptype[q] == "-") {
            #  keepindex <- which((jpoints$phi > phi1[q]) & (jpoints$phi < phi2[q]) & (phiactual > leftphi[q]))      # unique to this quad
            #   keepindex <- which((jpoints$phi > phi1[q]) && (jpoints$x > leftx[q]) && (phiactual > rightphi[q]))    # unique to this quad
            #}

            jkeep <- list("x" = jpoints$x[keepindex],
                          "y" = jpoints$y[keepindex],
                          "phi" = jpoints$phi[keepindex] )
            #points(jkeep$x, jkeep$y, pch = 16, col = "yellow", cex = 0.5)

            # update phiactual angles w.r.t. MLE to represent only jump-center points kept (jkeep variable)
            phiactual <- atan((jkeep$y - theta2.hat) / (jkeep$x - theta1.hat))    # unique to this quad

            # insert additional CR boundary points into existing list
            # sequence angles being kept in the proper order
            go <- c(which(jkeep$phi > phi1[q]), which(jkeep$phi <= phi1[q]))
            # ensure new points integrate into the combined confidence region plot at right location
            # these steps are unique to this quadrant
            #if (gaptype[q] == "|") {
            #options <- intersect(which(crlist$y < jumpxy[q, 2]) & which(crlist$phi <= pi / 2))                        # candidates for insertion after are below the jump-center
            #insertafter <- which(crlist$phi == min(crlist$phi[options])) - 1  # "new" points fall before the lowest phi value among those options
            #}
            #else if (gaptype[q] == "-") {
            #  options <- intersect(which(crlist$y > jumpxy[q, 2]), which(crlist$phi <= pi / 2))                         # candidates for insertion after are above the jump-center
            #  insertafter <- which(crlist$phi == max(crlist$phi[options]))      # "new" points fall after the highest phi value among those options
            #}
            crlist$x <- append(crlist$x, jkeep$x[go], after = insertafter)
            crlist$y <- append(crlist$y, jkeep$y[go], after = insertafter)
            crlist$phi <- append(crlist$phi, phiactual[go], after = insertafter)

          }

          #####################################################
          # Quadrant II (with respect to MLE) repairs
          q <- 2        # quad II
          if ((length(phinewstore[(phinewstore > pi / 2) & (phinewstore <= pi)]) > 0) && (done[q] != TRUE)) {
            #message("...repairing (Quad II relative to MLE)")
            repairinfo[12, q] <- TRUE      # annotate this quad is done
            #          #gaptype[q] <- "-"              # assumed b/c no "|" cases currently known, requires modification similar to Quad III otherwise

            # store the "new" centerpoint as theta1.hat & thata2.hat (albeit not an MLE; stays consistent with existing code)
            # maintain the MLE log-likelihood value as mleLLvalues because it is still needed in llrsolve
            # run smoothing search algorithm from jump point, and trim results to region of interest (per phi1, phi2 values)
            # note: since multiple roots available per phi, also check that y-values are not outside region of interest
            jinfo <- list("theta1.hat" = jumpxy[q, 1], "theta2.hat" = jumpxy[q, 2], "mleLLvalue" = mle.list$mleLLvalue)
            jpoints <- crsmooth(maxdeg = maxdeg, samp = mydata, cen = cen, alpha = alpha, mle.list = jinfo, ellipse_n = ellipse_n,
                                xrepair = xrepair, phinewstore = phinewstore, repairinfo = repairinfo, jumpxy = jumpxy, repairpass = TRUE,
                                repaircrlist = crlist, repairq = q)

            # algorithm for identifying insertafter index and keepindex adopted Nov 2018 (more robust) is:
            # - identify the midpoint of the edge bordering the inaccessible region,
            # - relative to this midpoint, identify CR points in its "furthest" quadrant (away from the MLE), and
            # - take the max or min of those points as appropriate for this quadrant (unique to current q value)
            # to identify one of the inaccessible region border points (& +/- 1 if appropriate to "jump" to other side)
            # - capture inaccessible region repair points based on which jpoints$phi angles occurs between the two boundary points
            midpt <- c(mean(c(leftx[q], rightx[q])), mean(c(lefty[q], righty[q])))
            hileft <- intersect(which(crlist$x < midpt[1]), which(crlist$y > midpt[2]))
            if (gaptype[q] == "|") {
              insertafter <- max(hileft)
            }
            else if (gaptype[q] == "-") {
              insertafter <- min(hileft) - 1
            }
            # identify boundary points produced using jump-center that are relevant to keep (occur in the previously inaccessible region)
            # conditional on what side of the inaccessible region border the uncharted region lies
            angle2bp <- rep(0, 2)     # initiatlize
            for (bpoint in 1:2) {
              if (bpoint == 1) {
                bp <- c(rightx[q], righty[q])
              }
              else if (bpoint == 2) {
                bp <- c(leftx[q], lefty[q])
              }
              if ((jumpx[q] < bp[1]) && (jumpy[q] < bp[2])) {        # bp is quad 1 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpy[q] - bp[2]) / (jumpx[q] - bp[1])))
              }
              else if ((jumpx[q] > bp[1]) && (jumpy[q] < bp[2])) {   # bp is quad 2 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpx[q] - bp[1]) / (jumpy[q] - bp[2]))) + pi / 2
              }
              else if ((jumpx[q] > bp[1]) && (jumpy[q] > bp[2])) {   # bp is quad 3 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpy[q] - bp[2]) / (jumpx[q] - bp[1]))) + pi
              }
              else if ((jumpx[q] < bp[1]) && (jumpy[q] > bp[2])) {   # bp is quad 4 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpx[q] - bp[1]) / (jumpy[q] - bp[2]))) + 3 * pi / 2
              }
            }
            #print(format(angle2bp, digits = 9))
            # capture points in inaccessible region, whose jpoints$phi angle occurs between the two boundary points (angle2bp[1], angle2bp[2])
            if (angle2bp[1] < angle2bp[2]) {
              #keepindex <- intersect(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2]))
              keepindex <- intersect(
                intersect(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2])),  # w/in range-fan of jump-center bounded by infeasible region border points
                union(which(jpoints$x < theta1.hat), which(jpoints$y > theta2.hat)))            # ID points NOT in diagonally opposite quadrant (for q2, q1/2/3 repair pts are valid)
            }
            else if (angle2bp[1] > angle2bp[2]) {    # spans 0 degree (or 2pi radians) region
              #keepindex <- c(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2]))
              keepindex <- intersect(
                c(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2])),  # w/in range-fan of jump-center bounded by infeasible region border points
                union(which(jpoints$x < theta1.hat), which(jpoints$y > theta2.hat)))    # ID points NOT in diagonally opposite quadrant (for q2, q1/2/3 repair pts are valid)
            }
            #keepindex <- which((jpoints$phi < phi1[q] | jpoints$phi > phi2[q]) & jpoints$y > righty[q])
            #points(crlist$x[hileft], crlist$y[hileft], col = "green", pch = 16)
            #points(crlist$x[insertafter], crlist$y[insertafter], col = "red", pch = 16)
            #points(jpoints$x[keepindex], jpoints$y[keepindex], col = "yellow", pch = 16)

            # identify boundary points produced using jump-center that are relevant to keep (occur in the previously inaccessible region)
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
            #options <- which(crlist$y > jumpxy[q, 2])                         # candidates for insertion after are above the jump-center
            #insertafter <- which(crlist$phi == min(crlist$phi[options])) - 1  # "new" points fall before the lowest phi value among those options
            crlist$x <- append(crlist$x, jkeep$x[go], after = insertafter)
            crlist$y <- append(crlist$y, jkeep$y[go], after = insertafter)
            crlist$phi <- append(crlist$phi, phiactual[go], after = insertafter)

          }

          #####################################################
          # Quadrant III (with respect to MLE) repairs
          q <- 3
          if ((length(phinewstore[(phinewstore > pi) & (phinewstore <= 3 * pi / 2)]) > 0) && (done[q] != TRUE)) {
            #message("...repairing (Quad III relative to MLE)")
            repairinfo[12, q] <- TRUE      # annotate this quad is done

            # store the "new" centerpoint as theta1.hat & thata2.hat (albeit not an MLE; stays consistent with existing code)
            # maintain the MLE log-likelihood value as mleLLvalues because it is still needed in llrsolve
            # run smoothing search algorithm from jump point, and trim results to region of interest (per phi1, phi2 values)
            # note: since multiple roots available per phi, also check that y-values are not outside region of interest
            jinfo <- list("theta1.hat" = jumpxy[q, 1], "theta2.hat" = jumpxy[q, 2], "mleLLvalue" = mle.list$mleLLvalue)
            jpoints <- crsmooth(maxdeg = maxdeg, samp = mydata, cen = cen, alpha = alpha, mle.list = jinfo, ellipse_n = ellipse_n,
                                xrepair = xrepair, phinewstore = phinewstore, repairinfo = repairinfo, jumpxy = jumpxy, repairpass = TRUE,
                                repaircrlist = crlist, repairq = q)

            # algorithm for identifying insertafter index and keepindex adopted Nov 2018 (more robust) is:
            # - identify the midpoint of the edge bordering the inaccessible region,
            # - relative to this midpoint, identify CR points in its "furthest" quadrant (away from the MLE), and
            # - take the max or min of those points as appropriate for this quadrant (unique to current q value)
            # to identify one of the inaccessible region border points (& +/- 1 if appropriate to "jump" to other side)
            # - capture inaccessible region repair points based on which jpoints$phi angles occurs between the two boundary points
            midpt <- c(mean(c(leftx[q], rightx[q])), mean(c(lefty[q], righty[q])))
            lowleft <- intersect(which(crlist$x < midpt[1]), which(crlist$y < midpt[2]))
            if (gaptype[q] == "|") {
              insertafter <- min(lowleft) - 1
            }
            else if (gaptype[q] == "-") {
              insertafter <- max(lowleft)
            }
            # identify boundary points produced using jump-center that are relevant to keep (occur in the previously inaccessible region)
            # conditional on what side of the inaccessible region border the uncharted region lies
            angle2bp <- rep(0, 2)     # initiatlize
            for (bpoint in 1:2) {
              if (bpoint == 1) {
                bp <- c(rightx[q], righty[q])
              }
              else if (bpoint == 2) {
                bp <- c(leftx[q], lefty[q])
              }
              if ((jumpx[q] < bp[1]) && (jumpy[q] < bp[2])) {        # bp is quad 1 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpy[q] - bp[2]) / (jumpx[q] - bp[1])))
              }
              else if ((jumpx[q] > bp[1]) && (jumpy[q] < bp[2])) {   # bp is quad 2 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpx[q] - bp[1]) / (jumpy[q] - bp[2]))) + pi / 2
              }
              else if ((jumpx[q] > bp[1]) && (jumpy[q] > bp[2])) {   # bp is quad 3 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpy[q] - bp[2]) / (jumpx[q] - bp[1]))) + pi
              }
              else if ((jumpx[q] < bp[1]) && (jumpy[q] > bp[2])) {   # bp is quad 4 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpx[q] - bp[1]) / (jumpy[q] - bp[2]))) + 3 * pi / 2
              }
            }
            #print(format(angle2bp, digits = 9))
            # capture points in inaccessible region, whose jpoints$phi angle occurs between the two boundary points (angle2bp[1], angle2bp[2])
            if (angle2bp[1] < angle2bp[2]) {
              #keepindex <- intersect(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2]))
              keepindex <- intersect(
                intersect(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2])),  # w/in range-fan of jump-center bounded by infeasible region border points
                union(which(jpoints$x < theta1.hat), which(jpoints$y < theta2.hat)))            # ID points NOT in diagonally opposite quadrant (for q3, q2/3/4 repair pts are valid)
            }
            else if (angle2bp[1] > angle2bp[2]) {    # spans 0 degree (or 2pi radians) region
              #keepindex <- c(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2]))
              keepindex <- intersect(
                c(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2])),          # w/in range-fan of jump-center bounded by infeasible region border points
                union(which(jpoints$x < theta1.hat), which(jpoints$y < theta2.hat)))            # ID points NOT in diagonally opposite quadrant (for q3, q2/3/4 repair pts are valid)
            }
            #points(crlist$x[lowleft], crlist$y[lowleft], col = "green", pch = 16)
            #points(crlist$x[insertafter], crlist$y[insertafter], col = "red", pch = 16)
            #points(jpoints$x[keepindex], jpoints$y[keepindex], col = "yellow", pch = 16)

            # identify phi angle corresponding to MLE to remain consistent
            # to get phi value in range [0, 2pi) need to adjust off tan which computes from -pi/2 to pi/2:
            # if (jpoints$y - theta2.hat) is positive, add pi/2
            # if (jpoints$y - theta2.hat) is negative, add 3pi/2
            phiactual <- (pi - (pi / 2) * sign(jpoints$y - theta2.hat)) +
              atan(-(jpoints$x - theta1.hat) / (jpoints$y - theta2.hat))

            # identify boundary points produced using jump-center that are relevant to keep (occur in the previously inaccessible region)
            # conditional on what side of the inaccessible region border the uncharted region lies
            #if (gaptype[q] == "|") {
            # keepindex <- which((jpoints$phi < phi1[q]) & (jpoints$x < rightx[q]) & (phiactual < rightphi[q]))
            #}
            #else if (gaptype[q] == "-") {
            #  keepindex <- which((jpoints$phi > phi1[q]) & (jpoints$x < leftx[q]) & (phiactual > rightphi[q]))
            #}
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
            # if (gaptype[q] == "|") {
            #   options <- which(crlist$x < jumpxy[q, 1])                         # candidates for insertion after are left of the jump-center
            #   insertafter <- which(crlist$phi == min(crlist$phi[options])) - 1  # "new" points fall before the lowest phi value among those options
            # }
            # else if (gaptype[q] == "-") {
            #   options <- which(crlist$y < jumpxy[q, 2])                         # candidates for insertion after are below the jump-center
            #   insertafter <- which(crlist$phi == max(crlist$phi[options]))      # "new" points fall after the highest phi value among those options
            # }
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

            # store the "new" centerpoint as theta1.hat & thata2.hat (albeit not an MLE; stays consistent with existing code)
            # maintain the MLE log-likelihood value as mleLLvalues because it is still needed in llrsolve
            # run smoothing search algorithm from jump point, and trim results to region of interest (per phi1, phi2 values)
            # note: since multiple roots available per phi, also check that y-values are not outside region of interest
            jinfo <- list("theta1.hat" = jumpxy[q, 1], "theta2.hat" = jumpxy[q, 2], "mleLLvalue" = mle.list$mleLLvalue)
            jpoints <- crsmooth(maxdeg = maxdeg, samp = mydata, cen = cen, alpha = alpha, mle.list = jinfo, ellipse_n = ellipse_n,
                                xrepair = xrepair, phinewstore = phinewstore, repairinfo = repairinfo, jumpxy = jumpxy, repairpass = TRUE,
                                repaircrlist = crlist, repairq = q)

            # algorithm for identifying insertafter index and keepindex adopted Nov 2018 (more robust) is:
            # - identify the midpoint of the edge bordering the inaccessible region,
            # - relative to this midpoint, identify CR points in its "furthest" quadrant (away from the MLE), and
            # - take the max or min of those points as appropriate for this quadrant (unique to current q value)
            # to identify one of the inaccessible region border points (& +/- 1 if appropriate to "jump" to other side)
            # - capture inaccessible region repair points based on which jpoints$phi angles occurs between the two boundary points
            midpt <- c(mean(c(leftx[q], rightx[q])), mean(c(lefty[q], righty[q])))
            lowright <- intersect(which(crlist$x > midpt[1]), which(crlist$y < midpt[2]))
            #          points(crlist$x[lowright], crlist$y[lowright], col = "green", pch = 16)
            if (gaptype[q] == "|") {
              insertafter <- max(lowright)
            }
            else if (gaptype[q] == "-") {
              insertafter <- min(lowright) - 1
            }
            #          points(crlist$x[insertafter], crlist$y[insertafter], col = "red", pch = 16)
            # identify boundary points produced using jump-center that are relevant to keep (occur in the previously inaccessible region)
            # conditional on what side of the inaccessible region border the uncharted region lies
            angle2bp <- rep(0, 2)     # initiatlize
            for (bpoint in 1:2) {
              if (bpoint == 1) {
                bp <- c(rightx[q], righty[q])
              }
              else if (bpoint == 2) {
                bp <- c(leftx[q], lefty[q])
              }
              if ((jumpx[q] < bp[1]) && (jumpy[q] < bp[2])) {        # bp is quad 1 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpy[q] - bp[2]) / (jumpx[q] - bp[1])))
              }
              else if ((jumpx[q] > bp[1]) && (jumpy[q] < bp[2])) {   # bp is quad 2 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpx[q] - bp[1]) / (jumpy[q] - bp[2]))) + pi / 2
              }
              else if ((jumpx[q] > bp[1]) && (jumpy[q] > bp[2])) {   # bp is quad 3 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpy[q] - bp[2]) / (jumpx[q] - bp[1]))) + pi
              }
              else if ((jumpx[q] < bp[1]) && (jumpy[q] > bp[2])) {   # bp is quad 4 relative to jumpxy
                angle2bp[bpoint] <- atan(abs((jumpx[q] - bp[1]) / (jumpy[q] - bp[2]))) + 3 * pi / 2
              }
            }
            #print(format(angle2bp, digits = 9))
            # capture points in inaccessible region, whose jpoints$phi angle occurs between the two boundary points (angle2bp[1], angle2bp[2])
            if (angle2bp[1] < angle2bp[2]) {
              #keepindex <- intersect(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2]))
              keepindex <- intersect(
                intersect(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2])),  # w/in range-fan of jump-center bounded by infeasible region border points
                union(which(jpoints$x > theta1.hat), which(jpoints$y < theta2.hat)))            # ID points NOT in diagonally opposite quadrant (for q4, q3/4/1 repair pts are valid)
            }
            else if (angle2bp[1] > angle2bp[2]) {    # spans 0 degree (or 2pi radians) region
              #keepindex <- c(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2]))
              keepindex <- intersect(
                c(which(jpoints$phi > angle2bp[1]), which(jpoints$phi < angle2bp[2])),          # w/in range-fan of jump-center bounded by infeasible region border points
                union(which(jpoints$x > theta1.hat), which(jpoints$y < theta2.hat)))            # ID points NOT in diagonally opposite quadrant (for q4, q3/4/1 repair pts are valid)
            }
            #          points(jpoints$x[keepindex], jpoints$y[keepindex], col = "yellow", pch = 16)

            # identify boundary points produced using jump-center that are relevant to keep (occur in the previously inaccessible region)
            #keepindex <- which((jpoints$phi < phi1[q] | jpoints$phi > phi2[q]) & jpoints$x > leftx[q])
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
            #options <- which(crlist$x > jumpxy[q, 1])                         # candidates for insertion after are right of jump-center
            #insertafter <- which(crlist$phi == max(crlist$phi[options]))      # "new" points fall after the highest phi value among those options
            crlist$x <- append(crlist$x, jkeep$x[go], after = insertafter)
            crlist$y <- append(crlist$y, jkeep$y[go], after = insertafter)
            crlist$phi <- append(crlist$phi, phiactual[go], after = insertafter)

          }
          #print(jumpxy)   # print jump-center locations corresponding to each respective quadrant (quadrants are relative to MLE)
        }                  # end repairpass
        #print("status of maxcurrent, maxrad, count:")
        #print(c(maxcurrent, maxrad, count))
        maxcurrent <- maxrad
        if (!repair) {
          #warn_maxdeg <- TRUE
          warning("max iteration tolerance hit; working confidence region plot is shown, however, maxdeg constraint is not met")
          #print("WARNING: max iteration tolerance hit; maxdeg constraint not met, however, working plot is shown.")
          #print("...this is attributable to either regions inaccessible via a radial azimuth from its MLE, or a uniroot tol argument too small.")
        }
        #else if ((repair) && (count == maxcount)) {
        #  warn_maxdeg <- TRUE
        #  warning("max iteration tolerance hit; working confidence region plot is shown, however, maxdeg constraint is not met")
        #}
      }                # end (if (count > maxcount) && repairpass != TRUE)

    }                # end (while (maxcurrent > maxrad))

    # ##### check if any uniroot solutions returned ~0 (i.e. identical location to MLE), and discard if necessary
    # # note: a numeric difficulty sometimes arises when phi angle is in {0, pi/2, pi, 3pi/2, 2pi} causing this error
    invalid_pt <- intersect(which(crlist$x == theta1.hat), which(crlist$y == theta2.hat))
    #print(invalid_pt)
    if (length(invalid_pt) != 0) {
      crlist$x <- crlist$x[-c(invalid_pt)]
      crlist$y <- crlist$y[-c(invalid_pt)]
      crlist$phi <- crlist$phi[-c(invalid_pt)]
      #print("REMOVED INDEX #(s):")
      #print(invalid_pt)
    }
    crlist$repair <- repairinfo

    invisible(crlist)        # return results of crsmooth function call
  }                          # end crsmooth


  ######################################################################################

  # nomadjust -------------------------------------------------------------------
  # Adjusts nominal value (1 - alpha) to compensate for negative bias inherent in small
  # sample sizes.  Uses tables assembled via 20 million Monte Carlo simulation replications
  # of the conf coversim() function.  Available for a limited range of alpha values
  # (roughly <= 0.25) and distributions (distn suffixes: weibull, llogis, norm).
  nomadjust = function(distn = distn, target_alpha = alpha, this_n = length(dataset)){

    # confirm distn is supported with MC sim table of results
    if (!(distn %in% c("weibull", "llogis", "norm"))) {
      return(invisible(alpha))     # stop function execution and return alpha value unchanged
      #print("cannot correct alpha for exact coverage; this distribution is not supported")
    }

    # when n >= 4 and n <= 50 use MC simulation table values
    # Monte Carlo simulation results:
    # Weibull MC results (coverage over 20 million replications)
    w_ac_a80 <- c(0.7035459, 0.72582845, 0.7398501, 0.7491972, 0.7560088, 0.76158625, 0.7655136, 0.77743695, 0.7833055, 0.786756, 0.7887199, 0.7902545, 0.79157055, 0.7925772, 0.79311415)
    w_ac_a85 <- c(0.7613011, 0.78229525, 0.79511275, 0.8039439, 0.81035855, 0.8150749, 0.81904215, 0.8297237, 0.83495145, 0.8380563, 0.84012485, 0.8414208, 0.84241, 0.8432222, 0.8441716)
    w_ac_a90 <- c(0.8239086, 0.84278945, 0.8538984, 0.86167525, 0.8670003, 0.87113315, 0.8742291, 0.8834281, 0.88758865, 0.89024785, 0.89195575, 0.89302095, 0.893891, 0.89452935, 0.89508435)
    w_ac_a95 <- c(0.8952844, 0.90955215, 0.91806025, 0.9236672, 0.9274081, 0.9303495, 0.9326231, 0.9387996, 0.9418163, 0.9434903, 0.94459105, 0.9454118, 0.9460034, 0.9464411, 0.9468038)
    w_ac_a99 <- c(0.96813145, 0.9749135, 0.9783814, 0.9805481, 0.9821339, 0.98329625, 0.98412775, 0.98633785, 0.98734695, 0.9879734, 0.98832405, 0.98857495, 0.98873545, 0.9888664, 0.9889654)
    mc_weibull <- c(w_ac_a80, w_ac_a85, w_ac_a90, w_ac_a95, w_ac_a99)
    # normal MC results (coverage over 20 million replications)
    n_ac_a80 <- c(0.7067501, 0.7287328, 0.74250665, 0.7517541, 0.75846435, 0.763348, 0.7673942, 0.7788696, 0.78439315, 0.7875951, 0.7897189, 0.79112585, 0.7921908, 0.79308715, 0.79377625)
    n_ac_a85 <- c(0.7640325, 0.78500925, 0.79780445, 0.80634185, 0.8126448, 0.8170626, 0.82075875, 0.83109545, 0.8361931, 0.83912065, 0.84089535, 0.84204465, 0.8431295, 0.843756, 0.8445106)
    n_ac_a90 <- c(0.82648165, 0.8450147, 0.85624205, 0.86372445, 0.86887655, 0.87280995, 0.8757599, 0.8843581, 0.888553, 0.8909364, 0.892516, 0.89367145, 0.89439125, 0.89508815, 0.89557665)
    n_ac_a95 <- c(0.89732245, 0.9115436, 0.91965495, 0.92503545, 0.92878295, 0.93150815, 0.93374385, 0.93971525, 0.9424231, 0.94405625, 0.94503705, 0.9457903, 0.94634395, 0.94663605, 0.9470519)
    n_ac_a99 <- c(0.9695351, 0.9757233, 0.9792197, 0.98130455, 0.98279225, 0.98369775, 0.98447425, 0.98661625, 0.98757105, 0.9881173, 0.98843675, 0.9887229, 0.9888389, 0.98894865, 0.98910405)
    mc_norm <- c(n_ac_a80, n_ac_a85, n_ac_a90, n_ac_a95, n_ac_a99)
    # log-logistic MC results (coverage over 20 million replications)
    ll_ac_a80 <- c(0.7185017, 0.73881295, 0.75085615, 0.75882015, 0.76458925, 0.76908515, 0.77242395, 0.7823134, 0.7866679, 0.7893375, 0.791194, 0.79257775, 0.79349415, 0.79407755, 0.79468935)
    ll_ac_a85 <- c(0.77326175, 0.7933564, 0.8049081, 0.81234855, 0.8177637, 0.821871, 0.82502725, 0.8341202, 0.83820135, 0.8404987, 0.8422028, 0.84331675, 0.84418595, 0.84473095, 0.84529455)
    ll_ac_a90 <- c(0.830719, 0.85042615, 0.8613384, 0.86850905, 0.87290365, 0.8764745, 0.8790613, 0.8868219, 0.8903339, 0.8921646, 0.89358425, 0.8945531, 0.89519575, 0.8957709, 0.89627705)
    ll_ac_a95 <- c(0, 0.9093955, 0.921372, 0.9271494, 0.930939, 0.93366255, 0.9356201, 0.9410063, 0.94345935, 0.94487485, 0.94571615, 0.94632005, 0.946788, 0.9471269, 0.9474715)
    ll_ac_a99 <- c(0, 0.59821, 0.898577, 0.969455, 0.980772, 0.9834805, 0.98458495, 0.98695035, 0.98782425, 0.98833965, 0.98859545, 0.98882365, 0.9889766, 0.9891222, 0.98919195)
    # zero-out invalid entries having prohibitive MC sim errors
    ll_ac_a90[1] <- 0
    ll_ac_a95[1:3] <- rep(0, 3)
    ll_ac_a99[1:5] <- rep(0, 5)
    mc_llogis <- c(ll_ac_a80, ll_ac_a85, ll_ac_a90, ll_ac_a95, ll_ac_a99)

    # parameters and labels corresponding to MC sim results
    lab_alpha <- c("0.8", "0.85", "0.9", "0.95", "0.99")
    lab_n <- c("4", "5", "6", "7", "8", "9", "10", "15", "20", "25", "30", "35", "40", "45", "50")
    ci_type <- c("Clopper-Pearson", "Wald", "Wilson-Score", "Jeffreys", "Agresti-Coull", "Arcsine", "Blaker")
    alpha_seq <- rep(c(seq(0.2, 0.05, by = -0.05), 0.01), each = length(lab_n))
    n_all <- c(4:10, seq(15, 50, by = 5))
    n_seq <- rep(n_all, length(lab_alpha))
    alpha_matrix <- matrix(alpha_seq, nrow = 5, byrow = TRUE, dimnames = list(lab_alpha, lab_n))
    n_matrix <- matrix(n_seq, nrow = length(lab_alpha), byrow = TRUE, dimnames = list(lab_alpha, lab_n))

    # assemble MC sim summaries in matrix form
    weibull_matrix <- matrix(mc_weibull, ncol = length(lab_n), byrow = TRUE, dimnames = list(lab_alpha, lab_n))
    llogis_matrix <- matrix(mc_llogis, ncol = length(lab_n), byrow = TRUE, dimnames = list(lab_alpha, lab_n))
    norm_matrix <- matrix(mc_norm, ncol = length(lab_n), byrow = TRUE, dimnames = list(lab_alpha, lab_n))

    # create confidence intervals
    mcreps <- rep(2 * 10 ^ 7, length(mc_weibull))     # each value represents 20 million replications
    ciw_max <- matrix(rep(0, 3 * length(ci_type)), ncol = 3, dimnames = list(ci_type, c("weibull", "llogis", "norm")))
    for (i in 1:3) {                                  # assess each distn (weibull, llogis, norm)
      if (i == 1) {
        use_this <- mc_weibull
        #mcreps <- mcreps_weibull
        this_distn <- "Weibull"
      }
      else if (i == 2) {
        use_this <- mc_llogis
        #mcreps <- mcreps_llogis
        this_distn <- "log-logistic"
      }
      else if (i == 3) {
        use_this <- mc_norm
        #mcreps <- mcreps_norm
        this_distn <- "normal"
      }
      mcnresults <- mcreps * use_this
      mc_ci <- matrix(rep(0, 2 * length(mcnresults)), ncol = 2)
      for (j in 1:length(ci_type)) {
        ci_t <- ci_type[j]
        a <- 0.05
        for (k in 1:length(mcnresults)) {
          if (as.integer(mcnresults[k]) != 0) {   # skip if no MC sim result recorded
            mc_ci[k, ] <- binomTest(mcreps[k], as.integer(mcnresults[k]), alpha = a, intervalType = ci_t)
          }
        }
        this_summary <- data.frame(cbind(n_seq, alpha_seq,
                                         mcnresults / mcreps, mc_ci,
                                         mcnresults, mcreps))
        colnames(this_summary) = c("n", "alpha", "coverage", "lower", "upper", "covered", "reps")
        this_summary <- this_summary[order(this_summary$n, this_summary$alpha),]   # sort by n, then alpha
        this_summary$reps[this_summary$coverage == 0] <- 0     # (no coversim data)
        this_summary[this_summary == 0] <- "NA"                # replace 0s with a blank (no coversim data)
        rownames(this_summary) <- 1:length(this_summary$n)
        if (i == 1) {
          summary_weibull <- this_summary
        }
        else if (i == 2) {
          summary_llogis <- this_summary
        }
        else if (i == 3) {
          summary_norm <- this_summary
        }
        ciw_max[j, i] <- max(mc_ci[,2] - mc_ci[,1])             # record max CI width
      }
      #if (i == 3) {           # print summary of maximum CI width and half-width among all parameterizations
      #  print("max CI widths:")
      #  print(ciw_max)
      #  print("max CI half-widths:")
      #  print(ciw_max / 2)
      #  print(paste0("Weibull min(max(CI width)): ", ci_type[which(ciw_max[,1] == min(ciw_max[,1]))]))
      #  print(paste0("llogis  min(max(CI width)): ", ci_type[which(ciw_max[,2] == min(ciw_max[,2]))]))
      #  print(paste0("normal  min(max(CI width)): ", ci_type[which(ciw_max[,3] == min(ciw_max[,3]))]))
      #}
    }    # end of create confidence intervals section

    # linear interpolation between "5s" to expand tables from (10, 15, 20, ..., 50) to (10, 11, 12, ..., 50)
    # so "get" matrices can be formed for any n (note, ideally MC sim later run on ALL values rather than interpolation route)
    # assemble "expanded" tables via linear interpolation between MC sim values & assemble in data.frame format
    norm_expanded <- as.data.frame(norm_matrix)
    llogis_expanded<- as.data.frame(llogis_matrix)
    for (k in 1:3) {
      if (k == 1) {
        this_matrix <- weibull_matrix
      } else if (k == 2) {
        this_matrix <- llogis_matrix
      } else if (k == 3) {
        this_matrix <- norm_matrix
      }
      this_expanded <- data.frame("n" = as.numeric(lab_n), "a0.8" = this_matrix[1,], "a0.85" = this_matrix[2,],
                                  "a0.9" = this_matrix[3,], "a0.95" = this_matrix[4,], "a0.99" = this_matrix[5,])
      rownames(this_expanded) <- 1:nrow(this_expanded)
      for (i in 0:7) {
        thisrow <- which(this_expanded$n == (10 + i * 5))
        slp <- as.numeric((this_expanded[thisrow + 1, 2:6] - this_expanded[thisrow, 2:6]) / 5)
        for (j in 1:4) {
          n <- as.numeric(10 + i * 5 + j)
          newrow <- as.numeric(this_expanded[thisrow, 2:6]) + slp * j
          this_expanded[nrow(this_expanded) + 1,] = list("n" = n, "a0.8" = newrow[1], "a0.85" = newrow[2],
                                                         "a0.9" = newrow[3], "a0.95" = newrow[4],
                                                         "a0.99" = newrow[5])
        }
      }
      this_expanded <- this_expanded[order(this_expanded$n),]   # sort by n
      rownames(this_expanded) <- 1:length(this_expanded$n)
      # store results in distn_expanded named data.frame
      if (k == 1) {
        weibull_get <- this_expanded
      } else if (k == 2) {
        llogis_get <- this_expanded
      } else if (k == 3) {
        norm_get <- this_expanded
      }
    }

    # calculate adjustment to alpha value necessary to compensate for negative bias and acheive actual coverage probability = 1 - alpha
    if (distn == "weibull") {
      this_get <- weibull_get
    } else if (distn == "llogis") {
      this_get <- llogis_get
    } else if (distn == "norm") {
      this_get <- norm_get
    }
    row_id <- which(this_get$n == this_n)
    # ensure target_alpha value is acheivable via interpolation of "known" points
    if ((length(as.numeric(which(this_get$n == this_n))) == 0) ||   # n value outside of table values
        (1 - target_alpha) < min(replace(this_get[row_id, 2:6], this_get[row_id, 2:6] == 0, 1)) || # alpha outside table values
        (distn == "llogis") &&                                      # llogis restrictions due to MC sim errors
        ((this_n <= 4) ||
         ((this_n == 5) && (target_alpha < 0.15)) ||
         ((this_n == 6) && (target_alpha < 0.15)) ||
         ((this_n == 7) && (target_alpha < 0.10)) ||
         ((this_n == 8) && (target_alpha < 0.10)))) {
      print("cannot correct alpha for exact coverage; extrapolation outside of known coverage values is necessary")
      return(invisible(alpha))      # return alpha value unchanged
    } else {
      #print("adjusting alpha to compensate for coverage bias")
      if (((this_n %% 5) != 0 ) && (n > 10) && (!silent)) {
        print(paste0("estimating n = ", this_n, " via linear interpolation of n = ", this_n - (this_n %% 5),
                     " and n = ", this_n + 5 - (this_n %% 5), " Monte Carlo simulation results"))
      }
      indx2 <- which.max(this_get[row_id, 2:6] > (1 - target_alpha)) + 1  # ID index 'after' target (+1 compensates for 2:6 range)
      if (all(FALSE == (this_get[row_id, 2:6] > (1 - target_alpha)))) {   # > 0.99 entry
        indx2 <- 7
        slp <- (1 - this_get[row_id, indx2 - 1]) / 0.01
        a <- 0.99
      }
      else if (indx2 == 6) {    # between 0.95 and 0.99 entries
        slp <- (this_get[row_id, indx2] - this_get[row_id, indx2 - 1]) / 0.04
        a <- 0.95
      }
      else {                    # between 0.8 and 0.95 entries
        slp <- (this_get[row_id, indx2] - this_get[row_id, indx2 - 1]) / 0.05
        a <- c(0, seq(0.8, 0.95, by = 0.05))[indx2 - 1]
      }
      miss <- (1 - target_alpha) - this_get[row_id, indx2 - 1]    # delta y (in nominal vs actual coverage plot)
      add <- miss / slp                                                # delta x
      get_target <- a + add
      get_target3 <- round(get_target, digits = 3)
      #if (!silent) {
      #  print(paste0("using nominal coverage of ", get_target3, " (alpha = ", 1 - get_target3, ") to achieve an actual coverage ", 1 - target_alpha))
      #}
      return(invisible(1 - get_target))   # return adjusted alpha value
    }
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

  # make adjustment for alpha value to compensate for coverage bias
  if (exact == TRUE) {
    alpha_target <- alpha      # store original alpha value, which will be adjusted in nomadjust() call
    alpha <- nomadjust(distn = distn, target_alpha = alpha, this_n = length(dataset))
  }

  # defaults are set to:
  #jumpshift <- 0.5                   # % along available shift range for phi to locate jump point
  #jumpuphill <- min(alpha, 0.01)     # percentage uphill from CR boarder
  if (length(jumpshift) == 1) {
    jumpshift <- rep(jumpshift, 4)
  }
  if (length(jumpuphill) == 1) {
    jumpuphill <- rep(jumpuphill, 4)
  }

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
    ase <- try(asesolve(x = dataset, cen, theta1.hat, theta2.hat), silent = TRUE)

    # circumstance where ASE is unatainable include computational singular observed info matrix (ase returned as a
    # "character" string with error message) or NaN result (both filtered with "if" statement below)
    if (is.list(ase) && !is.nan(ase$theta1.ase) && !is.nan(ase$theta2.ase)) {
      phi <- angles(a = ase$theta1.ase, b = ase$theta2.ase, npoints = ellipse_n)
    }
    else  {        # if ASE unavailable, estimate aspect ratio from MLE
      message("ASE calc unsuccessful; using MLE to estimate aspect ratio")
      phi <- angles(a = theta1.hat, b = theta2.hat, npoints = ellipse_n)
    }

    cr <- crsolve(samp = dataset, cen = cen, alpha = alpha, mle.list = mle.list, phi = phi)
    crlist <- list("x" = cr[, 1], "y" = cr[, 2], "phi" = phi)
  }

  # assemble conf region points and "close" the gap between first & last points graphically
  cr <- cbind(crlist$x, crlist$y)
  cr <- rbind(cr, cr[1, ])

  if (xyswap) {
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

  if (showplot) {
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
        if (!xyswap) {
          xlab = xaxislabel[which(disttype == distn)]
        }
        else if (xyswap) {
          xlab = yaxislabel[which(disttype == distn)]
        }
      }
    }
    if (!is.expression(ylab)) {
      if (ylab == "") {
        if (!xyswap) {
          ylab = yaxislabel[which(disttype == distn)]
        }
        else if (xyswap) {
          ylab = xaxislabel[which(disttype == distn)]
        }
      }
    }

    # axis tick-marks based on mlelab and origin
    xlabs <- c(min(cr[, 1]), max(cr[, 1]))
    ylabs <- c(min(cr[, 2]), max(cr[, 2]))
    if (mlelab) {
      xlabs <- c(xlabs, theta1.hat)
      ylabs <- c(ylabs, theta2.hat)
    }
    if (origin) {
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
      if (!xyswap) {
        xlim <- c(xmin, max(cr[,1]))
      }
      else if (xyswap) {
        xlim <- c(ymin, max(cr[,2]))
      }
    }
    else {
      xlabs <- c(xlim, xlabs)
      xlabs <- xlabs[xlabs >= xlim[1] & xlabs <= xlim[2]]
    }
    if(is.null(ylim)) {
      if (!xyswap) {
        ylim <- c(ymin, max(cr[,2]))
      }
      else if (xyswap) {
        ylim <- c(xmin, max(cr[,1]))
      }
    }
    else {
      ylabs <- c(ylim, ylabs)
      ylabs <- ylabs[ylabs >= ylim[1] & ylabs <= ylim[2]]
    }

    # plot
    if (!xyswap) {
      plot(cr, xlab = xlab, ylab = ylab, main = main, ylim = ylim, xlim = xlim,
           axes = FALSE, type = 'l')
      if (is.null(sf)) {
        my.atx <- xlabs
        my.aty <- ylabs
      }
      else {
        my.atx <- unique(round(xlabs, sf[1]))
        my.aty <- unique(round(ylabs, sf[2]))
      }
      axis(side = 1, at = my.atx, las = xlas)
      axis(side = 2, at = my.aty, las = ylas)
      if (mlelab) {
        points(theta1.hat, theta2.hat, pch = 3)
      }
      if (pts) points(unique(cr), lwd = 0.65)
    }
    else if (xyswap) {
      crswap <- cbind(cr[,2], cr[,1])
      plot(crswap, xlab = xlab, ylab = ylab, main = main, ylim = ylim, xlim = xlim,
           axes = FALSE, type = 'l')
      if (is.null(sf)) {
        my.atx <- ylabs
        my.aty <- xlabs
      }
      else {
        my.atx <- unique(round(ylabs, sf[1]))
        my.aty <- unique(round(xlabs, sf[2]))
      }
      axis(side = 1, at = my.atx, las = xlas)
      axis(side = 2, at = my.aty, las = ylas)
      if (mlelab) {
        points(theta2.hat, theta1.hat, pch = 3)
      }
      if (pts) points(unique(crswap), lwd = 0.65)
    }
    # enable plotting beyond margins for any post-processing user add-ons
    # ...or disable
    #par(xpd = TRUE)
    par(xpd = FALSE)
  } # end of:  if (showplot) {...}

  #print("dataset was: ")
  #print(dataset)
  # returned values (if requested), otherwise only output to screen # boundary points.
  if ((length(dataset) <= 5) && (distn != "unif") && (!silent)) {
    warning("small sample size is ill-suited to invoke the asymptotic properties assumed by the confidence region plot")
    #print("WARNING: small sample sizes can yield irregular confidence region shapes that are unatainable using crplot.")
    #print("WARNING: small sample sizes violate the asymptotic properties assumption used to produce the confidence region.")
  }
  # re-build crlist with info requested to return
  # assemble info with return_crlistx; if requested to return info or jumpinfo, then set return_crlist <- return_crlistx
  #if ((info) || (jumpinfo)) {
  #print(paste0("Confidence region plot complete; made using ", length(phi)," boundary points."))
  if (distn == "weibull") {
    #print(paste0("MLE value is: (kappa.hat = ", theta1.hat, ", lambda.hat = ", theta2.hat,")"))
    return_crlistx <- list("kappa" = crlist$x, "lambda" = crlist$y, "phi" = crlist$phi,
                           "kappahat" = mle.list$theta1.hat, "lambdahat" = mle.list$theta2.hat)
  }
  else if (distn == "invgauss") {
    #print(paste0("MLE value is: (mu.hat = ", theta1.hat, ", lambda.hat = ", theta2.hat,")"))
    return_crlistx <- list("mu" = crlist$x, "lambda" = crlist$y, "phi" = crlist$phi,
                           "muhat" = mle.list$theta1.hat, "lambdahat" = mle.list$theta2.hat)
  }
  else if (distn == "norm") {
    #print(paste0("MLE value is: (mu.hat = ", theta1.hat, ", sigma.hat = ", theta2.hat,")"))
    return_crlistx <- list("mu" = crlist$x, "sigma" = crlist$y, "phi" = crlist$phi,
                           "muhat" = mle.list$theta1.hat, "sigmahat" = mle.list$theta2.hat)
  }
  else if (distn == "lnorm") {
    #print(paste0("MLE value is: (mu.hat = ", theta1.hat, ", sigma.hat = ", theta2.hat,")"))
    return_crlistx <- list("mu" = crlist$x, "sigma" = crlist$y, "phi" = crlist$phi,
                           "muhat" = mle.list$theta1.hat, "sigmahat" = mle.list$theta2.hat)
  }
  else if (distn == "logis") {
    #print(paste0("MLE value is: (mu.hat = ", theta1.hat, ", sigma.hat = ", theta2.hat,")"))
    return_crlistx <- list("mu" = crlist$x, "sigma" = crlist$y, "phi" = crlist$phi,
                           "muhat" = mle.list$theta1.hat, "sigmahat" = mle.list$theta2.hat)
  }
  else if (distn == "llogis") {
    #print(paste0("MLE value is: (lambda.hat = ", theta1.hat, ", kappa.hat = ", theta2.hat,")"))
    return_crlistx <- list("lambda" = crlist$x, "kappa" = crlist$y, "phi" = crlist$phi,
                           "lambdahat" = mle.list$theta1.hat, "kappahat" = mle.list$theta2.hat)
  }
  else if (distn == "gamma") {
    #print(paste0("MLE value is: (theta.hat = ", theta1.hat, ", kappa.hat = ", theta2.hat,")"))
    return_crlistx <- list("theta" = crlist$x, "kappa" = crlist$y, "phi" = crlist$phi,
                           "thetahat" = mle.list$theta1.hat, "kappahat" = mle.list$theta2.hat)
  }
  else if (distn == "unif") {
    #print(paste0("MLE value is: (a.hat = ", theta1.hat, ", b.hat = ", theta2.hat,")"))
    return_crlistx <- list("a" = crlist$x, "b" = crlist$y, "phi" = crlist$phi,
                           "ahat" = mle.list$theta1.hat, "bhat" = mle.list$theta2.hat)
  }
  else if (distn == "cauchy") {
    #print(paste0("MLE value is: (a.hat = ", theta1.hat, ", s.hat = ", theta2.hat,")"))
    return_crlistx <- list("a" = crlist$x, "s" = crlist$y, "phi" = crlist$phi,
                           "ahat" = mle.list$theta1.hat, "shat" = mle.list$theta2.hat)
  }
  #return(return_crlist)
  #}   # end of "if (info == TRUE)" conditional
  if ((info) || (jumpinfo)) {
    return_crlist <- return_crlistx
  }

  # record jump-center information if it both was requested and a jump-center is incorporated in the plot
  #if ((jumpinfo) && (repair) && (TRUE %in% crlist$repairinfo["done",])) {
  # record jump-center information if a jump-center is incorporated in the plot
  if ((repair) && (TRUE %in% crlist$repair["done",])) {
    #print("here's repairinfo:")
    #print(crlist$repairinfo)
    a <- (crlist$repair)                      # 'a' is a placeholder to compile info on jump-center repair items
    keep <- as.numeric(which(a[12,] == TRUE))     # identify quadrants with repairs completed ("done" row is TRUE)
    if (xyswap) {
      a <- a[, c(1, 4, 3, 2)]    # q2 and q4 columns swap
      aswap <- matrix(a, ncol = 4)
      rownames(aswap) <- c("leftphi", "lefty", "leftx", "rightphi", "righty", "rightx", "jumpphi",
                           "jumpy", "jumpx", "phi1", "phi2", "done", "gaptype")
      # angles, after reflected through the x = y diagonal, are = (5 * pi) / 2 - (original phi angle)
      needadjust <- ceiling((pi / 2 - as.numeric(aswap[c(1, 4, 7, 10, 11), keep])) / 10)  # is 1 if angle > 2pi, 0 o.w.
      aswap[c(1, 4, 7, 10, 11), keep] <- (5 * pi) / 2 - as.numeric(aswap[c(1, 4, 7, 10, 11), keep]) - 2 * pi * needadjust
      #print(aswap)
      a <- aswap
    }
    dimnames(a) <- list(rownames(a), colnames(a, do.NULL = FALSE, prefix = "quad"))
    a <- a[, keep]
    #print(a)
    qlabel <- paste0("q", keep)
    qlabel <- paste0(qlabel, rep(c("jumpuphill", "jumpshift", "jumpxy",
                                   "jumpL", "jumpR",
                                   "gaptype"), each = length(qlabel)))
    mylist <- vector(mode = "list", length = length(qlabel))
    names(mylist) <- qlabel
    nkeep <- length(keep)
    if (nkeep > 1) {
      for (i in 1:nkeep) {
        mylist[ nkeep * 0 + i] <- jumpuphill[keep[i]]
        mylist[ nkeep * 1 + i] <- jumpshift[keep[i]]
        mylist[[nkeep * 2 + i]] <- t(as.numeric(c(a["jumpx", i], a["jumpy", i])))
        mylist[[nkeep * 3 + i]] <- t(as.numeric(c(a["leftx", i], a["lefty", i])))
        mylist[[nkeep * 4 + i]] <- t(as.numeric(c(a["rightx", i], a["righty", i])))
        mylist[ nkeep * 5 + i] <- a["gaptype", i]
      }
    }
    else if (nkeep == 1) {
      mylist[1]   <- jumpuphill[keep]
      mylist[2]   <- jumpshift[keep]
      mylist[[3]] <- t(as.numeric(c(a["jumpx"], a["jumpy"])))
      mylist[[4]] <- t(as.numeric(c(a["leftx"], a["lefty"])))
      mylist[[5]] <- t(as.numeric(c(a["rightx"], a["righty"])))
      mylist[6]   <- a["gaptype"]
    }
    return_crlistx$repair <- mylist

    # add repairinfo to returned attributes when requested via jumpinfo
    if (jumpinfo) {
      return_crlist$repair <- mylist
    }
    #return_crlist$repairinfo <- crlist$repairinfo
    #return_crlist$jumpuphill <- jumpuphill
    #return_crlist$jumpshift <- jumpshift

    # plot jump-center repair references when requested
    if (showplot && showjump) {
      # establish pL, pR, and pJ as matricies for jump-center reference points (left, right, and jump-center)
      # with columns representing respective x & y coordinates, and rows for each respective repair quadrant
      if (nkeep == 1) {
        pL <- matrix(as.numeric(c(a["leftx"], a["lefty"])), ncol = 2)
        pR <- matrix(as.numeric(c(a["rightx"], a["righty"])), ncol = 2)
        pJ <- matrix(as.numeric(c(a["jumpx"], a["jumpy"])), ncol = 2)
      }
      else if (nkeep > 1) {
        pL <- matrix(as.numeric(c(a["leftx", 1:nkeep], a["lefty", 1:nkeep])), ncol = 2)
        pR <- matrix(as.numeric(c(a["rightx", 1:nkeep], a["righty", 1:nkeep])), ncol = 2)
        pJ <- matrix(as.numeric(c(a["jumpx", 1:nkeep], a["jumpy", 1:nkeep])), ncol = 2)
      }
      cptheta1 <- unlist(return_crlistx[1], use.names = FALSE)    # copy theta1 into vector double values
      cptheta2 <- unlist(return_crlistx[2], use.names = FALSE)    # copy theta1 into vector double values
      for (i in 1:nkeep) {                                        # ID & plot jump-center repair points in each quad
        m1R <- abs(cptheta1 - pR[i, 1]) / diff(range(cptheta1))   # match theta1 values
        m2R <- abs(cptheta2 - pR[i, 2]) / diff(range(cptheta2))   # match theta2 values
        m1L <- abs(cptheta1 - pL[i, 1]) / diff(range(cptheta1))   # match theta1 values
        m2L <- abs(cptheta2 - pL[i, 2]) / diff(range(cptheta2))   # match theta2 values
        fm <- which.min(m1R + m2R) + 1                            # "right" points (nearest right boundary point)
        to <- which.min(m1L + m2L) - 1                            # "to" points (nearest left boundary point)
        points(cptheta1[fm:to], cptheta2[fm:to], pch = 16, cex = 0.4, col = "blue")
        #points(cptheta1[fm:to], cptheta2[fm:to], col = "blue")
      }
      points(pL, pch = 24, bg = "green")    # points on left-side of inaccessible region
      points(pR, pch = 24, bg = "yellow")   # points on right-side of inaccessible region
      points(pJ, pch = 24, bg = "red")      # jump-center
      if (pts) {
        legend("topright", legend = c("CR boundary points", "jump-center (JC)", "JC left boundary", "JC right boundary", "JC repair points"),
               pch = c(1, rep(24, 3), 21), pt.bg = c("black", "red", "green", "yellow", "blue"), pt.cex = c(1, 1, 1, 1, 0.7),
               bty = "n", cex = 0.8)
      }
      else {
        legend("topright", legend = c("jump-center (JC)", "JC left boundary", "JC right boundary", "JC repair points"),
               pch = c(rep(24, 3), 21), pt.bg = c("red", "green", "yellow", "blue"), pt.cex = c(1, 1, 1, 0.7),
               bty = "n", cex = 0.8)
      }
    }
  }   # end of:  if ((repair) && (TRUE %in% crlist$repairinfo["done",]))

  if ((info) && (exact) && (alpha != alpha_target)) {
    return_crlist$alpha_adjusted <- alpha
    return_crlist$alpha_target <- alpha_target
  }

  if ((info) || (jumpinfo)) {
    if (!silent) {
      if ((exact) && (alpha != alpha_target)) {
        newnom <- format(100 * (1 - alpha), digits = 3)
        print(paste0(100 * (1 - alpha_target), "% exact confidence region plot complete; made using a ",
                     newnom, "% nominal coverage value and ", length(phi)," boundary points."))
      }
      else {
        print(paste0(100 * (1 - alpha), "% confidence region plot complete; made using ", length(phi)," boundary points."))
      }
    }
    return(return_crlist)
  }
  else if ((!info) && (!jumpinfo) && (!silent)) {
    if ((exact) && (alpha != alpha_target)) {
      newnom <- format(100 * (1 - alpha), digits = 3)
      return(paste0(100 * (1 - alpha_target), "% exact confidence region plot complete; made using a ",
                    newnom, "% nominal coverage value and ", length(phi)," boundary points."))
    }
    else {
      return(paste0(100 * (1 - alpha), "% confidence region plot complete; made using ", length(phi)," boundary points."))
    }
  }
}
