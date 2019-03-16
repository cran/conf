#' Confidence Region Coverage
#'
#' @description
#' Creates a confidence region and determines coverage results for a corresponding point of interest.
#' Iterates through a user specified number of trials.
#' Each trial uses a random dataset with user-specified parameters (default) or a user specified dataset
#' matrix (\code{'n'} samples per column, \code{'iter'} columns) and returns the corresponding actual coverage results.
#' See the CRAN website https://CRAN.R-project.org/package=conf for a link to a \code{coversim} vignette.
#'
#' @param alpha significance level; scalar or vector; resulting plot illustrates a 100(1 - \code{alpha})\% confidence region.
#' @param distn distribution to fit the dataset to; accepted values: \code{'cauchy'}, \code{'gamma'}, \code{'invgauss'},
#' \code{'logis'}, \code{'llogis'}, \code{'lnorm'}, \code{'norm'}, \code{'unif'}, \code{'weibull'}.
#' @param n trial sample size (producing each confidence region); scalar or vector; needed if a dataset is not given.
#' @param iter iterations (or replications) of individual trials per parameterization; needed if a dataset is not given.
#' @param dataset a \code{'n'} x \code{'iter'} matrix of dataset values, or a vector of length \code{'n'} (for a
#' single iteration).
#' @param point coverage is assessed relative to this point.
#' @param seed random number generator seed.
#' @param a distribution parameter (when applicable).
#' @param b distribution parameter (when applicable).
#' @param kappa distribution parameter (when applicable).
#' @param lambda distribution parameter (when applicable).
#' @param mu distribution parameter (when applicable).
#' @param sigma distribution parameter (when applicable).
#' @param s distribution parameter (when applicable).
#' @param theta distribution parameter (when applicable).
#' @param heuristic numeric value selecting method for plotting: 0 for elliptic-oriented point distribution, and
#' 1 for smoothing boundary search heuristic.
#' @param maxdeg maximum angle tolerance between consecutive plot segments in degrees.
#' @param ellipse_n number of roughly equidistant confidence region points to plot using the
#' elliptic-oriented point distribution (must be a multiple of four because its algorithm
#' exploits symmetry in the quadrants of an ellipse).
#' @param pts displays confidence region boundary points if \code{TRUE} (applies to confidence region plots in which \code{showplot = TRUE}).
#' @param mlelab logical argument to include the maximum likelihood estimate coordinate point (default is \code{TRUE},
#' applies to confidence region plots when \code{showplot = TRUE}).
#' @param sf significant figures in axes labels specified using sf = c(x, y), where x and y represent the optional digits argument
#' in the R function \code{\link{round}} as it pertains the horizontal and vertical labels.
#' @param mar specifies margin values for \code{par(mar = c( ))} (see \code{mar} in \code{\link{par}}).
#' @param xlab string specifying the horizontal axis label (applies to confidence region plots when \code{showplot = TRUE}).
#' @param ylab string specifying the vertical axis label (applies to confidence region plots when \code{showplot = TRUE}).
#' @param main string specifying the plot title (applies to confidence region plots when \code{showplot = TRUE}).
#' @param xlas numeric in {0, 1, 2, 3} specifying the style of axis labels (see \code{las} in \code{\link{par}},
#' applies to confidence region plots when \code{showplot = TRUE}).
#' @param ylas numeric in {0, 1, 2, 3} specifying the style of axis labels (see \code{las} in \code{\link{par}},
#' applies to confidence region plots when \code{showplot = TRUE}).
#' @param origin logical argument to include the plot origin (applies to confidence region plots when \code{showplot = TRUE}).
#' @param xlim two element vector containing horizontal axis minimum and maximum values (applies to confidence region plots
#' when \code{showplot = TRUE}).
#' @param ylim two element vector containing vertical axis minimum and maximum values (applies to confidence region plots
#' when \code{showplot = TRUE}).
#' @param tol the \code{\link{uniroot}} parameter specifying its required accuracy.
#' @param info logical argument to return coverage information in a list; includes \code{alpha} value(s), \code{n} value(s), coverage
#' and error results per iteration, and \code{returnsamp} and/or \code{returnquant} when requested.
#' @param returnsamp logical argument; if \code{TRUE} returns random samples used in a matrix with \code{n} rows, \code{iter} cols.
#' @param returnquant logical argument; if \code{TRUE} returns random quantiles used in a matrix with \code{n} rows, \code{iter} cols.
#' @param repair logical argument to repair regions inaccessible using a radial angle from its MLE (multiple root azimuths).
#' @param exact logical argument specifying if alpha value is adjusted to compensate for negative coverage bias in order to achieve
#' (1 - alpha) coverage probability using previously recorded Monte Carlo simulation results; available for limited values of
#' alpha (roughly <= 0.2--0.3), n (typically n = 4, 5, ..., 50) and distributions (distn suffixes: weibull, llogis, norm).
#' @param showplot logical argument specifying if each coverage trial produces a plot.
#' @param delay numeric value of delay (in seconds) between trials so its plot can be seen (applies when \code{showplot = TRUE}).
#' @import stats
#' @import graphics
#' @importFrom SDMTools pnt.in.poly
#' @importFrom STAR rllogis pllogis
#' @importFrom utils capture.output
#' @importFrom statmod dinvgauss pinvgauss qinvgauss rinvgauss
#' @export
#' @return If the optional argument \code{info = TRUE} is included then a list of coverage results is returned.  That list
#' includes \code{alpha} value(s), \code{n} value(s), coverage and error results per iteration.  Additionally, \code{returnsamp = TRUE}
#' and/or \code{returnquant = TRUE} will result in an \code{n} row, \code{iter} column maxtix of sample and/or sample cdf values.
#' @concept confidence region plot
#' @keywords Graphical Methods, Parameter Estimation, Numerical Optimization
#' @references Weld, C., Loh, A., Leemis, L. (in press), "Plotting Likelihood-Ratio Based Confidence Regions for
#' Two-Parameter Univariate Probability Models", The American Statistician.
#' @seealso \code{\link{crplot}}, \code{\link{uniroot}}
#' @author Christopher Weld (\email{ceweld@email.wm.edu})
#' @author Lawrence Leemis (\email{leemis@math.wm.edu})
#' @keywords confidence region, confidence intervals, statistical graphics, data visualization, coverage
#'
#' @usage
#' coversim(alpha, distn,
#'                 n         = NULL,
#'                 iter      = NULL,
#'                 dataset   = NULL,
#'                 point     = NULL,
#'                 seed      = NULL,
#'                 a         = NULL,
#'                 b         = NULL,
#'                 kappa     = NULL,
#'                 lambda    = NULL,
#'                 mu        = NULL,
#'                 s         = NULL,
#'                 sigma     = NULL,
#'                 theta     = NULL,
#'                 heuristic = 1,
#'                 maxdeg    = 5,
#'                 ellipse_n = 4,
#'                 pts       = FALSE,
#'                 mlelab    = TRUE,
#'                 sf        = c(5, 5),
#'                 mar       = c(4, 4.5, 2, 1.5),
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
#'                 returnsamp  = FALSE,
#'                 returnquant = FALSE,
#'                 repair    = TRUE,
#'                 exact     = FALSE,
#'                 showplot  = FALSE,
#'                 delay     = 0 )
#'
#' @details
#' Parameterizations for supported distributions are given following
#' the default axes convention in use by \code{crplot} and \code{coversim}, which are:
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
#' Each respective distribution is defined below.
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
#'
#' @examples
#' ## assess actual coverage at various alpha = {0.5, 0.1} given n = 30 samples,  completing
#' ## 10 trials per parameterization (iter) for a normal(mean = 2, sd = 3) rv
#' coversim(alpha = c(0.5, 0.1), "norm", n = 30, iter = 10, mu = 2, sigma = 3)
#'
#' ## show plots for 5 iterations of 30 samples each from a Weibull(2, 3)
#' coversim(0.5, "weibull", n = 30, iter = 5, lambda = 1.5, kappa = 0.5, showplot = TRUE,
#' origin = TRUE)
#'


#######################################################################################
# This R function plots confidence regions.
#
# Name             : crplot.R
# Authors          : Chris Weld & Larry Leemis
# Language         : R (part of conf package)
# Latest Revision  : May 2018
#######################################################################################
coversim <- function(alpha,
                     distn,
                     n = NULL,
                     iter = NULL,
                     dataset = NULL,
                     point =   NULL,
                     seed =    NULL,
                     a =       NULL,
                     b =       NULL,
                     kappa =   NULL,
                     lambda =  NULL,
                     mu =      NULL,
                     s =       NULL,
                     sigma =   NULL,
                     theta =   NULL,
                     heuristic = 1,
                     maxdeg = 5,
                     ellipse_n = 4,
                     pts = FALSE,
                     mlelab = TRUE,
                     sf = c(5, 5),
                     mar = c(4, 4.5, 2, 1.5),
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
                     returnsamp = FALSE,
                     returnquant = FALSE,
                     repair = TRUE,
                     exact = FALSE,
                     showplot = FALSE,
                     delay = 0) {


  # parameter error checking ###########################################################

  if (missing(alpha)) stop ("argument 'alpha' is missing, with no default")
  if (missing(distn)) stop ("argument 'distn' is missing, with no default")

  if ((!is.null(iter)) && ((length(iter) != 1) || (floor(iter) != iter))) {
    stop("'iter' must be a scalar integer value")
  }

  if (is.null(dataset) && (is.null(n) || is.null(iter)))
    stop ("both 'n' and 'iter', or 'dataset' are required to parameterize the simulation")

  if (!is.null(dataset) && (length(dataset) == 1 || !is.numeric(dataset)))
    stop("'dataset' must be a numeric vector with length > 1, or an 'n' by 'iter' matrix")

  if (!is.null(dataset) && (length(alpha) != 1))
    stop("when using 'dataset', a scalar value for 'alpha' argument is required")

  if (!is.null(dataset) && (!is.null(c(lambda, kappa, theta, mu, sigma))))
    warning("distribution parameters are ignored when 'dataset' is given")

  if (!is.null(dataset) && !is.null(n)) {
    if (is.matrix(dataset) && (nrow(dataset) != n)) {
      stop("when using a matrix 'dataset', 'n' must correspond to nrow(dataset) if given")
    }
    if (is.vector(dataset) && (length(dataset) != n)) {
      stop("when using a 'dataset' vector, 'n' must correspond to length(dataset) if given")
    }
  }

  if (!is.null(dataset) && !is.null(iter)) {
    if  (is.matrix(dataset) && (ncol(dataset) != iter)) {
      stop("when using a 'dataset' matrix, 'iter' must correspond to ncol(dataset) if given")
    }
    if (is.vector(dataset) && (iter != 1)) {
      stop("when using a 'dataset' vector, 'iter' must be 1 if given")
    }
  }

  if (!is.null(dataset) && is.null(point))
    stop("'point' is required when using 'dataset'")

  if (!is.null(point) && (length(point) !=2 || !is.numeric(point)))
    stop("'point' must be a numeric vector with length 2")

  if (!is.element(distn, c("weibull", "invgauss", "norm", "lnorm", "logis", "llogis", "gamma", "unif", "cauchy")))
    stop("'distn' invalid; supported are: 'cauchy', 'gamma', 'invgauss', 'lnorm', 'logis', 'llogis', 'norm', 'unif', 'weibull'")

  if ((distn == "cauchy") && is.null(dataset)) {
    if (is.null(a) || is.null(s)) {
      stop("gamma distribution requires 'a' and 's' arguments (use ?crplot to view corresponding pdf)")
    }
  }

  if ((distn == "gamma") && is.null(dataset)) {
    if (is.null(kappa) || is.null(theta)) {
      stop("gamma distribution requires 'kappa' and 'theta' arguments (use ?crplot to view corresponding pdf)")
    }
  }

  if (distn == "gamma" && !is.null(dataset)) {
    if (min(dataset) <= 0) {
      stop("'dataset' parameter contains infeasible gamma distribution outcome(s) <= 0")
    }
  }

  if ((distn == "invgauss") && is.null(dataset)) {
    if (is.null(mu) || is.null(lambda)) {
      stop("inverse Gaussian distribution requires 'mu' and 'lambda' arguments (use ?crplot to view corresponding pdf)")
    }
  }

  if (distn == "invgauss" && !is.null(dataset)) {
    if (min(dataset) <= 0) {
      stop("'dataset' parameter contains infeasible inverse Gaussian distribution outcome(s) <= 0")
    }
  }

  if ((distn == "logis") && is.null(dataset)) {
    if (is.null(mu) || is.null(sigma)) {
      stop("logistic distribution requires 'mu' and 'sigma' arguments (use ?crplot to view corresponding pdf)")
    }
  }

  if ((distn == "llogis") && is.null(dataset)) {
    if (is.null(kappa) || is.null(lambda)) {
      stop("log logistic distribution requires 'lambda' and 'kappa' arguments (use ?crplot to view corresponding pdf)")
    }
  }

  if (distn == "llogis" && !is.null(dataset)) {
    if (min(dataset) <= 0) {
      stop("'dataset' parameter contains infeasible loglogistic distribution outcome(s) <= 0")
    }
  }

  if ((distn == "lnorm") && is.null(dataset)) {
    if (is.null(mu) || is.null(sigma)) {
      stop("log normal distribution requires 'mu' and 'sigma' arguments (use ?crplot to view corresponding pdf)")
    }
  }

  if ((distn == "norm") && is.null(dataset)) {
    if (is.null(mu) || is.null(sigma)) {
      stop("normal distribution requires 'mu' and 'sigma' arguments (use ?crplot to view corresponding pdf)")
    }
  }

  if ((distn == "unif") && is.null(dataset)) {
    if (is.null(a) || is.null(b)) {
      stop("uniform distribution requires 'a' and 'b' arguments (use ?crplot to view corresponding pdf)")
    }
  }

  if ((distn == "unif") && is.null(dataset)) {
    if (a >= b) {
      stop("uniform distribution requires that a < b (use ?crplot to view corresponding pdf)")
    }
  }

  #if ((distn == "unif") && (sum(cen) != length(dataset))) {
  #  if (max(as.numeric(dataset[which(cen == 0)] %in% max(dataset))) == 1) {
  #    stop("undefined 'unif' confidence region when max(dataset) corresponds to a (cen = 0) censored value")
  #  }
  #}

  if ((distn == "weibull") && is.null(dataset)) {
    if (is.null(kappa) || is.null(lambda)) {
      stop("Weibull distribution requires 'lambda' and 'kappa' arguments (use ?crplot to view corresponding pdf)")
    }
  }

  if (distn == "weibull" && !is.null(dataset)) {
    if (min(dataset) <= 0) {
      stop("'dataset' parameter contains infeasible Weibull distribution outcome(s) <= 0")
    }
  }

  if (is.null(alpha) || all(alpha <= 0) || all(alpha >= 1) || !is.numeric(alpha))
    stop("'alpha' numeric significance level parameter required such that 0 < alpha < 1")

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

  if (length(sf) != 2 || !is.numeric(sf) || floor(sf)[1] != sf[1] || floor(sf)[2] != sf[2] || min(sf) <= 0)
    stop("'sf' must be a numeric vector with positive integer values of length two")

  if (length(mar) != 4 || !is.numeric(mar) || min(mar) < 0)
    stop("'mar' must be a vector of length four with positive numeric entries")

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
    stop("'tol' numeric parameter given is invalid (default .Machine$double.eps^0.3)")

  if (!is.logical(returnsamp) || length(returnsamp) != 1)
      stop("'returnsamp' must be a single logical parameter")

  if (!is.logical(returnquant) || length(returnquant) != 1)
    stop("'returnquant' must be a single logical parameter")

  if (returnsamp && (length(n) != 1 || length(alpha) != 1))
    warning("'returnsamp' only returns samples from the final parameterization of 'n' and 'alpha'")

  if (returnquant && (length(n) != 1 || length(alpha) != 1))
    warning("'returnquant' only returns quantiles from the final parameterization of 'n' and 'alpha'")

  if (!is.logical(repair) || length(repair) != 1)
    stop("'repair' must be a single logical parameter")

  if (!is.logical(exact) || length(exact) != 1)
    stop("'exact' must be a single logical parameter")

  if (exact && !(distn %in% c("weibull", "llogis", "norm")))
    warning("'exact' not available for this distn; proceeding without it")

  if (!is.logical(showplot) || length(showplot) != 1)
    stop("'showplot' must be a single logical parameter")

  if (!is.numeric(delay) || length(delay) != 1 || delay < 0 )
    stop("'delay' must be a non-negative numeric scalar value")

  if (delay && !showplot)
    message("requested delay ignored because showplot = FALSE")

  ######################################################################################

  # sequence of code given below:
  #
  # This file contains two parts:
  #
  # 1. A function named checkpoint.  It accepts two arguments: a point c(x, y), and a
  # variable x containing information from crplot (from x <- crplot(..., info = TRUE)).  It
  # returns 1 when the point is contained within the confidence region, and 0 otherwise.
  #
  # 2. Code to iterate crplot execution.
  # This section of code uses the user inputs and functions (both given above) to iterate
  # through desired parameterizations, model confidence regions, determine coverage,
  # and return results.  Note: some features of #3 are commented out to streamline code
  # for typical coverage analysis
  #

  ######################################################################################

  # 1. checkpoint function.  Accepts two arguments: a point c(x, y), and a
  # variable x containing information from crplot (from x <- crplot(..., info = TRUE)).
  # It returns 1 when the point is contained within the confidence region, and 0 o.w.
  # It also adds the point of interst in red (missed CR) or green (in CR) to the
  # existing plot.  It leverages the pnt.in.poly function from SDMTools.

  checkpoint <- function(p = c(0, 0), x = 0) {
    if (distn == "cauchy") {
      xy <- c(x$a, x$s)
    }
    else if (distn == "gamma") {
      xy <- c(x$theta, x$kappa)
    }
    else if (distn == "invgauss") {
      xy <- c(x$mu, x$lambda)
    }
    else if (distn == "llogis") {
      xy <- c(x$lambda, x$kappa)
    }
    else if (distn %in% c("lnorm", "norm", "logis")) {
      xy <- c(x$mu, x$sigma)
    }
    else if (distn == "unif") {
      xy <- c(x$a, x$b)
    }
    else if (distn == "weibull") {
      xy <- c(x$kappa, x$lambda)
    }
    if (distn == "unif") {
      pgon <- matrix(xy, nrow = 3, ncol = 2)
    }
    else {
      pgon <- matrix(xy, nrow = length(x$phi), ncol = 2)
    }
    p <- matrix(p, nrow = 1, ncol = 2)
    result <- SDMTools::pnt.in.poly(p, pgon)
    if (showplot == TRUE) {
      if (result$pip == 1) {
        points(p, pch = 21, bg = 'green')
      }
      else {
        points(p, pch = 21, bg = 'red')
        lines(pgon, col = 'red')
      }
      Sys.sleep(delay)                 # pause to view plot
    }
    return(result$pip)
  }


  ######################################################################################
  # 2. Code to iterate crplot execution.
  # This section of code uses the user inputs and functions (both given above) to iterate
  # through desired parameterizations, model confidence regions, determine coverage,
  # and return results.

  if (!is.null(seed)) {                    # for the option of repeatable results
    set.seed(seed)
  }

  if (!is.null(dataset)) {                 # ensure datatset in matrix form and ID n, iter
    if (is.vector(dataset)) {
      dataset <- matrix(dataset, ncol = 1)
    }
    if (is.matrix(dataset)) {
      if (is.null(n)) {
        n <- nrow(dataset)
      }
      if (is.null(iter)) {
        iter <- ncol(dataset)
      }
    }
  }

  # initialize counter & variables to store results:
  count <- 0
  addtosummary <- FALSE
  allcoverage <- alpha_adjust_record <- rep(0, length(n) * length(alpha))   # store % coverage per parameterization
  allerrors <- allresults <- matrix(rep(0, length(n) * length(alpha) * iter), nrow = length(n) * length(alpha))

  # iterate through parameterization permutations and record results next
  # note: some functionallity is commented-out to speed calculations when those features are note needed
  for (samp in n) {
    print(paste0("------------------ ", samp, " samples in each of ", iter, " iterations ---------------------"))
    for (this_alpha in alpha) {
      cat("\n")
      print(paste0("...now using alpha: ", this_alpha))
      count <- count + 1
      errors <- 0
      incount <- 0

      # initialize storage for results
      coverage_vector <- rep(0, iter)
      errors_vector <- rep(0, iter)      # will record the index numbers of any itterations with errors

      # identify samples in n rows and iter columns (trial given in one column)
      if (!is.null(dataset)) {
        if (dim(dataset)[2] == 1) {
          samples <- matrix(dataset, nrow = length(dataset))
        }
        else {
          samples <- dataset
        }
      }
      else {
        if (distn == "cauchy") {
         samples <- matrix(rcauchy(samp * iter, location = a, scale = s), nrow = samp)
         ru01 <- pcauchy(samples, location = a, scale = s)
        }
        else if (distn == "gamma") {
          samples <- matrix(rgamma(samp * iter, shape = kappa, scale = theta), nrow = samp)
          ru01 <- pgamma(samples, shape = kappa, scale = theta)
        }
        else if (distn == "invgauss") {
          samples <- matrix(statmod::rinvgauss(samp * iter, mean = mu, shape = lambda), nrow = samp)
          ru01 <- statmod::pinvgauss(samples, mean = mu, shape = lambda)
        }
        else if (distn == "llogis") {
          samples <- matrix(STAR::rllogis(samp * iter, location = log(1 / lambda), scale = 1 / kappa), nrow = samp)
          ru01 <- STAR::pllogis(samples, location = log(1 / lambda), scale = 1 / kappa)
        }
        else if (distn == "logis") {
          samples <- matrix(rlogis(samp * iter, location = mu, scale = sigma), nrow = samp)
          ru01 <- plogis(samples, location = mu, scale = sigma)
        }
        else if (distn == "norm") {
          samples <- matrix(rnorm(samp * iter, mean = mu, sd = sigma), nrow = samp)
          ru01 <- pnorm(samples, mean = mu, sd = sigma)
        }
        else if (distn == "lnorm") {
          samples <- matrix(rlnorm(samp * iter, meanlog = mu, sdlog = sigma), nrow = samp)
          ru01 <- pnorm(samples, mean = mu, sd = sigma)
        }
        else if (distn == "unif") {
          samples <- matrix(runif(samp * iter, min = a, max = b), nrow = samp)
          ru01 <- punif(samples, min = a, max = b)
        }
        else if (distn == "weibull") {
          samples <- matrix(rweibull(samp * iter, shape = kappa, scale = 1 / lambda), nrow = samp)
          ru01 <- pweibull(samples, shape = kappa, scale = 1 / lambda)
        }
      }

      assessed <- 0                 # initialize counter to record number of m.c. sims assessed (does not count errors)
      exact2 <- exact   # set exact2 to user input value (exact2 is subsequently subject to change)
      for (i in 1:iter) {
        invisible(utils::capture.output(
          x <- suppressMessages(try(crplot(dataset = samples[, i],
                          alpha = this_alpha,
                          distn = distn,
                          heuristic = heuristic,
                          maxdeg = maxdeg,
                          ellipse_n = ellipse_n,
                          pts = pts,
                          mlelab = mlelab,
                          sf = sf,
                          mar = mar,
                          xlab = xlab,
                          ylab = ylab,
                          main = main,
                          xlas = xlas,
                          ylas = ylas,
                          origin = origin,
                          xlim = xlim,
                          ylim = ylim,
                          tol = tol,
                          repair = repair,
                          exact = exact2,
                          showplot = showplot,
                          silent = TRUE,
                          info = TRUE),
                   silent = TRUE))
        ))

        # record adjusted alpha value when 'exact = TRUE' (done within first iteration only)
        if ((typeof(x) == "list") && (i == 1) && (exact) && !is.null(x$alpha_adjusted)) {    # adjust alpha value during first iteration only to reduce runtime
          alphadj3dig <- round(x$alpha_adjusted, digits = 3)
          print(paste0("...adjusting alpha to ", alphadj3dig, " to achieve an actual coverage ", 1 - this_alpha))
          this_alpha <- x$alpha_adjusted            # record adjusted alpha value to use next iter
          alpha_adjust_record[count] <- this_alpha  # record to return in summary matrix
          exact2 <- FALSE                     # turn-off exact for subsequent iter b/c it is recorded in this_alpha
          addtosummary <- TRUE
        }
        else if ((typeof(x) == "list") && (i == 1) && (exact) && is.null(x$alpha_adjusted)) {
          print("...cannot adjust alpha for exact coverage ({n, alpha} combination and/or distn not supported); proceeding anyhow")
          alpha_adjust_record[count] <- this_alpha  # record to return in summary matrix
          exact2 <- FALSE                     # turn-off exact for subsequent iter b/c it is recorded in this_alpha
        }
        else if ((typeof(x) != "list") && (i == 1) && (exact)) {
          print("error within the first plot when 'exact = TRUE' adjustments are made; therefore terminating so seed can re-set")
          stop()
        }

        # if error, record as error, otherwise assess actual coverage of "point"
        # (if given) or population parameters o.w.
        if("try-error" %in% class(x)) {
          allerrors[count, i] <- 1
        }
        else if (!is.null(point)) {
          result <- checkpoint(p = point, x = x)
          allresults[count, i] <- result
        }
        else {
          if (distn == "cauchy") {
            p <- c(a, s)
          }
          else if (distn == "gamma") {
            p <- c(theta, kappa)
          }
          else if (distn == "invgauss") {
            p <- c(mu, lambda)
          }
          else if (distn %in% c("lnorm", "norm", "logis"))  {
            p = c(mu, sigma)
          }
          else if (distn == "llogis") {
            p = c(lambda, kappa)
          }
          else if (distn == "unif") {
            p = c(a, b)
          }
          else if (distn == "weibull") {
            p <- c(kappa, lambda)
          }
          result <- checkpoint(p = p, x = x)
          allresults[count, i] <- result
        }
        if ((((i / iter) * 100) %% 1) == 0) {
          pworking <- format(100 * sum(allresults[count, ]) / max(1, (i - sum(allerrors[count, ]))), digits = 2)
          cat(paste0("\r", ((i / iter) * 100), "% complete"))
                     # ; covered thus far: ", pworking, "%"))
          }
      }                                # end iter for-loop

      # record results and print to screen
      cat("\n")
      incount <- sum(allresults[count, ])
      print(paste0(iter, " replications complete."))
      print(paste0(incount, " confidence regions contained the true parameters."))
      print(paste0(sum(allerrors[count, ]), " total errors."))
      #print(paste0(incount, " of the ", iter - sum(allerrors[count, ]), " confidence regions (",
      #             incount / (iter - sum(allerrors[count, ])), ") contained the true parameters."))
      #print(paste0(incount / iter, " contained the true parameters
      #             (assumes errors associate with misses, as typically the case)."))
      allcoverage[count] <- incount / (iter - sum(allerrors[count, ]))
    }                # end this_alpha for-loop
    cat("\n")
  }                  # end n for-loop

  # assemble a summary matrix to print to the screen
  alab <- rep(alpha, length(n))
  nlab <- rep(n, each = length(alpha))
  notexact <- unique(alpha[alpha_adjust_record == alpha])
  if ((exact) && (length(notexact) > 1)) {
    cat("Note: adjusted alpha for exact coverage is not currently available for alpha value(s):", notexact, "\n")
  }
  if ((exact) && (!is.null(alpha_adjust_record)) && (addtosummary)) {
    matrixsummary <- matrix(c(alab, alpha_adjust_record, nlab, allcoverage, rowSums(allerrors)), ncol = 5)
    colnames(matrixsummary) <- c("alpha", "adjusted_alpha", "n", "coverage", "errors")
  }
  else {
    matrixsummary <- matrix(c(alab, nlab, allcoverage, rowSums(allerrors)), ncol = 4)
    colnames(matrixsummary) <- c("alpha", "n", "coverage", "errors")
  }
  print(paste0("coverage simulation summary (", iter, " replications per parameterization):"))
  print(matrixsummary)

  if (main == "") {
    if (length(n) <= length(alpha)) {
      main <- paste0("samp: ", n, " (iter: ", iter, ")")
    }
    else {
      main <- paste0("alpha: ", alpha, " (iter: ", iter, ")")
    }
  }

  # plot results if > 1 datapoint
  if (length(allcoverage) > 1) {
    nplots <- min(length(n), length(alpha))      # necessary number of plots to display all results
    if (length(n) == 1) {
      plot(1 - alpha, allcoverage, main = main,
           ylab = "actual coverage", xlab = "nominal coverage")
      lines(c(min(c(1-alpha, allcoverage)), max(c(1-alpha, allcoverage))), c(min(c(1-alpha, allcoverage)), max(c(1-alpha, allcoverage))), lty = 3)
    }
    else if (length(alpha) == 1) {
      plot(n, allcoverage, main = main, ###paste0("alpha: ", a, " (iter: ", iter, ")"),
           ylab = "actual coverage", xlab = "sample size")
      lines(c(min(n), max(n)), c(1 - alpha, 1 - alpha), lty = 3)
    }
    else if (nplots <= 9) {
      if (nplots <= 3) {
        par(mfrow = c(1, nplots))
      }
      else if (nplots == 4) {
        par(mfrow = c(2, 2))
      }
      else if (nplots <= 6) {
        par(mfrow = c(2, 3))
      }
      else {
        par(mfrow = c(3, 3))
      }
    }
    else {
      par(mfrow = c(3, 3))
      message("note: only first 9 summary plots will display")
    }
    # plot multiple summaries if 1 < nplots <= 9
    if ((nplots > 1) && (length(n) <= length(alpha))) {      # plot nominal vs actual coverage
      marker <- 1       # use to annotate the start location of results for "current" parameterization
      for (i in 1:min(length(n), 9)) {
        plot(1 - alpha, allcoverage[marker:(marker + length(alpha) - 1)], main = main[i], ###paste0("n: ", n[i], " (iter: ", iter, ")"),
             ylab = "actual coverage", xlab = "nominal coverage", xlim = c(0, 1), ylim = c(0, 1))
        lines(c(0, 1), c(0, 1), lty = 3)
        marker <- marker + length(alpha)
      }
      par(mfrow = c(1, 1))
    }
    else if ((nplots > 1) && (length(n) > length(alpha))) {   # plot n vs actual coverage
      marker <- seq(1, (length(alpha) * length(n)), by = length(alpha))     # use to annotate the start location of results for "current" parameterization
      for (i in 1:min(length(alpha), 9)) {
        plot(n, allcoverage[marker + i - 1], main = main[i],   ###paste0("alpha: ", alpha[i], " (iter: ", iter, ")"),
             ylab = "actual coverage", xlab = "sample size")
        lines(c(min(n), max(n)), c(1 - alpha[i], 1 - alpha[i]), lty = 3)
      }
      par(mfrow = c(1, 1))
    }
  }

  if (info) {
    # return coverage results; include samples and runif(0, 1) quantiles upon request
    #alab <- rep(alpha, length(n))
    #nlab <- rep(n, each = length(alpha))
    if (returnsamp && returnquant) {
      return(list("alab" = alab, "nlab" = nlab, "results" = allresults, "errors" = allerrors, "coverage" = allcoverage,
                  "quantiles" = ru01, "samples" = samples))
    }
    else if (returnquant) {
      return(list("alab" = alab, "nlab" = nlab, "results" = allresults, "errors" = allerrors, "coverage" = allcoverage,
                  "quantiles" = ru01))
    }
    else if (returnsamp) {
      return(list("alab" = alab, "nlab" = nlab, "results" = allresults, "errors" = allerrors, "coverage" = allcoverage,
                  "samples" = samples))
    }
    else {
      return(list("alab" = alab, "nlab" = nlab, "results" = allresults, "errors" = allerrors, "coverage" = allcoverage))
    }
  }
  else if (returnsamp && returnquant) {
    return(list("quantiles" = ru01, "samples" = samples))
  }
  else if (returnsamp) {
    return(list("samples" = samples))
  }
  else if (returnquant) {
    return(list("quantiles" = ru01))
  }

}  # end coversim function



