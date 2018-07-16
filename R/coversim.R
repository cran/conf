#' Confidence Region Coverage
#'
#' @description
#' Creates a confidence region and determines coverage results for a corresponding point of interest.
#' Iterates through a user specified number of trials.
#' Each trial uses a random dataset with user-specified parameters (default) or a user specified dataset
#' matrix (\code{'n'} samples per column, \code{'iter'} columns) and returns the corresponding actual coverage results.
#'
#' @param alpha significance level; resulting plot illustrates a 100(1 - alpha)\% confidence region.
#' @param distn distribution to fit the dataset to; accepted values: \code{'gamma'}, \code{'invgauss'},
#' \code{'lnorm'}, \code{'llogis'}, \code{'norm'}, \code{'weibull'}.
#' @param n trial sample size (producing each confidence region); needed if a dataset is not given.
#' @param iter iterations of individual trials per parameterization; needed if a dataset is not given.
#' @param dataset a \code{'n'} x \code{'iter'} matrix of dataset values.
#' @param point coverage is assessed relative to this point.
#' @param seed random number generator seed.
#' @param kappa distribution parameter (when applicable).
#' @param lambda distribution parameter (when applicable).
#' @param mu distribution parameter (when applicable).
#' @param sigma distribution parameter (when applicable).
#' @param theta distribution parameter (when applicable).
#' @param heuristic numeric value selecting method for plotting: 0 for elliptic-oriented point distribution, and
#' 1 for smoothing boundary search heuristic.
#' @param maxdeg maximum angle tolerance between consecutive plot segments in degrees.
#' @param ellipse_n number of roughly equidistant confidence region points to plot using the
#' elliptic-oriented point distribution (must be a multiple of four because its algorithm
#' exploits symmetry in the quadrants of an ellipse).
#' @param pts displays confidence region boundary points if \code{TRUE}.
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
#' @param info logical argument to return coverage information in a list; includes alpha value(s), n value(s), coverage
#' and error results per iteration, and \code{returnsamp} and/or \code{returnquant} when requested.
#' @param returnsamp logical argument; if \code{TRUE} returns random samples used in a matrix with n rows, iter cols.
#' @param returnquant logical argument; if \code{TRUE} returns random quantiles used in a matrix with n rows, iter cols.
#' @param repair logical argument to repair regions inaccessible using a radial angle from its MLE (multiple root azimuths).
#' @param showplot logical argument specifying if each coverage trial produces a plot.
#' @param delay numeric value of delay (in seconds) between trials so its plot can be seen.
#' @import stats
#' @import graphics
#' @importFrom SDMTools pnt.in.poly
#' @importFrom STAR rllogis pllogis
#' @importFrom utils capture.output
#' @importFrom statmod dinvgauss pinvgauss qinvgauss rinvgauss
#' @export
#' @return if the optional argument \code{info = TRUE} is included then a list of coverage results is returned.
#' @concept confidence region plot
#' @keywords Graphical Methods, Parameter Estimation, Numerical Optimization
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
#'                 kappa     = NULL,
#'                 lambda    = NULL,
#'                 mu        = NULL,
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
#'                 main      = paste0("alpha: ", a, "  n: ", samp),
#'                 xlas      = 1,
#'                 ylas      = 2,
#'                 origin    = FALSE,
#'                 xlim      = NULL,
#'                 ylim      = NULL,
#'                 tol       = .Machine$double.eps ^ 0.5,
#'                 info      = FALSE,
#'                 returnsamp   = FALSE,
#'                 returnquant   = FALSE,
#'                 repair    = TRUE,
#'                 showplot  = FALSE,
#'                 delay     = 0 )
#'
#' @details
#' This package uses the following parameterization for its supported distributions, and illustrates
#' the corresponding confidence regions accordingly:
#'
#' \itemize{
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
#' The confidence region horizontal and vertical axis convention in use by \code{crplot} for each
#' distribution is:
#' \itemize{
#' \item The gamma distribution confidence region plot shows \eqn{\theta} on its horizontal axis,
#' and \eqn{\kappa} on its vertical axis.
#'
#' \item The inverse Gaussian distribution confidence region plot shows \eqn{\mu} on its horizontal
#' axis, and \eqn{\lambda} on its vertical axis.
#'
#' \item The log normal distribution confidence region plot shows \eqn{\mu} on its horizontal axis,
#' and \eqn{\sigma} on its vertical axis.
#'
#' \item The log logistic distribution confidence region plot shows \eqn{\lambda} on its
#' horizontal axis, and \eqn{\kappa} on its vertical axis.
#'
#' \item The normal distribution confidence region plot shows \eqn{\mu} on its horizontal axis, and
#' \eqn{\sigma} on its vertical axis.
#'
#' \item The uniform distribution confidence region plot shows \eqn{a} on its horizontal axis, and
#' \eqn{b} on its vertical axis.
#'
#' \item The Weibull distribution confidence region plot shows \eqn{\kappa} on its horizontal axis,
#' and \eqn{\lambda} on its vertical axis.
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
                     kappa =   NULL,
                     lambda =  NULL,
                     mu =      NULL,
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
                     main = paste0("alpha: ", a, "  n: ", samp),
                     xlas = 1,
                     ylas = 2,
                     origin = FALSE,
                     xlim = NULL,
                     ylim = NULL,
                     tol = .Machine$double.eps ^ 0.5,
                     info = FALSE,
                     returnsamp = FALSE,
                     returnquant = FALSE,
                     repair = TRUE,
                     showplot = FALSE,
                     delay = 0) {


  # parameter error checking ###########################################################

  if (missing(alpha)) stop ("argument 'alpha' is missing, with no default")
  if (missing(distn)) stop ("argument 'distn' is missing, with no default")

  if (is.null(dataset) && (is.null(n) || is.null(iter)))
    stop ("both 'n' and 'iter', or 'dataset' are required to parameterize the simulation")

  if (!is.null(dataset) && (length(dataset) == 1 || !is.numeric(dataset)))
    stop("'dataset' must be a numeric vector with length > 1")

  if (!is.null(dataset) && (length(alpha) != 1))
    stop("when using 'dataset', a scalar value for 'alpha' argument is required")

  if (!is.null(dataset) && (!is.null(c(lambda, kappa, theta, mu, sigma))))
    warning("distribution parameters are ignored when 'dataset' is given")

  if ((!is.null(dataset) && !is.null(n)) && (nrow(dataset) != n))
    stop("when using 'dataset', 'n' must correspond to nrow(dataset) if given")

  if ((!is.null(dataset) && !is.null(iter)) && (ncol(dataset) != iter))
    stop("when using 'dataset', 'iter' must correspond to ncol(dataset) if given")

  if (!is.null(dataset) && is.null(point))
    stop("'point' is required when using 'dataset'")

  if (!is.null(point) && (length(point) !=2 || !is.numeric(point)))
    stop("'point' must be a numeric vector with length 2")

  if (!is.element(distn, c("weibull", "invgauss", "norm", "lnorm", "llogis", "gamma")))
    stop("'distn' invalid; only 'gamma', 'invgauss', 'lnorm', 'llogis', 'norm', 'weibull' are supported")

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

  if (is.null(alpha) || alpha <= 0 || alpha >= 1 || !is.numeric(alpha))
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
    if (distn == "gamma") {
      xy <- c(x$theta, x$kappa)
    }
    else if (distn == "invgauss") {
      xy <- c(x$mu, x$lambda)
    }
    else if (distn == "llogis") {
      xy <- c(x$lambda, x$kappa)
    }
    else if (distn %in% c("lnorm", "norm")) {
      xy <- c(x$mu, x$sigma)
    }
    else if (distn == "weibull") {
      xy <- c(x$kappa, x$lambda)
    }
    pgon <- matrix(xy, nrow = length(x$phi), ncol = 2)
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

  if (!is.null(seed))  set.seed(seed)             # for the option of repeatable results
  if ((!is.null(dataset)) && (is.null(n)))  n <- nrow(dataset)
  if ((!is.null(dataset)) && (is.null(iter)))  iter <- ncol(dataset)

  # initialize counter & variables to store results:
  count <- 0
  allcoverage <- rep(0, length(n) * length(alpha))           # store % coverage per parameterization
  allerrors <- allresults <- matrix(rep(0, length(n) * length(alpha) * iter), nrow = length(n) * length(alpha))

  # iterate through parameterization permutations and record results next
  # note: some functionallity is commented-out to speed calculations when those features are note needed
  for (samp in n) {
    print(paste0("------------------ ", samp, " samples in each of ", iter, " iterations ---------------------"))
    for (a in alpha) {
      cat("\n")
      print(paste0("...now using alpha: ", a))
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
        if (distn == "gamma") {
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
        else if (distn == "norm") {
          samples <- matrix(rnorm(samp * iter, mean = mu, sd = sigma), nrow = samp)
          ru01 <- pnorm(samples, mean = mu, sd = sigma)
        }
        else if (distn == "lnorm") {
          samples <- matrix(rlnorm(samp * iter, meanlog = mu, sdlog = sigma), nrow = samp)
          ru01 <- pnorm(samples, mean = mu, sd = sigma)
        }
        else if (distn == "weibull") {
          samples <- matrix(rweibull(samp * iter, shape = kappa, scale = 1 / lambda), nrow = samp)
          ru01 <- pweibull(samples, shape = kappa, scale = 1 / lambda)
        }
      }

      assessed <- 0           # initialize counter to record number of m.c. sims assessed (does not count errors)
      for (i in 1:iter) {
        invisible(utils::capture.output(
          x <- try(crplot(dataset = samples[, i],
                          alpha = a,
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
                          origin = origin,
                          xlim = xlim,
                          ylim = ylim,
                          tol = tol,
                          repair = repair,
                          showplot = showplot,
                          info = TRUE))
        ))


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
          if (distn == "gamma") {
            p <- c(theta, kappa)
          }
          else if (distn == "invgauss") {
            p <- c(mu, lambda)
          }
          else if (distn %in% c("lnorm", "norm"))  {
            p = c(mu, sigma)
          }
          else if (distn == "llogis") {
            p = c(lambda, kappa)
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
    }                # end alpha for-loop
    cat("\n")
  }                  # end n for-loop

  # plot results
  #par(mfrow = c(1, 1))
  if (distn == "weibull") {
    # use if cycling alpha values but keeping sample size constant:
    if (length(alpha) > 1) {
      plot(1 - alpha, allcoverage, main = paste0("kappa: ", kappa, ", lambda: ", lambda, ", samp: ", n, ", iter: ", iter),
           ylab = "coverage", xlab = "stated coverage (1 - alpha)")
    }
    # use if cycling through sample sizes but keeping alpha constant
    if (length(n) > 1) {
      plot(n, allcoverage, main = paste0("kappa: ", kappa, ", lambda: ", lambda, ", alpha: ", a, ", iter: ", iter),
           ylab = "coverage", xlab = "sample size")
    }
  }

  if (distn == "invgauss") {
    # use if cycling alpha values but keeping sample size constant:
    if ((length(alpha) > 1) && (length(n) == 1)) {
      plot(1 - alpha, allcoverage, main = paste0("mu: ", mu, ", lambda: ", lambda, ", samp: ", n, ", iter: ", iter),
           ylab = "coverage", xlab = "stated coverage (1 - alpha)")
    }
    # use if cycling through sample sizes but keeping alpha constant:
    if ((length(alpha) == 1) && (length(n) > 1)) {
      plot(n, allcoverage, main = paste0("mu: ", mu, ", lambda: ", lambda, ", runs per sample size: ", iter,
                                        " alpha: ", a), ylab = "coverage", xlab = "sample size")
    }
  }

  if (info) {
    # return coverage results; include samples and runif(0, 1) quantiles upon request
    alab <- rep(alpha, length(n))
    nlab <- rep(n, each = length(alpha))
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
