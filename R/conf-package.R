#' conf: Visualization and Analysis of Statistical Measures of Confidence
#'
#' @description
#' Enables:
#' \enumerate{
#' \item confidence region plots in two-dimensions corresponding to a user given dataset,
#' level of significance, and parametric probability distribution (supported distribution suffixes:
#' cauchy, gamma, invgauss, lnorm, llogis, logis, norm, unif, weibull),
#' \item coverage simulations (if a point of interest is within or outside of a confidence region
#' boundary) for either random samples drawn from a user-specified parametric distribution or for a
#' user-specified dataset and point of interest, and
#' \item calculating confidence intervals and the associated actual coverage for binomial proportions.
#' }
#'
#' \bold{Request from authors}:  Please properly cite any use of this package and/or its algorithms,
#' which are detailed in the corresponding publication by Weld (2018)
#' <doi:10.1080/00031305.2018.1564696>.  Additionally, we welcome and appreciate your feedback and
#' insights as to how this resource is being leveraged to improve whatever it is you do.  Please
#' include your name and academic and/or business affiliation in your correspondance.
#'
#' @details
#' This package includes the functions:
#' \itemize{
#' \item confidence region plots: \code{\link{crplot}},
#' \item confidence region coverage analysis: \code{\link{coversim}},
#' \item confidence intervals for binomial proportions: \code{\link{binomTest}},
#' \item actual coverage calculation for binomial proportions: \code{\link{binomTestCoverage}},
#' \item actual coverage plots for binomial proportions: \code{\link{binomTestCoveragePlot}}, and
#' \item ensemble confidence intervals for binomial proportions: \code{\link{binomTestEnsemble}}.
#' }
#'
#' @section Vignettes:
#' The CRAN website https://CRAN.R-project.org/package=conf contains links for vignettes on the
#' \code{\link{crplot}} and \code{\link{coversim}} functions.
#'
#' @section Acknowledgments:
#' The lead author thanks The Omar Bradley Fellowship for Research in Mathematics for funding that partially
#' supported this work.
#'
#' @author
#' Christopher Weld, Hayeon Park, Larry Leemis
#'
#' Maintainer: Christopher Weld <ceweld@email.wm.edu>
#'
#' @docType package
#' @name conf

NULL
