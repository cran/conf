#' @export
# ////////////////////////////////////////////////////////////////////////////////////////////
# BEGINNING OF CHRISTOPHE POUZAT'S STAR PACKAGE CODE
# Code within this .R file for the functions:
# invgaussMLE, llogisMLE, gammaMLE, and variations of invgauss and llogis
# are attributed to Christophe Pouzat's (now archived) STAR Package
# conf had reverse dependencies with STAR; CRAN is archiving STAR on 2022-05-23
# references:
# https://cran.microsoft.com/snapshot/2014-12-11/web/packages/STAR/index.html
# https://cran.microsoft.com/snapshot/2014-12-11/web/packages/STAR/STAR.pdf

# ############################################################################################
# ############################################################################################
# invgaussMLE code below from Christophe Pouzat (taken from now archived STAR package)
invgaussMLE <- function (yi, ni = numeric(length(yi)) + 1,
                         si = numeric(length(yi)) + 1,
                         parameterization = "sigma2"
) {

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  if (any(yi < 0))
    stop("yi elements must be non-negative")
  yi <- as.numeric(yi)
  if (!identical(length(yi), length(ni)))
    stop("yi and ni should have the same length")
  if (any(ni < 0))
    stop("ni elements must be non-negative")
  if (!identical(class(ni), "integer") && !identical(ni, round(ni)))
    stop("ni should be a vector of positive integers")
  if (!identical(length(yi), length(si)))
    stop("yi and si should have the same length")
  if (any(si < 0))
    stop("si elements must be non-negative")
  if (!identical(si, round(si)))
    stop("si should be a vector of positive integers")
  if (any(si > ni))
    stop("si elements should not be greater than ni elements")
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the log of the mu parameter,\n",
                   "  component 2 is the log of the sigma2 parameter,\n",
                   "using the 'sigma2' parameterization of the inverse Gaussian.\n")
      cat(txt)
    } else {
      mu <- exp(p[1])
      sigma2 <- exp(p[2])
      -(ifelse(s.dot > 0,
               sum(dinvgauss(yi[si > 0], mu, sigma2, log = TRUE) * si[si > 0], 0)
      ) +
        ifelse(c.dot > 0,
               sum(pinvgauss(yi[ci > 0], mu, sigma2, lower.tail = FALSE,log.p = TRUE) * ci[ci > 0]), 0))
    }
  }
  s.dot <- sum(si)
  n.dot <- sum(ni)
  ci <- ni - si
  c.dot <- sum(ci)
  if (s.dot == n.dot) {
    mu.hat <- weighted.mean(yi, ni)
    inv.mean <- weighted.mean(1/yi, ni)
    sigma2.hat <- inv.mean - 1/mu.hat
    estimate <- c(mu.hat, sigma2.hat)
    observedI <- matrix(c(n.dot/(sigma2.hat * mu.hat^3),
                          0, 0, n.dot/sigma2.hat^2/2), nrow = 2, byrow = TRUE)
    se <- sqrt(diag(solve(observedI)))
    l <- -minusLogLik(log(estimate))
  } else {
    if (s.dot >= 10) {
      mu.hat <- weighted.mean(yi, si)
      inv.mean <- weighted.mean(1/yi, si)
      sigma2.hat <- inv.mean - 1/mu.hat
    } else {
      mu.hat <- weighted.mean(yi, ni)
      inv.mean <- weighted.mean(1/yi, ni)
      sigma2.hat <- inv.mean - 1/mu.hat
    }
    mleFit <- optim(par=log(c(mu.hat, sigma2.hat)),
                    fn=minusLogLik,
                    method="BFGS",
                    hessian = TRUE)
    estimate <- exp(mleFit$par)
    newVar <- (1/estimate) %o% (1/estimate)
    observedI <- mleFit$hessian * newVar
    se <- sqrt(diag(solve(observedI)))
    l <- -mleFit$value
  }
  if (parameterization == "sigma2") {
    names(estimate) <- c("mu", "sigma2")
    names(se) <- c("mu", "sigma2")
    # rFct <- function(mu, sigma2) -minusLogLik(log(c(mu, sigma2))) - l
  } else {
    boundary.hat <- (sigma2.hat)^(-0.5)
    mu.hat <- mu.hat/boundary.hat
    if (s.dot == n.dot) {
      estimate <- c(mu.hat, boundary.hat)
      observedI <- observedI * matrix(c(boundary.hat^2,
                                        -2/boundary.hat^2, -2/boundary.hat^2, 4/boundary.hat^6),
                                      nrow = 2, byrow = TRUE)
      se <- se * c(1/boundary.hat, boundary.hat^3/2)
    } else {
      estimate <- c(estimate[1] * sqrt(estimate[2]), 1/sqrt(estimate[2]))
      newVar <- newVar * (c(estimate[2], -2/estimate[2]^3) %o%
                            c(estimate[2], -2/estimate[2]^3))
      observedI <- mleFit$hessian * newVar
      se <- sqrt(diag(solve(observedI)))
    }
    names(estimate) <- c("mu", "boundary")
    names(se) <- c("mu", "boundary")
    # rFct <- function(mu, boundary) -minusLogLik(log(c(mu * boundary, 1/boundary^2))) - l
  }

  # (added rFct function below to remedy multiple local function definitions for 'rFct')
  rFct <- function(mu, myparameter) {
    if (parameterization == "sigma2"){
      myresult <- -minusLogLik(log(c(mu, myparameter))) - l   # myparameter is sigma2
    } else {
      myresult <- -minusLogLik(log(c(mu * myparameter, 1/myparameter^2))) - l  # myparameter is boundary
    }
    return(myresult)
  }
  # (end of addition)

  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
  )
  class(result) <- "durationFit"
  return(result)

}

#' @export
# ############################################################################################
# dinvgauss code below from Christophe Pouzat (taken from now archived STAR package)
dinvgauss <- function(x,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL,
                      log = FALSE
){

  if( any(x <= 0) ) stop("y must contain positive values")
  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }

  tmp <- -(x-mu)^2/(2*x*sigma2*mu^2)-(log(2*pi*sigma2)+3*log(x))/2
  if(!log) tmp <- exp(tmp)
  tmp

}

#' @export
# ############################################################################################
# pinvgauss code below from Christophe Pouzat (taken from now archived STAR package)
pinvgauss <- function(q,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL,
                      lower.tail = TRUE,
                      log.p = FALSE
){

  if( any(q <= 0) ) stop("q must contain positive values")
  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }
  t <- q/mu
  v <- sqrt(q*sigma2)

  ## Use Eq. 4 of Whitemore GA and Seshadri V (1987)
  ## The American Statistician 41:280-281
  if (lower.tail & !log.p)
    return(pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v))
  if (!lower.tail & !log.p)
    return(1 - (pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v)))
  if (lower.tail & log.p)
    return(log(pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v)))
  if (!lower.tail & log.p)
    return(log(1 - (pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v))))

}

# #' @export
# # ############################################################################################
# # hinvgauss code below from Christophe Pouzat (taken from now archived STAR package)
# hinvgauss <- function(x,
#                       mu = 1,
#                       sigma2 = 1,
#                       boundary = NULL,
#                       log = FALSE) {
#
#   if( any(mu <= 0) ) stop("mu must be positive")
#   if (all(!is.null(sigma2)))
#     if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
#   if (all(!is.null(boundary)))
#     if( any(boundary <= 0) ) stop("boundary must be positive")
#   if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
#     stop("One of sigma2 or boundary must be specified, not both")
#
#   if (all(!is.null(boundary))) {
#     ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
#     ## We convert it in term of mu and sigma2
#     sigma2 <- (1/boundary)^2
#     mu <- boundary*mu
#   }
#
#   t <- x/mu
#   v <- sqrt(x*sigma2)
#   cutOff <- 1e-12
#   bigEnough <- pinvgauss(x, mu, sigma2, lower.tail = FALSE) > cutOff
#   logIntensity <- -( (t[bigEnough]-1)^2/(x[bigEnough]*sigma2) + log(2*sigma2*pi*x[bigEnough]^3) )/2 -
#     log(1-pnorm((t[bigEnough]-1)/v[bigEnough])-exp(2/(mu*sigma2))*pnorm(-(t[bigEnough]+1)/v[bigEnough])
#     )
#   if (sum(bigEnough) < length(x)) {
#     cutOffQ <- qinvgauss(1 - cutOff, mu, sigma2)
#     cutOffValue <- dinvgauss(cutOffQ, mu, sigma2, log = TRUE) - log(cutOff)
#     logIntensity <- c(logIntensity,rep(cutOffValue,length(x) - sum(bigEnough)))
#   }
#   if (!log) return(exp(logIntensity))
#   else return(logIntensity)
#
# }

#' @export
# ############################################################################################
# qinvgauss code below from Christophe Pouzat (taken from now archived STAR package)
qinvgauss <- function(p,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL
){

  if( any(p < 0 | p > 1) ) stop("p must lie between 0 and 1")
  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }

  len <- max(length(p),length(mu),length(sigma2))

  if(length(p) != len) {
    if(length(p) == 1) p <- rep(p,len)
    else stop("length of p incorrect")
  }
  if(length(mu) != len) {
    if(length(mu) == 1) mu <- rep(mu,len)
    else stop("length of m incorrect")
  }
  if(length(sigma2) != len) {
    if(length(sigma2) == 1) sigma2 <- rep(sigma2,len)
    else stop("length of sigma2 incorrect")
  }

  ## Use Whitemore and Yalovky (1978, Technometrics, 20:207-208)
  ## approximation to get starting value for the numerical
  ## inversion of the cumulative distribution function.
  theta <- 1/mu/sigma2
  approx <- mu * exp(qnorm(p)*sqrt(1/theta)-0.5/theta)
  sapply(1:len, function(idx) {
    if (identical(p[idx],0)) return(0)
    if (identical(p[idx],1)) return(Inf)
    interval <- approx[idx]*c(0.95,1.05)
    h <- function(q) pinvgauss(q, mu[idx], sigma2[idx]) - p[idx]
    while (h(interval[1])*h(interval[2]) > 0)
      interval <- interval*c(0.9,1.1)
    uniroot(h,interval)$root
  }
  )

}

#' @export
# ############################################################################################
# rinvgauss code below from Christophe Pouzat (taken from now archived STAR package)
rinvgauss <- function(n = 1,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL
){

  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }

  ## Use method of Michael JR, Schucany WR and Haas RW (1976)
  ## The American Statistician, 30:88-90
  v0 <- rchisq(n,1)
  x1 <- mu + 0.5*mu^2*sigma2*v0 - 0.5*mu*sigma2*sqrt(4*mu*v0/sigma2+mu^2*v0^2)
  ifelse(rbinom(length(x1),1,mu/(mu+x1)) == 1,x1,mu^2/x1)
}

coef.durationFit <- function(object,...) object$estimate
logLik.durationFit <- function(object,...) object$logLik

#' @export
# ############################################################################################
# ############################################################################################
# gammaMLE code below from Christophe Pouzat (taken from now archived STAR package)
gammaMLE <- function(yi,
                     ni = numeric(length(yi))+1,
                     si = numeric(length(yi))+1,
                     scale = TRUE
) {
  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi < 0)) stop("yi elements must be non-negative")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")

  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the log of the shape parameter,\n",
                   "  component 2 is the log of the scale parameter.\n"
      )
      cat(txt)
    } else {
      shape <- exp(p[1])
      scale <- exp(p[2])
      -(ifelse(s.dot>0,sum(dgamma(yi[si>0],shape=shape,scale=scale,log=TRUE)*si[si>0],0))+
          ifelse(c.dot>0,sum(pgamma(yi[ci>0],shape=shape,
                                    scale=scale,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
      )
    }
  }

  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)

  ## Define function logProfileLik returning the log of the
  ## profiled likelihood for the shape parameter
  logProfileLik <- function(shape) {
    if (shape <= 0) return(-Inf)
    term1 <- s.dot*shape*(log(shape) - log(scale.hat) -1)
    term2 <- (shape - 1) * sum(ni * log(yi))
    term3 <- - s.dot * lgamma(shape)
    return(term1 + term2 + term3)
  }

  if (s.dot == n.dot) {
    ## no censored event
    ## Get the empirical mean which is the MLE of parameter scale
    scale.hat <- weighted.mean(yi, ni)
    s2 <- weighted.mean(yi^2, ni) - scale.hat^2
    shape.hat <- scale.hat^2 / s2
    shape.lower <- floor(shape.hat)
    shape.upper <- ceiling(shape.hat)
    shape.hat <- optimize(logProfileLik,
                          lower = shape.lower,
                          upper = shape.upper,
                          maximum = TRUE)$maximum
    ## mleFit <- nlm(minusLogLik,log(c(shape.hat,scale.hat)),hessian=TRUE)
    mleFit <- optim(par=log(c(shape.hat,scale.hat)),
                    fn=minusLogLik,
                    method="BFGS",
                    hessian=TRUE)
  } else {
    ## some censored events
    ## if more than 10 events are uncesored get inital guess from them
    ## otherwise use all events
    if (s.dot >= 10) {
      scale.hat <- weighted.mean(yi, si)
      s2 <- weighted.mean(yi^2, si) - scale.hat^2
      shape.hat <- scale.hat^2 / s2
      shape.lower <- floor(shape.hat)
      shape.upper <- ceiling(shape.hat)
      shape.hat <- optimize(logProfileLik,
                            lower = shape.lower,
                            upper = shape.upper,
                            maximum = TRUE)$maximum
    } else {
      scale.hat <- weighted.mean(yi, ni)
      s2 <- weighted.mean(yi^2, ni) - scale.hat^2
      shape.hat <- scale.hat^2 / s2
    }
    ## mleFit <- nlm(minusLogLik,log(c(shape.hat,scale.hat)),hessian=TRUE)
    mleFit <- optim(par=log(c(shape.hat,scale.hat)),
                    fn=minusLogLik,
                    method="BFGS",
                    hessian=TRUE)
  } ## End of conditional on s.dot == n.dot
  ## estimate <- exp(mleFit$estimate)
  estimate <- exp(mleFit$par)
  newVar <- (1/estimate) %o% (1/estimate)
  observedI <- mleFit$hessian * newVar
  se <- sqrt(diag(solve(observedI)))
  ## l <- -mleFit$minimum
  l <- -mleFit$value


  if (scale) {
    names(estimate) <- c("shape","scale")
    names(se) <- c("shape","scale")
    # rFct <- function(shape,scale) -minusLogLik(log(c(shape,scale))) - l
  } else {
    estimate <- c(estimate[1],1/estimate[2])
    names(estimate) <- c("shape","rate")
    se <- se * c(1,estimate[2]^2)
    names(se) <- c("shape","rate")
    # rFct <- function(shape,rate) -minusLogLik(log(c(shape,1/rate))) - l
  }

  # (added rFct function below to remedy multiple local function definitions for 'rFct')
  rFct <- function(shape, myparameter) {
    if (scale){
      myresult <- -minusLogLik(log(c(shape,myparameter))) - l     # myparameter is scale
    } else {
      myresult <- -minusLogLik(log(c(shape,1/myparameter))) - l   # myparameter is rate
    }
    return(myresult)
  }
  # (end of addition)

  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
  )
  class(result) <- "durationFit"
  return(result)

}

#' @export
# ############################################################################################
# ############################################################################################
# llogisMLE code below from Christophe Pouzat (taken from now archived STAR package)
llogisMLE <- function(yi,
                      ni = numeric(length(yi))+1,
                      si = numeric(length(yi))+1
) {

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi <= 0)) stop("yi elements must be positive or null")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")

  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the location parameter,\n",
                   "  component 2 is the log of the scale parameter.\n"
      )
      cat(txt)
    } else {
      location <- p[1]
      scale <- exp(p[2])
      -(ifelse(s.dot>0,sum(dllogis(yi[si>0],location,scale,log=TRUE)*si[si>0],0))+
          ifelse(c.dot>0,sum(pllogis(yi[ci>0],location,scale,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
      )
    }
  }

  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)

  lyi <- log(yi)
  if (s.dot == n.dot) {
    ## no censored event
    ## Get the empirical mean which is the MLE of parameter mu
    location.moment <- weighted.mean(lyi, si)
    scale.moment <- sqrt(3*(weighted.mean(lyi^2, si) - location.moment^2))/pi
  } else {
    ## some censored events
    ## if more than 10 events are uncesored get inital guess from them
    ## otherwise use all events
    if (s.dot >= 10) {
      location.moment <- weighted.mean(lyi, si)
      scale.moment <- sqrt(3*(weighted.mean(lyi^2, si) - location.moment^2))/pi
    } else {
      location.moment <- weighted.mean(lyi, ni)
      scale.moment <- sqrt(3*(weighted.mean(lyi^2, ni) - location.moment^2))/pi
    }
  }
  ## mleFit <- nlm(minusLogLik,c(location.moment,log(scale.moment)),hessian=TRUE)
  mleFit <- optim(fn=minusLogLik,
                  par=c(location.moment,log(scale.moment)),
                  method="BFGS",
                  hessian=TRUE)
  ## estimate <- c(mleFit$estimate[1],exp(mleFit$estimate[2]))
  estimate <- c(mleFit$par[1],exp(mleFit$par[2]))
  newVar <- c(1,1/estimate[2]) %o% c(1,1/estimate[2])
  se <- sqrt(diag(solve(mleFit$hessian * newVar)))
  ## l <- -mleFit$minimum
  l <- -mleFit$value
  names(estimate) <- c("location","scale")
  names(se) <- c("location","scale")
  rFct <- function(location,scale) -minusLogLik(c(location,log(scale))) - l

  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
  )
  class(result) <- "durationFit"
  return(result)

}

#' @export
# ############################################################################################
# dllogis code below from Christophe Pouzat (taken from now archived STAR package)
dllogis <- function(x,
                    location = 0,
                    scale = 1,
                    log = FALSE) {

  ## Check that x elements are strictly positive
  if (any(x <= 0)) stop("x elements must be strictly positive")

  if (!log) dlogis(log(x),location,scale)/x
  else dlogis(log(x),location,scale,log=TRUE) - log(x)

}

#' @export
# ############################################################################################
# pllogis code below from Christophe Pouzat (taken from now archived STAR package)
pllogis <- function(q,
                    location = 0,
                    scale = 1,
                    lower.tail = TRUE,
                    log.p = FALSE) {

  plogis(log(q),location,scale,lower.tail,log.p)

}

#' @export
# ############################################################################################
# qllogis code below from Christophe Pouzat (taken from now archived STAR package)
qllogis <- function(p,
                    location = 0,
                    scale = 1,
                    lower.tail = TRUE,
                    log.p = FALSE) {

  exp(qlogis(p,location,scale,lower.tail,log.p))

}

#' @export
# ############################################################################################
# rllogis code below from Christophe Pouzat (taken from now archived STAR package)
rllogis <- function(n,
                    location = 0,
                    scale = 1) {

  exp(rlogis(n,location,scale))

}

# #' @export
# # ############################################################################################
# # hllogis code below from Christophe Pouzat (taken from now archived STAR package)
# hllogis <- function(x,
#                     location = 0,
#                     scale = 1,
#                     log = FALSE) {
#
#   ## Check that x elements are strictly positive
#   if (any(x <= 0)) stop("x elements must be strictly positive")
#
#   if (!log) pllogis(x,location,scale)/(scale*x)
#   else pllogis(x,location,scale,log.p=TRUE) - log(scale) - log(x)
#
# }

# END OF CHRISTOPHE POUZAT'S STAR PACKAGE CODE CITATION
# ////////////////////////////////////////////////////////////////////////////////////////////
