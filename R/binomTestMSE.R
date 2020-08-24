#' @import rootSolve
#' @importFrom utils flush.console
#' @export
#######################################################################################
#
# This R function calculates RMSE-minimizing confidence interval bounds for the
# binomial proportion
#
# Name             : binomTestMSE
# Authors          : Kexin Feng & Heather Sasinowska & Larry Leemis
# Language         : R (part of conf package)
# Latest Revision  : July 2020
#
#######################################################################################

binomTestMSE <- function (n, x, alpha = 0.05, smooth = 1, showRMSE = TRUE, showAll = FALSE) {

  X <- x
  smo <- smooth

  if (!is.numeric(n) || length(n) != 1 || floor(n) != n || n <= 0)
    stop("'n' must be a positive integer")
  if (missing(n))
    stop("argument 'n' is missing, with no default")
  if (!is.numeric(x) || length(x) != 1 || floor(x) != x || x < 0 || x > n)
    stop("'x' must be a positive integer between 0 and n")
  if (missing(x))
    stop("argument 'x' is missing, with no default")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha >= 0.5 || alpha < 0.001)
    stop("'alpha' must be a numeric value in the interval [0.001, 0.5)")
  if (!is.numeric(smooth) || length(smooth) != 1 || (smooth != 0 && smooth != 1))
    stop("'smooth' must be 0 or 1")
  if (n == 1 && ( alpha <= 0 || alpha >= 0.25 ))
    stop("'alpha' must be a numeric value in the interval [0.001, 0.25) for n = 1")
  if (!identical(showRMSE, TRUE) && !identical(showRMSE, FALSE))
    stop("'showRMSE' must be TRUE or FALSE ")
  if (!identical(showAll, TRUE) && !identical(showAll, FALSE))
    stop("'showAll' must be TRUE or FALSE ")

  if (n >= 10 & n<= 12) {
    cat("It takes a few seconds to produce the interval(s).\n")
    cat(paste0(choose(n, floor(n / 2)), " symmetric Dyck paths are enumerated.\n"))
  }

  if (n == 13) {
    cat("It takes several seconds to produce the interval(s).\n")
    cat(paste0(choose(n, floor(n / 2)), " symmetric Dyck paths are enumerated.\n"))
  }

  if (n == 14) {
    cat ("It takes about a minute to produce the interval(s).\n" )
    cat(paste0(choose(n, floor(n / 2)), " symmetric Dyck paths are enumerated.\n"))
  }

  if (n == 15) {
    cat ("It takes several minutes to produce the interval(s).\n")
    cat(paste0(choose(n, floor(n / 2)), " symmetric Dyck paths are enumerated.\n"))
  }

  if (n > 15) {
    cat(paste0("There are ", choose(n, floor(n / 2)), " symmetric Dyck paths.\n"))
    cat("Instead of enumerating all of the symmetric Dyck paths,\nonly the Dyck paths of Wilson-score, Jeffreys, Arcsine,\nand Agresti-Coull are used in the procedure to save considerable time.\n")
    cat("Therefore, the procedure does NOT guarantee a lowest RMSE.\nHowever, it does achieve a lower RMSE than the four existing confidence intervals.\n")

    if (n >= 400 & n<= 4000) {
      cat("It takes a few seconds to produce the interval(s).\n")
    }

    if (n > 4000 & n <= 20000) {
      cat("It takes several seconds to produce the interval(s).\n")
    }

    if (n > 20000 & n <= 30000) {
      cat ("It takes about a minute to produce the interval(s).\n" )
    }

    if (n > 30000) {
      cat ("It takes several minutes to produce the interval(s).\n")
    }
  }

  if (n > 100 && smooth == 0) {
    smooth <- 1
    smo <- 1
    cat ("\nThe smooth parameter is changed to 1,\nsince non-smoothed interval is only allowed for n < 100.\n\n")
  }

  utils::flush.console()

  #library(rootSolve)  # placeholder prior to incorporating binomTestMSE in conf
  #library(conf)       # placeholder prior to incorporating binomTestMSE in conf


  if (n == 1) {
    X <- x
    f1 <- function(p) {2 * (1 - p) + p - 2 * (1 - alpha)}
    p1 <- uniroot(f1, c(0, 1), tol = 0.00001)$root
    p2 <- 1 - p1
    LB <- c(0, p1)
    UB <- c(p2, 1)
    d1 <- c(0, 0, 1)
    d2 <- c(0, 1, 1)
    xx <- rbind(d1, d2)
    p <- rbind(LB, UB)
    pp <- sort(p)
    m <- 0
    for (i in 1:(2 * n + 1)) {
      for (x in xx[1, i]:xx[2, i]) {
        k <- 0:(n - x)
        m <- m + choose(n, x) * sum(choose(n - x, k) * (-1) ^ k * (pp[i + 1] ^ (k + x + 1) - pp[i] ^ (k + x + 1)) / (k + x + 1))
      }
    }
    v <- 0
    for (i in 1:(2 * n + 1)) {
      for (x in xx[1, i]:xx[2, i]) {
        for (y in xx[1, i]:xx[2, i]) {
          k <- 0:(2 * n - x - y)
          v <- v + choose(n, x) * choose(n, y) * sum(choose(2 * n - x - y, k) * (-1) ^ k * (pp[i + 1] ^ (k + x + y + 1) - pp[i] ^ (k + x + y + 1)) / (k + x + y + 1))
        }
      }
    }
    v <- v - m ^ 2
    rmse <- sqrt(v + (m - (1 - alpha)) ^ 2)
    if (showRMSE) {
      cat(paste0("RMSE = ", rmse, "\n"))
    }
    if (!showAll) {
      return(c(LB[X + 1], UB[X + 1]))
    }
    FBL <- rbind(LB, UB)
    rownames(FBL)[1] <- "Lower Bounds"
    rownames(FBL)[2] <- "Upper Bounds"
    rn <- 0:n
    colnames(FBL) <- paste("x = ", rn)
    return(FBL)
  }
  if (n > 1 && n <= 15) {
    alpha <- alpha
    n <- n
    npt <- choose(n, floor(n / 2))
    DyckWord <- c()
    dyck <- function(x, i, n0, n1) {
      if ((n0 + n1 < n) && (n0 > n1)) {
        i <- i + 1
        x[1, i] <- 0
        n0 <- n0 + 1
        dyck(x, i, n0, n1)
        n0 <- n0 - 1
        x[1, i] <- 1
        n1 <- n1 + 1
        dyck(x, i, n0, n1)
        n1 <- n1 - 1
      }
      if ((n0 == n1) && (n0 + n1 < n)) {
        i <- i + 1
        x[1, i] <- 0
        n0 <- n0 + 1
        dyck(x, i, n0, n1)
        n0 <- n0 - 1
      }
      if (n0 + n1 == n) {
        DyckWord <<- rbind(DyckWord, x)
      }
      return()
    }
    x <- matrix(c(0, rep(3, n - 1)), 1, n)
    dyck(x, 1, 1, 0)
    DyckWord2 <- 1 - DyckWord[ , n:1]
    DyckWord <- cbind(DyckWord, DyckWord2)
    DP1 <- apply(DyckWord, 1, cumsum)
    DP2 <- apply(DyckWord[ , (2 * n):1], 1, cumsum)
    DyckPath <- matrix(0, 2 * npt, 2 * n)
    DyckPath[seq(1, (2 * npt) - 1, 2), ] <- t(DP1)
    DyckPath[seq(2, 2 * npt, 2), ] <- t(DP2)
    DyckPath <- cbind(rep(0, 2 * npt), DyckPath)
    k <- 0
    prob <- c(0)
    list <- c()
    rpt <- 0
    point <- 0
    disc <- c()
    for (i in 1:npt) {
      rpt <- 0
      for (j in 1:n) {
        lim <- point
        f1 <- function(p) {pbinom(DyckPath[2 * i, j], size = n, p) - pbinom(DyckPath[2 * i - 1, j] - 1, size = n, p)}
        f2 <- function(p) {pbinom(DyckPath[2 * i, j + 1], size = n, p) - pbinom(DyckPath[2 * i - 1, j + 1] - 1, size = n, p)}
        f <- function(p) {f1(p) + f2(p) - 2 * (1 - alpha)}
        if (rpt > 0) {
          prob <- cbind(prob, c(point))
          rpt <- rpt - 1
          f <- function(p) {0}
        }
        if (f(lim) * f(0.5) < 0) {
          po <- uniroot.all(f, c(lim, 0.5), tol = 0.00001)
          point <- po[1]
          prob <- cbind(prob, c(point))
        }
        if (f(lim) * f(0.5) > 0) {
          ab1 <- abs(f1(lim) - 1 + alpha)
          ab2 <- abs(f2(lim) - 1 + alpha)
          if (ab1 > ab2) {point <- lim}
          if (ab1 < ab2) {
            rpt <- 0
            k <- 2
            point <- 1
            while ((point == 0 || point == 1) && {j + k <= (n + 1)}) {
              rpt <- rpt + 1
              f3 <- function(p) {pbinom(DyckPath[2 * i, j + k], size = n, p) - pbinom(DyckPath[2 * i - 1, j + k] - 1, size = n, p)}
              f4 <- function(p) {abs(f1(p) - 1 + alpha) - abs(f3(p) - 1 + alpha)}
              point1 <- uniroot.all(f4, c(lim, 0.50001), tol = 0.00001)
              point <- 1
              if (length(point1) != 0) {point <- point1[1]}
              k <- k + 1
  	    }
            if ((j + k > (n + 1)) && (point == 1)) { point <- 0.5 }
          }
          prob <- cbind(prob, c(point))
          k <- 0
        }
      }
      point <- 0
      list <- rbind(list, prob)
      prob <- c(0)
    }
    list2 <- t(apply(1 - list, 1, sort))
    finalist <- cbind(list, list2)
    RMSE <- matrix(0, npt, 1)
    m <- matrix(0, nrow(finalist), (ncol(finalist) - 1))
    v <- matrix(0, nrow(finalist), (ncol(finalist) - 1))
    for (i in 1:nrow(finalist)) {
      for (j in 1:(ncol(finalist) - 1)) {
        fx <- function(p) {pbinom(DyckPath[2 * i, j], size = n, p) - pbinom(DyckPath[2 * i - 1, j] - 1, size = n, p)}
        m[i, j] <- integrate(fx, lower = finalist[i, j], upper = finalist[i, j + 1])$value
        fy <- function(p) {(pbinom(DyckPath[2 * i, j], size = n, p) - pbinom(DyckPath[2 * i - 1, j] - 1, size = n, p)) ^ 2}
        v[i, j] <- integrate(fy, lower = finalist[i, j], upper = finalist[i, j + 1])$value
      }
    }
    RMSE[, 1] <- rowSums(v) - 2 * (1 - alpha) * rowSums(m) + (1 - alpha) ^ 2
    a1 <- min(RMSE[ , 1])
    a2 <- which.min(RMSE[ , 1])
    if (smo == 0) {
      LB <- c(0, finalist[a2, (c(0, diff(DyckPath[2 * a2, ])) == 1)])
      LU <- sort(1 - LB)
      if (showRMSE) {
        cat(paste0("RMSE = ", sqrt(a1), "\n"))
      }
      if (!showAll) {
        return(as.vector(c(LB[X + 1], LU[X + 1])))
      }
      FBL <- rbind(LB, LU)
      rownames(FBL)[1] <- "Lower Bounds"
      rownames(FBL)[2] <- "Upper Bounds"
      rn <- 0:n
      colnames(FBL) <- paste("x = ", rn)
      return(FBL)
    }

    for (z in 1:nrow(finalist)) {
      for (s in 3:(ncol(finalist) - 1)) {
        if (finalist[z,s] == finalist[z, s - 1] && finalist[z, s] == finalist[z, s + 1]) {
          disc <- c(disc, z)
        }
        if (finalist[z,s] == finalist[z, s - 1] && DyckPath[2 * z, s - 1] - DyckPath [2 * z,s - 2] == DyckPath[2 * z , s] - DyckPath[2 * z, s - 1]) {
          disc <- c(disc, z)
        }
      }
    }
    finalist[disc, ] <- NA
    RMSE[disc, ] <- 100
    LLBB <- matrix(0, npt, (n + 1))
    for (i in 1:npt) {
      LLBB[i,] <- c(0, finalist[i, (c(0, diff(DyckPath[2 * i, ])) == 1)])
    }
    for (i in 1:npt) {
      if (is.na(LLBB[i, 2]) == FALSE) {
         diffv <- diff(LLBB[i, ])
         h <- 2:n
         si <- min(diffv[h] / diffv[h - 1])
         if (si <= 1) {RMSE[i, 1] <- 100}
      }
    }
    b1 <- min(RMSE[ , 1])
    b2 <- which.min(RMSE[ , 1])
    opt1 <- b2
    opv <- sqrt(b1)
    opbounds1 <- LLBB[b2, ]
    al <- alpha
    A <- matrix(0, 8, n + 1)
    for (x in 0:n) {
      A[1:2, x + 1] <- binomTest(n, x, alpha = al, intervalType = "Wilson-Score")
      A[3:4, x + 1] <- binomTest(n, x, alpha = al, intervalType = "Jeffreys")
      A[5:6, x + 1] <- binomTest(n, x, alpha = al, intervalType = "Arcsine")
      A[7:8, x + 1] <- binomTest(n, x, alpha = al, intervalType = "Agresti-Coull")
    }
    E <- matrix(0, 4, n)
    i <- 1:4
    j <- 2:(n + 1)
    E[i, j - 1] <- A[2 * i - 1, j] - A[2 * i - 1, j - 1]
    F <- colMeans(E)
    F <- c(0, F)
    G <- matrix(0, 2, n)
    i <- 1:(n - 1)
    G[1, i] <- F[i + 1] - (F[i + 1] - F[i]) / 2
    G[2, i] <- F[i + 1] + (F[i + 2] - F[i + 1]) / 2
    G[1, n] <- G[2, n - 1]
    G[2, n] <- 1 - sum(G[2, ])
    LB <- matrix(0, npt, n + 1)
    for (i in 1:npt) {
      lc <- 1
      for (j in 1:(2 * n - 1)) {
        cc <- DyckPath[2 * i, j + 1] - DyckPath[2 * i, j]
        if (cc == 1) {
          lc <- lc + 1
          f1 <- function(p) {pbinom(DyckPath[2 * i, j], size = n, p) - pbinom(DyckPath[2 * i - 1, j] - 1, size = n, p)}
          f2 <- function(p) {pbinom(DyckPath[2 * i, j + 1], size = n, p) - pbinom(DyckPath[2 * i - 1, j + 1] - 1, size = n, p)}
          f <- function(p) {f1(p) + f2(p) - 2 * (1 - al)}
          po <- uniroot.all(f, c(0, 1), tol = 0.00001)
          if (length(po) == 0) {
            p1 <- LB[i, lc - 1] + G[1, lc - 1]
            p2 <- LB[i, lc - 1] + G[2, lc - 1]
            if (abs(f1(p1) - 1 + al) > abs(f2(p1) - 1 + al)) pr <- p1
            else pr <- p2
          }
          else if (length(po) > 0) {
            for (l in 1:length(po)) {
              if ((LB[i, lc - 1] + G[1, lc - 1]) < po[l] && (po[l] < LB[i, lc - 1] + G[2, lc - 1])) {
                pr <- po[l]
                break()
              }
              else if (abs(po[1] - LB[i, lc - 1] - G[1, lc - 1]) < abs(po[1] - LB[i, lc - 1] - G[2, lc - 1])) {
                pr <- LB[i, lc - 1] + G[1, lc - 1]
  	      }
              else {pr <- LB[i,lc - 1] + G[2, lc - 1]}
            }
          }
          LB[i, lc] <- pr
        }
      }
    }

    LB1 <- LB
    LU <- matrix(0, npt, n + 1)
    LU <-  1 - LB
    LB <- cbind(LB, LU)
    NDW <- matrix(0, npt, 2 * n)
    for (p in 1:npt) {
      l <- matrix(0, 2, n + 1)
      l[2, ] <- 10
      u <- matrix(0, 2, n + 1)
      u[2, ] <- 20
      l[1, ] <- LB1[p, ]
      u[1, ] <- 1 - LB1[p, ]
      bounds <- cbind(l, u)
      bounds <- bounds[ , order(bounds[1, ], decreasing = FALSE)]
      d1 <- rep(0, 2 * n)
      d1 <- as.integer(bounds[2, 2:(2 * n + 1)] == 20)
      NDW[p, ] <- d1
    }
    NDP1 <- apply(NDW, 1, cumsum)
    NDP2 <- apply(NDW[ , (2 * n):1], 1, cumsum)
    NDP <- matrix(0, 2 * npt, 2 * n)
    NDP[seq(1, (2 * npt) - 1, 2), ] <- t(NDP1)
    NDP[seq(2, 2 * npt, 2), ] <- t(NDP2)
    NDP <- cbind(rep(0, 2 * npt), NDP)
    RMSE <- matrix(0, npt, 1)
    for (j in 1:npt) {
      xx <- c()
      p <- c()
      xx <- rbind(NDP[2 * j - 1, ], NDP[2 * j, ])
      p <- rbind(LB1[j, ], LU[j, ])
      pp <- sort(p)
      m <- 0
      for (i in 1:(2 * n + 1)) {
        for (x in xx[1, i]:xx[2, i]) {
          k <- 0:(n - x)
          m <- m + choose(n, x) * sum(choose(n - x, k) * (-1) ^ k * (pp[i + 1] ^ (k + x + 1) - pp[i] ^ (k + x + 1)) / (k + x + 1))
        }
      }
      v <- 0
      for (i in 1:(2 * n + 1)) {
        for (x in xx[1, i]:xx[2, i]) {
          for (y in xx[1, i]:xx[2, i]) {
            k <- 0:(2 * n - x - y)
            v <- v + choose(n, x) * choose(n, y) * sum(choose(2 * n - x - y, k) * (-1) ^ k * (pp[i + 1] ^ (k + x + y + 1) - pp[i] ^ (k + x + y + 1)) / (k + x + y + 1))
          }
        }
      }
      v <- v - m ^ 2
      rmse <- sqrt(v + (m - (1 - al)) ^ 2)
      RMSE[j, 1] <- rmse
    }
    b1 <- min(RMSE[ , 1])
    b2 <- which.min(RMSE[ , 1])
    if (opv <= b1) {
      FRMSE <- opv
      FLB <- opbounds1
      dd <- rbind(DyckPath[opt1 * 2 - 1, ], DyckPath[opt1 * 2, ])
    }
    if (opv > b1) {
      FRMSE <- b1
      FLB <- LB1[b2, ]
      dd <- rbind(NDP[b2 * 2 - 1, ], NDP[b2 * 2, ])
    }
    FUB <- sort(1 - FLB)
    if (showRMSE) {
      cat(paste0("RMSE = ", FRMSE, "\n"))
    }
    if (!showAll) {
      return(c(FLB[X + 1], FUB[X + 1]))
    }
    FBL <- rbind(FLB, FUB)
    rownames(FBL)[1] <- "Lower Bounds"
    rownames(FBL)[2] <- "Upper Bounds"
    rn <- 0:n
    colnames(FBL) <- paste("x = ", rn)
    return(FBL)
  }

  if (n > 15 & smo == 1) {
    n <- n
    al <- alpha
    M1 <- matrix(0, 8, 2 * n + 2)
    list <- c()
    for (i in 0:n) {
      ci <- binomTest(n, i, intervalType = "Wilson-Score")
      ci2 <- c(10, 20)
      B <- rbind(ci, ci2)
      list <- cbind(list, B)
    }
    list <- list[ , order(list[1, ], decreasing = FALSE)]
    M1[1:2, ] <- list
    list <- c()
    for (i in 0:n) {
      ci <- binomTest(n, i, intervalType = "Jeffreys")
      ci2 <- c(10, 20)
      B <- rbind(ci, ci2)
      list <- cbind(list, B)
    }
    list <- list[ , order(list[1, ], decreasing = FALSE)]
    M1[3:4, ] <- list
    list <- c()
    for (i in 0:n) {
      ci <- binomTest(n, i, intervalType = "Arcsine")
      ci2 <- c(10, 20)
      B <- rbind(ci, ci2)
      list <- cbind(list, B)
    }
    list <- list[ , order(list[1, ], decreasing = FALSE)]
    M1[5:6, ] <- list
    list <- c()
    for (i in 0:n) {
      ci <- binomTest(n, i, intervalType = "Agresti-Coull")
      ci2 <- c(10, 20)
      B <- rbind(ci, ci2)
      list <- cbind(list, B)
    }
    list <- list[ , order(list[1, ], decreasing = FALSE)]
    M1[7:8, ] <- list
    dd <- matrix(0, 4, (n + 1))
    for (j in 1:4) {
      d1 <- c()
      for (i in 1:(n + 1)) {
        if (M1[2 * j, i] == 10) { d1 <- c(d1, 0) }
        if (M1[2 * j, i] == 20) { d1 <- c(d1, 1) }
      }
      dd[j, ] <- d1
    }
    dd1 <- 1 - dd[ , (n + 1):1]
    dd <- cbind(dd, dd1)
    DyckPath <- matrix(0, 8, (2 * n + 1))
    for (i in 1:4) {
      for (j in 2:(2 * n + 1)) {
        if (dd[i, j] == 0) {
          DyckPath[2 * i - 1, j] <- DyckPath[2 * i - 1, j - 1]
          DyckPath[2 * i, j] <- DyckPath[2 * i, j - 1] + 1
        }
        if (dd[i, j] == 1) {
          DyckPath[2 * i - 1, j] <- DyckPath[2 * i - 1, j - 1] + 1
          DyckPath[2 * i, j] <- DyckPath[2 * i, j - 1]
        }
      }
    }
    A <- matrix(0, 8, n + 1)
    for (x in 0:n) {
      A[1:2, x + 1] <- binomTest(n, x, alpha = al, intervalType = "Wilson-Score")
      A[3:4, x + 1] <- binomTest(n, x, alpha = al, intervalType = "Jeffreys")
      A[5:6, x + 1] <- binomTest(n, x, alpha = al, intervalType = "Arcsine")
      A[7:8, x + 1] <- binomTest(n, x, alpha = al, intervalType = "Agresti-Coull")
    }
    E <- matrix(0, 4, n)
    i <- 1:4
    j <- 2:(n + 1)
    E[i, j - 1] <- A[2 * i - 1, j] - A[2 * i - 1, j - 1]
    F <- colMeans(E)
    F <- c(0, F)
    G <- matrix(0, 2, n)
    i <- 1:(n - 1)
    G[1, i] <- F[i + 1] - (F[i + 1] - F[i]) / 2
    G[2, i] <- F[i + 1] + (F[i + 2] - F[i + 1]) / 2
    G[1, n] <- G[2, n - 1]
    G[2, n] <- 1 - sum(G[2, ])
    LB <- matrix(0, 4, n + 1)
    for (i in 1:4) {
      lc <- 1
      for (j in 1:(2 * n - 1)) {
        cc <- DyckPath[2 * i, j + 1] - DyckPath[2 * i, j]
        if (cc == 1) {
          lc <- lc + 1
          f1 <- function(p) {pbinom(DyckPath[2 * i, j], size = n, p) - pbinom(DyckPath[2 * i - 1, j] - 1, size = n, p)}
          f2 <- function(p) {pbinom(DyckPath[2 * i, j + 1], size = n, p) - pbinom(DyckPath[2 * i - 1, j + 1] - 1, size = n, p)}
          f <- function(p) {f1(p) + f2(p) - 2 * (1 - al)}
          po <- uniroot.all(f, c(0, 1), tol = 0.00001)
          if (length(po) == 0) {
            p1 <- LB[i, lc - 1] + G[1, lc - 1]
            p2 <- LB[i, lc - 1] + G[2, lc - 1]
            if (abs(f1(p1) - 1 + al) > abs(f2(p1) - 1 + al)) pr <- p1
            else pr <- p2
          }
          else if (length(po) > 0) {
            for (l in 1:length(po)) {
              if ((LB[i, lc - 1] + G[1, lc - 1]) < po[l] && (po[l] < LB[i, lc - 1] + G[2, lc - 1])) {
                pr <- po[l]
                break()
              }
              else if (abs(po[1] - LB[i, lc - 1] - G[1, lc - 1]) < abs(po[1] - LB[i, lc - 1] - G[2, lc - 1])) {
                pr <- LB[i, lc - 1] + G[1, lc - 1]
              }
              else {pr <- LB[i, lc - 1] + G[2, lc - 1]}
            }
          }
          LB[i, lc] <- pr
        }
      }
    }
    LB1 <- LB
    LU <- matrix(0, 4, n + 1)
    LU <-  1 - LB
    finalist <- cbind(LB, LU)
    finalist <- t(apply(finalist, 1, sort))
    NDW <- matrix(0, 4, 2 * n)
    for (p in 1:4) {
      l <- matrix(0, 2, n + 1)
      l[2, ] <- 10
      u <- matrix(0, 2, n + 1)
      u[2, ] <- 20
      l[1, ] <- LB1[p, ]
      u[1, ] <- 1 - LB1[p, ]
      bounds <- cbind(l, u)
      bounds <- bounds[ , order(bounds[1, ], decreasing = FALSE)]
      d1 <- rep(0, 2 * n)
      d1 <- as.integer(bounds[2, 2:(2 * n + 1)] == 20)
      NDW[p, ] <- d1
    }
    NDP1 <- apply(NDW, 1, cumsum)
    NDP2 <- apply(NDW[ , (2 * n):1], 1, cumsum)
    NDP <- matrix(0, 8 , 2 * n)
    NDP[seq(1, 7, 2), ] <- t(NDP1)
    NDP[seq(2, 8, 2), ] <- t(NDP2)
    NDP <- cbind(rep(0, 8), NDP)
    RMSE <- matrix(0, 4, 1)
    for (i in 1:4) {
      m <- 0
      v <- 0
      for (j in 1:(ncol(finalist) - 1)) {
        fx <- function(p) {pbinom(NDP[2 * i,j], size = n, p) - pbinom(NDP[2 * i - 1, j] - 1, size = n, p)}
        a1 <- integrate(fx, lower = finalist[i, j], upper = finalist[i, j + 1])$value
        m <- a1 + m
        fy <- function(p) {(pbinom(NDP[2 * i, j], size = n, p) - pbinom(NDP[2 * i - 1, j] - 1, size = n, p)) ^ 2}
        a2 <- integrate(fy, lower = finalist[i, j], upper=finalist[i, j + 1])$value
        v <- a2 + v
      }
      rmse <- sqrt(v - 2 * (1 - al) * m + (1 - al) ^ 2)
      RMSE[i, 1] <- rmse
    }
    b1 <- min(RMSE[, 1])
    b2 <- which.min(RMSE[, 1])
    FL <- LB1[b2, ]
    FU <- sort(LU[b2, ])
    if (showRMSE) {
      cat(paste0("RMSE = ", b1, "\n"))
    }

    if (!showAll) {
      return(c(FL[X + 1], FU[X + 1]))
    }
    FBL <- rbind(FL, FU)
    rownames(FBL)[1] <- "Lower Bounds"
    rownames(FBL)[2] <- "Upper Bounds"
    rn <- 0:n
    colnames(FBL) <- paste("x = ", rn)
    return(FBL)
  }
  if (n > 15 & smo == 0) {
    n <- n
    alpha <- alpha
    npt <- 4
    M1 <- matrix(0, 8, 2 * n + 2)
    list <- c()
    for (i in 0:n) {
      ci <- binomTest(n, i, intervalType = "Wilson-Score")
      ci2 <- c(10, 20)
      B <- rbind(ci, ci2)
      list <- cbind(list, B)
    }
    list <- list[ , order(list[1, ], decreasing = FALSE)]
    M1[1:2, ] <- list
    list <- c()
    for (i in 0:n) {
      ci <- binomTest(n, i, intervalType = "Jeffreys")
      ci2 <- c(10, 20)
      B <- rbind(ci ,ci2)
      list <- cbind(list, B)
    }
    list <- list[ , order(list[1, ], decreasing = FALSE)]
    M1[3:4, ] <- list
    list <- c()
    for (i in 0:n) {
      ci <- binomTest(n, i, intervalType = "Arcsine")
      ci2 <- c(10, 20)
      B <- rbind(ci, ci2)
      list <- cbind(list, B)
    }
    list <- list[ , order(list[1, ], decreasing = FALSE)]
    M1[5:6, ] <- list
    list <- c()
    for (i in 0:n) {
      ci <- binomTest(n, i, intervalType = "Agresti-Coull")
      ci2 <- c(10, 20)
      B <- rbind(ci, ci2)
      list <- cbind(list, B)
    }
    list <- list[ , order(list[1, ], decreasing = FALSE)]
    M1[7:8, ] <- list
    dd <- matrix(0, 4, (n + 1))
    for (j in 1:4) {
      d1 <- c()
      for (i in 1:(n + 1)) {
        if (M1[2 * j, i] == 10) { d1 <- c(d1, 0) }
        if (M1[2 * j, i] == 20) { d1 <- c(d1, 1) }
      }
      dd[j, ] <- d1
    }
    dd1 <- 1 - dd[ , (n + 1):1]
    dd <- cbind(dd, dd1)
    DyckPath <- matrix(0, 8, (2 * n + 1))
    for (i in 1:4) {
      for (j in 2:(2 * n + 1)) {
        if (dd[i, j] == 0) {
          DyckPath[2 * i - 1, j] <- DyckPath[2 * i - 1, j - 1]
          DyckPath[2 * i, j] <- DyckPath[2 * i, j - 1] + 1
        }
        if (dd[i, j] == 1) {
          DyckPath[2 * i - 1, j] <- DyckPath[2 * i - 1, j - 1] + 1
          DyckPath[2 * i, j] <- DyckPath[2 * i, j - 1]
        }
      }
    }
    k <- 0
    prob <- c(0)
    list <- c()
    rep <- 0
    point <- 0
    disc <- c()
    for (i in 1:npt) {
      rep <- 0
      for (j in 1:n) {
        lim <- point
        f1 <- function(p) {pbinom(DyckPath[2 * i, j], size = n, p) - pbinom(DyckPath[2 * i - 1, j] - 1, size = n, p)}
        f2 <- function(p) {pbinom(DyckPath[2 * i, j + 1], size = n, p) - pbinom(DyckPath[2 * i - 1, j + 1] - 1, size = n, p)}
        f <- function(p) {f1(p) + f2(p) - 2 * (1 - alpha)}
        if (rep > 0) {
          prob <- cbind(prob, c(point))
          rep <- rep - 1
          f <- function(p) {0}
        }
        if (f(lim) * f(0.5) < 0) {
          po <- uniroot.all(f, c(lim, 0.5), tol = 0.00001)
          point <- po[1]
          prob <- cbind(prob, c(point))}
        if (f(lim) * f(0.5) > 0) {
          ab1 <- abs(f1(lim) - 1 + alpha)
          ab2 <- abs(f2(lim) - 1 + alpha)
          if (ab1 > ab2) { point <- lim }
          if (ab1 < ab2) {
            rep <- 0
            k <- 2
            point <- 1
            while ({point == 0} | {point == 1} && {j + k <= (n + 1)}) {
              rep <- rep + 1
              f3 <- function(p) {pbinom(DyckPath[2 * i, j + k], size = n, p) - pbinom(DyckPath[2 * i - 1, j + k] - 1, size = n, p)}
              f4 <- function(p) {abs(f1(p) - 1 + alpha) - abs(f3(p) - 1 + alpha)}
              point1 <- uniroot.all(f4, c(lim, 0.50001), tol = 0.00001)
              point <- 1
              if (length(point1) != 0) { point <- point1[1] }
              k = k + 1}
            if ((j + k > (n + 1) ) & (point == 1)) { point <- 0.5 }
          }
          prob <- cbind(prob, c(point))
          k <- 0
        }
      }
      point <- 0
      list <- rbind(list, prob)
      prob <- c(0)
    }
    list2 <- matrix (0, npt, n + 1)
    for (i in 1:npt) { list2[i, ] <- sort(1 - list[i, ]) }
    finalist <- cbind(list, list2)
    RMSE <- matrix(0, 4, 1)
    for (i in 1:4) {
      m <- 0
      v <- 0
      for (j in 1:(ncol(finalist) - 1)) {
        fx <- function(p) {pbinom(DyckPath[2 * i,j], size = n, p) - pbinom(DyckPath[2 * i - 1, j] - 1, size = n, p)}
        a1 <- integrate(fx, lower = finalist[i, j], upper = finalist[i, j + 1])$value
        m <- a1 + m
        fy <- function(p) {(pbinom(DyckPath[2 * i, j], size = n, p) - pbinom(DyckPath[2 * i - 1, j] - 1, size = n, p)) ^ 2}
        a2 <- integrate(fy, lower = finalist[i, j], upper=finalist[i, j + 1])$value
        v <- a2 + v
      }
      rmse <- sqrt(v - 2 * (1 - alpha) * m + (1 - alpha) ^ 2)
      RMSE[i, 1] <- rmse
    }
    b1 <- min(RMSE[ , 1])
    b2 <- which.min(RMSE[ , 1])
    LB <- c(0, finalist[b2, (c(0, diff(DyckPath[2 * b2, ])) == 1)])
    LU <- sort(1 - LB)
    if (showRMSE) {
      cat(paste0("RMSE = ", b1, "\n"))
    }
    if (!showAll) {
      return(as.vector(c(LB[X + 1], LU[X + 1])))
    }
    FBL <- rbind(LB, LU)
    rownames(FBL)[1] <- "Lower Bounds"
    rownames(FBL)[2] <- "Upper Bounds"
    rn <- 0:n
    colnames(FBL) <- paste("x = ", rn)
    return(FBL)
  }

}

