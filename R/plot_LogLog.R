#Log-log
plot_LogLog <- function(n, lambdaT, lambdaC, alpha, show, table, value, new_plot, col){
  graycol = 'gray'
  p  = seq(0,1,by = 0.001)
  h  <- lambdaT / (lambdaT + lambdaC)
  t  <- -log(p)/ lambdaT
  Fx <- 1-exp(-(lambdaT+lambdaC)*t)

  outcome <- km.outcomes(n)

  newcol <- rep(-1, nrow(outcome))
  outcome <- cbind(outcome, newcol, newcol, newcol)
  colnames(outcome)[n + 5] <- "LB"
  colnames(outcome)[n + 6] <- "UB"
  colnames(outcome)[n + 7] <- "Contribution"
  z = qnorm(1-alpha/2)

  if (new_plot) {
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = 'p', ylab = 'c(p)', las = 1, font.axis = 5, cex.axis = 0.8)
  }
  # Add a horizontal line at 1 - alpha with red dashed line style
  abline(h = 1 - alpha, lty = 2, col = 'red')
  # Monte Carlo test
  # loglog

  for (i in 2:nrow(outcome)){
    if (!is.na(outcome[i,n+2])){
      bigs <- sum(outcome[i,2:n]/(c(n:2)*c((n-1):1)), na.rm = TRUE)
      v <- (outcome[i,n+2]^2)*bigs
      outcome[i,n+5] <- max(c(outcome[i,n+2]^exp(-z*sqrt(v)/log(outcome[i,n+2])),0)) # Lower bound
      outcome[i,n+6] <- min(c(outcome[i,n+2]^exp(z*sqrt(v)/log(outcome[i,n+2])),1)) # Upper bound
    }
  }

  outcome[is.na(outcome[,n+2]),n+5] <- NA
  outcome[is.na(outcome[,n+2]),n+6] <- NA

  lim <- sort(unique(c(outcome[,n+5], outcome[,n+6])))

  if (show == 1){
    for (j in 1:(length(lim) - 1)){
      a <- 1*(lim[j] >= outcome[,n+5] & lim[j+1] <= outcome[,n+6])
      a <- replace(a, is.na(a),0)
      include <- which(a == 1)
      poly <- rep(0, 1001)
      for (k in include){
        l <- outcome[k, 1]
        d <- sum(outcome[k, 2:(l + 1)])
        cont <- choose(n,l) * h^d * (1-h)^(l-d) * Fx^l * (1-Fx)^(n-l)
        poly <- poly + cont
      }
      lines(p, poly, col = 'gray75', lwd = 0.3)
    }
  }

  outcome[(outcome[,n+2] == 1 ),n+5] <- 1
  outcome[(outcome[,n+2] == 1 ),n+6] <- 1
  outcome[(outcome[,n+2] == 0 ),n+5] <- 0
  outcome[(outcome[,n+2] == 0 ),n+6] <- 0

  lim <- sort(unique(c(outcome[,n+5], outcome[,n+6])))

  for (j in 1:(length(lim) - 1)){
    a <- 1*(lim[j] >= outcome[,n+5] & lim[j+1] <= outcome[,n+6])
    a <- replace(a, is.na(a),0)
    include <- which(a == 1) # fall in the interval
    poly <- rep(0, 1001)
    for (k in include){
      l <- outcome[k, 1]
      d <- sum(outcome[k, 2:(l + 1)])
      cont <- choose(n,l) * h^d * (1-h)^(l-d) * Fx^l * (1-Fx)^(n-l)
      poly <- poly + cont
    }
    lines(p[sum(lim[j] > p):sum(lim[j+1] > p)], poly[sum(lim[j] > p):sum(lim[j+1] > p)], #DELETE + 0.0041
          col = col)
  }

  if(table == TRUE){
    p  <- value
    h  <- lambdaT / (lambdaT + lambdaC)
    t  <- -log(p)/ lambdaT
    Fx <- 1-p^(1/h)

    for (i in 2:nrow(outcome)){
      if (!is.na(outcome[i,n+2])){
        bigs <- sum(outcome[i,2:n]/(c(n:2)*c((n-1):1)), na.rm = TRUE)
        v <- (outcome[i,n+2]^2)*bigs
        outcome[i,n+5] <- max(c(outcome[i,n+2]^exp(-z*sqrt(v)/log(outcome[i,n+2])),0)) # Lower bound
        outcome[i,n+6] <- min(c(outcome[i,n+2]^exp(z*sqrt(v)/log(outcome[i,n+2])),1)) # Upper bound
      }
    }

    outcome[is.na(outcome[,n+2]),n+5] <- NA
    outcome[is.na(outcome[,n+2]),n+6] <- NA
    outcome[(outcome[,n+2] == 1 & !is.na(outcome[,5])),n+5] <- 1
    outcome[(outcome[,n+2] == 1 & !is.na(outcome[,5])),n+6] <- 1

    for(i in 2:2^(n+1)-1){
      l = outcome[i,1]
      d <- sum(outcome[i, 2:(l + 1)])
      outcome[i , n+7] <- choose(n,l) * h^d * (1-h)^(l-d) * Fx^l * (1-Fx)^(n-l)
    }
    outcome[1, n+7] <- choose(n,0)* (1-Fx)^(n)
    print(outcome)
  }
}
