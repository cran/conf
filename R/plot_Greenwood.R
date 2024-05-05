# Greenwood plot
plot_Greenwood <- function(n, lambdaT, lambdaC, alpha, show, table, value, new_plot, col){
  graycol <- 'gray'
  p  <- seq(0,1,by = 0.001)
  h  <- lambdaT / (lambdaT + lambdaC)
  t  <- -log(p)/ T # dependent on p
  Fx <- 1-p^(1/h)

  outcome <- km.outcomes(n)

  newcol <- rep(-1, nrow(outcome))
  outcome <- cbind(outcome, newcol, newcol, newcol)
  colnames(outcome)[n + 5] <- "LB"
  colnames(outcome)[n + 6] <- "UB"
  colnames(outcome)[n + 7] <- "Contribution"
  z = qnorm(1-alpha/2)

  # Create an empty plot
  if (new_plot) {
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = 'p', ylab = 'c(p)', las = 1, font.axis = 5, cex.axis = 0.8)
  }

  # Add a horizontal line at 1 - alpha with red dashed line style
  abline(h = 1 - alpha, lty = 2, col = 'red')
  #Monte Carlo test
  #abline(h = 0.1503, lty = 2, col = 'blue')
  #abline(v = 0.7, lty = 2, col = 'blue') #test for the theoritical result p = 0.7 c(p) = 0.1503

  for (i in 2:nrow(outcome)){
    if (!is.na(outcome[i,n+2])){

      bigs <- sum(outcome[i,2:n]/(c(n:2)*c((n-1):1)), na.rm = TRUE)
      v <- (outcome[i,n+2]^2)*bigs

      outcome[i,n+5] <- max(c(outcome[i,n+2] - z*sqrt(v),0)) # Lower bound
      outcome[i,n+6] <- min(c(outcome[i,n+2] + z*sqrt(v),1)) # Upper bound

    }
  }

  outcome[is.na(outcome[,n+2]),n+5] <- NA
  outcome[is.na(outcome[,n+2]),n+6] <- NA

  # Find unique values of lower and upper bounds
  lim <- sort(unique(c(outcome[,n+5], outcome[,n+6])))
  if (show == 1){
    # Plot the lines
    for (j in 1:(length(lim) - 1)){
      a <- 1*(lim[j] >= outcome[,n+5] & lim[j+1] <= outcome[,n+6])
      a <- replace(a, is.na(a),0)
      include <- which(a == 1)  # covered CI/ row
      poly <- rep(0, 1001)
      for (k in include){
        l <- outcome[k, 1]
        d <- sum(outcome[k, 2:(l + 1)])
        cont <- choose(n,l) * h^d * (1-h)^(l-d) * Fx^l * (1-Fx)^(n-l)
        poly <- poly + cont
      }
      lines(p, poly, col = 'gray', lwd = 0.3)
    }
  }

  for (i in 2:nrow(outcome)){
    if (!is.na(outcome[i,n+2])){
      bigs <- sum(outcome[i,2:n]/(c(n:2)*c((n-1):1)), na.rm = TRUE)
      v <- (outcome[i,n+2]^2)*bigs
      outcome[i,n+5] <- max(c(outcome[i,n+2] - z*sqrt(v),0)) # Lower bound
      outcome[i,n+6] <- min(c(outcome[i,n+2] + z*sqrt(v),1)) # Upper bound
    }
  }


  outcome[(outcome[,n+2] == 1 ),n+5] <- 1
  outcome[(outcome[,n+2] == 1 ),n+6] <- 1
  outcome[(outcome[,n+2] == 0 ),n+5] <- 0
  outcome[(outcome[,n+2] == 0 ),n+6] <- 0

  lim <- sort(unique(c(outcome[,n+5], outcome[,n+6])))

  for (j in 1:(length(lim) - 1)){
    a <- 1*(lim[j] >= outcome[,n+5] & lim[j+1] <= outcome[,n+6])
    a <- replace(a, is.na(a),0) # replace na in a with 0
    include <- which(a == 1)
    poly <- rep(0, 1001) # initializes with 1001 elements, all set to 0
    for (k in include){
      l <- outcome[k, 1]
      d <- sum(outcome[k, 2:(l + 1)])
      cont <- choose(n,l) * h^d * (1-h)^(l-d) * Fx^l * (1-Fx)^(n-l)
      poly <- poly + cont
    }
    lines(p[sum(lim[j] > p):sum(lim[j+1] > p)], poly[sum(lim[j] > p):sum(lim[j+1] > p)], # extract a part from p and a part from poly
          col = col) # extracts the portion of p/poly that falls between lim[j] and lim[j+1]
  }

  if(table == TRUE){
    p <- value
    h  <- lambdaT / (lambdaT + lambdaC)
    t  <- -log(p)/ lambdaT
    Fx <- 1-p^(1/h)

    for (i in 2:nrow(outcome)){
      if (!is.na(outcome[i,n+2])){

        bigs <- sum(outcome[i,2:n]/(c(n:2)*c((n-1):1)), na.rm = TRUE)
        v <- (outcome[i,n+2]^2)*bigs

        outcome[i,n+5] <- max(c(outcome[i,n+2] - z*sqrt(v),0)) # Lower bound
        outcome[i,n+6] <- min(c(outcome[i,n+2] + z*sqrt(v),1)) # Upper bound
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
