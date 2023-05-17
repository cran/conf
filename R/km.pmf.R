#' @export
#######################################################################
# R function to calculate the probability mass function of the support
# values of the Kaplan-Meier product-limit estimator, with an option to
# plot the graph
######################################################################

km.pmf <- function(n, h, perc, plot = TRUE, sep = TRUE, xfrac = TRUE,
                   cex.lollipop = 1)
{
######################################################################
  
  # first, some parameter checking...
  
  if (!is.numeric(n) || length(n) != 1 || floor(n) != n  || n <= 0)
    stop("'n' must be a positive integer")
  if (missing(n))
    stop("argument 'n' is missing, with no default")
  if (!is.numeric(h) || length(h) != 1 || h > 1  || h < 0)
    stop("'h' must be a numeric on the interval [0, 1]")
  if (missing(h))
    stop("argument 'h' is missing, with no default")
  if (!is.numeric(perc) || length(perc) != 1 || perc > 1  || perc < 0)
    stop("'perc' must be a numeric on the interval [0, 1]")
  if (missing(perc))
    stop("argument 'perc' is missing, with no default")
  if (!is.logical(plot))
    stop("'plot' must be logical type")
  if (!is.logical(sep))
    stop("'sep' must be logical type")
  if (!is.logical(xfrac))
    stop("'xfrac' must be logical type")
  if (!is.numeric(cex.lollipop) || length(cex.lollipop) != 1)
    stop("'cex.lollipop' must be a numeric")
  
######################################################################
  
  # next user warnings based on the value of n
  if (n >= 24) stop("modify code for n >= 24 because of CPU restrictions")

######################################################################
  
  outcome <- km.outcomes(n)
  newcol <- rep(-1, nrow(outcome))
  outcome <- cbind(outcome, newcol)
  colnames(outcome)[n + 5] <- "prob"
  
  outcome[1, n + 5] <- (1 - perc) ^ n
  for (i in 2:nrow(outcome)) {
    l <- outcome[i, 1]
    d <- sum(outcome[i, 2:(l + 1)])
    prob <- choose(n, l) * (h ^ d) * ((1 - h) ^ (l - d)) * (perc ^ l) * ((1 - perc) ^ (n - l))
    outcome[i, n + 5] = prob
  }
  
  df <- data.frame(outcome[ , (n + 2)], outcome[ , (n + 5)])
  colnames(df) <- c("S", "P")
  totalp <- aggregate(P~S, df, sum)
  
  totalp[(nrow(totalp) + 1), 2] <- 1 - sum(totalp['P'])
  totalp[nrow(totalp), 1] <- NA
  
  if (plot == TRUE){
    s <- km.support(n)
    i <- order(s$num / s$den)
    m <- length(s$num)
    sn <- s$num
    sd <- s$den
    lab <- c("0", paste(sn[i][2:(m - 1)], "/", sd[i][2:(m - 1)]), "1", "NA")

    sup <- c(totalp[1:(nrow(totalp) - 1), 1], 1.1)
    pmf <- totalp[1:nrow(totalp), 2]
    plot(sup, pmf, las = 1, xlab = "", ylab = "", xaxt = 'n', 
    ylim = c(0, max(pmf) + 0.01), type = "h")
    points(sup, pmf, cex = cex.lollipop, pch = 20)
    title(xlab = "s", line = 2)
    title(ylab = expression(P~"("~hat(S)~"("~t[0]~")" == s~")"), line = 2.7)
    if (xfrac == TRUE) {
      axis(1, at = c(sup[1:(length(sup) - 1)], 1.1), labels = lab)
    }
    else{
      axis(1, at = seq(0, 1.1, by = 0.1), 
           labels = c("0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1", "NA"))
    }
    if (sep == TRUE) {
      cat('Additional delay due to printing hash marks.\n')
      df2 <- df[order(df$S), ]
      df2[(nrow(df2) + 1 - 2 ^ (n - 1)):nrow(df2), 1] <- 1.1
      cs <- ave(df2$P,df2$S, FUN = cumsum)
      points(df2[, 1], cs, pch = "-")
    }
  }
  
  if (h == 1){
    totalp[nrow(totalp), 2] <- 0
  }
  
  return(totalp)
  
}
