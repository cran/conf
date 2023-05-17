#' @export
######################################################################
# R function to plot the probability mass functions for the support 
# values of the Kaplan-Meier product-limit estimator for a given sample 
# size n with a probability of observing  a failure h at various times 
# of interest expressed as the cumulative probability associated 
# with X =  min(T, C).
######################################################################
km.surv <- function(n, h, lambda = 10, ev = FALSE, line = FALSE, graydots = FALSE,
                    gray.cex = 1, gray.outline = TRUE, xfrac = TRUE)
{
  ######################################################################
  
  # first, some parameter checking...
  
  if (!is.numeric(n) || length(n) != 1 || floor(n) != n  || n <= 0)
    stop("'n' must be a positive integer")
  if (missing(n))
    stop("argument 'n' is missing, with no default")
  
  ######################################################################
  
  # next user warnings based on the value of n
  if (n > 35) stop("modify code for n > 35 because of CPU restrictions")
  if (n > 9 && n <= 27) cat("please wait, this might take a few seconds\n")
  if (n > 27) cat("please wait, this might take a few minutes\n")
  
  ######################################################################
  
  # memory.limit(size = 160000)

  if (lambda < 10) stop("lambda must be greater than 10")
  
  plot(1, type = "n", xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.1),
       yaxt = 'n', cex.axis=0.7, xlab = 'Cumulative Probability of X',
       ylab = 's')
  
  # establish gray scale
  g <- rep(c(0:9,LETTERS[1:6]),each = 16)
  G <- paste(g,c(0:9,LETTERS[1:6]), recycle0 = T, sep = '')
  gray <- rev(paste('#',G,G,G, sep = ''))
  mm <- seq(0,1,length.out = 256)
  
  s <- km.support(n)
  sup <- c(sort(s$num/s$den), 1.05)
  
  expected <- c()

  for (i in seq(0, 1, by = 1/lambda)){
    pmf <- km.pmf(n, h, i, plot = FALSE)[['P']]
    
    # use exponential to decrease difference
    pmf2 <- exp(2*pmf)/sum(exp(2*pmf))
    pmf3 <- (pmf^2)/sum((pmf^2))
    
    # gray color
    if (graydots){
      loc <- ceiling(pmf*256)
      loc[loc == 0] <- 1
      if (gray.outline == TRUE){points(rep(i, length(sup)), sup, col = 'gray90', pch = 1, cex = gray.cex + 0.05)}
      points(rep(i, length(sup)), sup, col = gray[loc], pch = 16, cex = gray.cex)
    }
    
    else{
      # points(rep(i, length(sup)), sup, cex = pmf * 8, pch = 16)
      # points(rep(i, length(sup)), sup, cex = pmf2 * 10, pch = 16)
      points(rep(i, length(sup)), sup, cex = pmf3 * 2, pch = 16)
    }
    
    if (ev) {
      pmfnorm <- c(pmf[1:(length(pmf) - 1)])/(1 - pmf[length(pmf)])
      supnorm <- sort(s$num/s$den)
      exp <- (pmfnorm %*% supnorm)[1, 1]
      expected <- c(expected, exp)
    }
  
  }
  
  if (ev){
    points(seq(0, 1, by = 1/lambda), expected, pch = 16, col = 'red', cex = 0.5)
  }
  # needs to calculate all the expected values first and plot them
  # at the very end in order to prevent overlapping
  
  if (line){
    lines(seq(0, 1, by = 1/lambda), expected, col = 'red')
  }
  
  sn <- s$num
  sd <- s$den
  m <- length(s$num)
  i <- order(s$num / s$den)
  lab <- c("0", paste(sn[i][2:(m-1)],"/", sd[i][2:(m-1)]), "1", "NA")
  if (n == 1){lab = c('0', '1', 'NA')}
  if (xfrac == TRUE){
    axis(2, at = sup, labels = lab, cex.axis=0.7, las = 1)
  }
  else{
    axis(2, at = c(seq(0, 1, by = 0.1), 1.05), 
         labels = c("0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1", "NA"),
         cex.axis=0.7, las = 1)
  }
}
