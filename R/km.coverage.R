#' @export

km.coverage <- function(n, lambdaT, lambdaC, alpha = 0.1, interval = c("Greenwood"), show=TRUE, table = FALSE, value = 0) {
  if (!is.numeric(n) || n != as.integer(n) || n < 1 || n > 15) {
    stop("Error: n must be an integer value between 1 and 15.")
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("Error: alpha must be between 0 and 1.")
  }
  if (!is.numeric(value) || value < 0 || value > 1) {
    stop("Error: selected p value must be between 0 and 1.")
  }
  # Check if other numeric inputs are actually numeric
  if (!is.numeric(lambdaT) || !is.numeric(lambdaC) || !is.numeric(alpha)) {
    stop("Error: lambda T, lambda C, and alpha must be numeric values.")
  }
  # Check if 'show' is a logical value
  if (!is.logical(show)) {
    stop("Error: show must be a logical value (TRUE or FALSE).")
  }
  # Check if the lambda values are valid
  if (lambdaC<0|| lambdaT<0) {
    stop("Error: T, C, must be positive values.")
  }
  # Check if the interval type is valid
  interval_options <- c("Greenwood", "Exp-Greenwood", "Log-Log", "Log-log", "Arcsine", "Peto", "all")
  if (!all(interval %in% interval_options)) {
    stop("Error: Interval type not recognized.")
  }
  if (length(interval) > 1 & ("all" %in% interval)) {
    stop("Error: Interval can not be 'all' and a specific interval.")
  }
  if (length(interval) > 1 & table == TRUE) {
    warning("Warning: A table can not be produced for multiple intervals.")
    table = FALSE
  }
  ci_colors <- c("green4", "royalblue", "red3", "red3","orange1", "purple3")
  for (i in interval){
    if (i == "Greenwood"){
      plot_Greenwood(n, lambdaT, lambdaC, alpha, show, table, value, new_plot = TRUE, col = ci_colors[1])
      par(new = TRUE)
    }
    if (i == "Exp-Greenwood"){
      plot_ExpGreenwood(n, lambdaT, lambdaC, alpha, show, table, value, new_plot = TRUE, col = ci_colors[2])
      par(new = TRUE)
    }
    if (i == "Log-Log" || i == "Log-log"){
      plot_LogLog(n, lambdaT, lambdaC, alpha, show, table, value, new_plot = TRUE, col = ci_colors[3])
      par(new = TRUE)
    }
    if (i == "Arcsine"){
      plot_Arcsine(n, lambdaT, lambdaC, alpha, show, table, value, new_plot = TRUE, col = ci_colors[5])
      par(new = TRUE)
    }
    if (i == "Peto"){
      plot_Peto(n, lambdaT, lambdaC, alpha, show, table, value, new_plot = TRUE, col = ci_colors[6])
      par(new = TRUE)
    }
    if (i == "all"){
      plot_all(n,lambdaT,lambdaC, alpha, show)
    }
  }
  #To prevent gray curve being plotted over color interval curves, plot interval curves without gray curves on same plot
  if (length(interval)>1 & show == TRUE){
    for (i in interval){
      if (i == "Greenwood"){
        plot_Greenwood(n, lambdaT, lambdaC, alpha, show = FALSE, table, value, new_plot = TRUE, col = ci_colors[1])
        par(new = TRUE)
      }
      if (i == "Exp-Greenwood"){
        plot_ExpGreenwood(n, lambdaT, lambdaC, alpha, show = FALSE, table, value, new_plot = TRUE, col = ci_colors[2])
        par(new = TRUE)
      }
      if (i == "Log-Log" || i == "Log-log"){
        plot_LogLog(n, lambdaT, lambdaC, alpha, show = FALSE, table, value, new_plot = TRUE, col = ci_colors[3])
        par(new = TRUE)
      }
      if (i == "Arcsine"){
        plot_Arcsine(n, lambdaT, lambdaC, alpha, show = FALSE, table, value, new_plot = TRUE, col = ci_colors[5])
        par(new = TRUE)
      }
      if (i == "Peto"){
        plot_Peto(n, lambdaT, lambdaC, alpha, show = FALSE, table, value, new_plot = TRUE, col = ci_colors[6])
        par(new = TRUE)
      }
    }
  }
  #Legend
  if (length(interval) > 1) {
    color = ci_colors[match(interval, interval_options)]
    legend("topleft", legend=interval, col=ci_colors[match(interval, interval_options)], lty=1, cex=0.8, bg = "white")
  }
  par(new = FALSE)
}

#pdf("name")
#km.coverage(2,1,2,0.2)
#dev.off()
#km.coverage(6, 9, 1,show = T, interval = c("Greenwood", "Exp-Greenwood", "Log-Log", "Arcsine"))

