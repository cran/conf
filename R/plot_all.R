plot_all <- function(n, lambdaT, lambdaC, alpha, show) {
  # Initialize the plot
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = 'p', ylab = 'c(p)', las = 1, font.axis = 5, cex.axis = 0.8)

  abline(h = 1 - alpha, lty = 2, col = 'red')  # Common reference
  ci_colors <- c("green4", "royalblue", "red3", "orange1", "purple3")
  # Plot each method on the same graph
  plot_Greenwood(n, lambdaT, lambdaC, alpha, show, table = FALSE, value = 0.7, new_plot = FALSE, col=ci_colors[1])
  plot_ExpGreenwood(n, lambdaT, lambdaC, alpha, show, table = FALSE, value = 0.7, new_plot = FALSE, col=ci_colors[2])
  plot_LogLog(n, lambdaT, lambdaC, alpha, show, table = FALSE, value = 0.7, new_plot = FALSE, col=ci_colors[3])
  plot_Arcsine(n, lambdaT, lambdaC, alpha, show, table = FALSE, value = 0.7, new_plot = FALSE, col=ci_colors[4])
  plot_Peto(n, lambdaT, lambdaC, alpha, show, table = FALSE, value = 0.7, new_plot = FALSE, col=ci_colors[5])

  if (show == TRUE){
    # Plot each method on the same graph over the already drawn gray lines
    plot_Greenwood(n, lambdaT, lambdaC, alpha, show = FALSE, table = FALSE, value = 0.7, new_plot = FALSE, col=ci_colors[1])
    plot_ExpGreenwood(n, lambdaT, lambdaC, alpha, show = FALSE, table = FALSE, value = 0.7, new_plot = FALSE, col=ci_colors[2])
    plot_LogLog(n, lambdaT, lambdaC, alpha, show = FALSE, table = FALSE, value = 0.7, new_plot = FALSE, col=ci_colors[3])
    plot_Arcsine(n, lambdaT, lambdaC, alpha, show = FALSE, table = FALSE, value = 0.7, new_plot = FALSE, col=ci_colors[4])
    plot_Peto(n, lambdaT, lambdaC, alpha, show = FALSE, table = FALSE, value = 0.7, new_plot = FALSE, col=ci_colors[5])
  }
    # Add a legend
  legend("topleft", legend=c("Greenwood", "Exp-Greenwood", "Log-Log", "Arcsine", "Peto"),
         col=ci_colors, lty=1, cex=0.8, bg = "white")
  text(0.53, -0.15, 'xlab')
  text(-0.2, 0.5, 'ylab')
}
