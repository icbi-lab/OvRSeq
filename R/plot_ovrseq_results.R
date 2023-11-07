#' Plot Marginal Distributions for OvRSeq Results
#'
#' This function takes the results from OvRSeq, extracts specified variables from
#' the column data, and creates a ggplot with marginal histograms to show the
#' distribution of x and y variables separately.
#'
#' @param se A SummarizedExperiment object containing results from OvRSeq.
#' @param x_var The name of the variable in colData(se) to use for the x-axis.
#' @param y_var The name of the variable in colData(se) to use for the y-axis.
#' @param color_var The name of the variable in colData(se) to use for coloring the points.
#' @return A ggplot object with added marginal histograms.
#' @importFrom ggplot2 ggplot aes aes_string geom_point theme theme_bw
#' @importFrom ggExtra ggMarginal
#' @examples
#' # Assuming `se` is a SummarizedExperiment object with relevant data
#' # plot_ggmarginal(se, "variable1", "variable2", "groupingVar")
#'
#' @export
plot_ggmarginal <- function(se, x_var, y_var, color_var) {
  # Check that the provided variables are in colData
  required_vars <- c(x_var, y_var, color_var)
  if (!all(required_vars %in% colnames(colData(se)))) {
    stop("Not all specified variables are present in colData of the SummarizedExperiment object.")
  }

  # Convert the desired column data to a data frame for ggplot
  plot_data <- as.data.frame(colData(se)[, required_vars])
  colnames(plot_data) <- required_vars
  # Create the base ggplot
  # Convert string to symbols
  x_sym <- sym(x_var)
  y_sym <- sym(y_var)
  color_sym <- sym(color_var)

  # Use aes() with !! to force evaluation of the symbols
  p <- ggplot(plot_data, aes(x = !!x_sym, y = !!y_sym, color = !!color_sym)) +
    geom_point() +
    theme(legend.position = "left")

  # Add theme
  p <- p + theme_bw()

  # Add marginal histograms
  p <- ggMarginal(p, type="histogram", groupColour=TRUE, groupFill=TRUE)

  # Return the plot
  return(p)
}
