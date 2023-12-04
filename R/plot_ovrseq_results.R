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
#' @importFrom ggplot2 ggplot aes aes_string geom_point theme theme_bw guides guide_legend
#' @importFrom ggExtra ggMarginal
#' @importFrom rlang sym
#' @examples
#' # Assuming `se` is a SummarizedExperiment object with relevant data
#' # plot_ggmarginal(se, "variable1", "variable2", "groupingVar")
#'
#' @export
plot_ggmarginal <- function(se, x_var, y_var, color_var = NA) {
  # Check that the provided variables are in colData
  required_vars <- c(x_var, y_var, if (!is.na(color_var)) color_var)
  if (!all(required_vars %in% c(rownames(assay(se)),colnames(colData(se))))) {
    stop("Not all specified variables are present in colData of the SummarizedExperiment object.")
  }

  # Convert the desired column data to a data frame for ggplot
  data <- cbind(colData(se), t(assay(se)))
  plot_data <- as.data.frame(data[, required_vars])
  colnames(plot_data) <- required_vars
  # Create the base ggplot
  # Convert string to symbols
  x_sym <- sym(x_var)
  y_sym <- sym(y_var)

  if (!is.na(color_var)){
    color_sym <- sym(color_var)
  }

  # Use aes() with !! to force evaluation of the symbols
  if (!is.na(color_var)){
    p <- ggplot(plot_data, aes(x = !!x_sym, y = !!y_sym, color = !!color_sym))
  } else{
    p <- ggplot(plot_data, aes(x = !!x_sym, y = !!y_sym))
  }
  p <- p + geom_point() +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(nrow=1, byrow=TRUE))


  # Add theme
  p <- p + theme_bw() #+ coord_fixed(ratio = 1.5)


  # Add marginal histograms
  if (!is.na(color_var)){
    p <- p + theme(legend.position = "bottom") +
      guides(color=guide_legend(nrow=1, byrow=TRUE))
    p <- ggMarginal(p, type="histogram", groupColour=TRUE, groupFill=TRUE)
  } else {
    p <- ggMarginal(p, type="histogram", groupColour=F, groupFill=F)

  }

  # Return the plot
  return(p)
}


#' Plot one sample for report
#'
#' This function creates a scatter plot with marginal histograms for specified variables from a `SummarizedExperiment` object. It combines the provided `se` object with the TCGA_OV dataset, plots specified variables, and highlights a specific sample from `se` with a distinct style.
#'
#' @param se A `SummarizedExperiment` object.
#' @param x_var A string representing the variable from `se` to be plotted on the x-axis (default is "C1QA").
#' @param y_var A string representing the variable from `se` to be plotted on the y-axis (default is "CD8A").
#' @param color_var A string representing the variable from `se` to be used for coloring points (default is "BRCAness_Prob").
#'
#' @details
#' The function first combines the provided `se` object with the TCGA_OV dataset. It then creates a scatter plot for the specified x and y variables, with points colored based on the specified `color_var`. A specific sample from `se` is highlighted in the plot with a black border and white fill, distinct from other data points.
#'
#' The function assumes that the TCGA_OV dataset and `se` object have compatible structures and the required columns. The user needs to ensure that the TCGA_OV dataset is loaded and accessible.
#'
#' @return A `ggplot` object representing the scatter plot with marginal histograms.
#'
#' @examples
#' # se is a pre-loaded SummarizedExperiment object
#' plot_ggmarginal(se, "C1QA", "CD8A", "BRCAness_Prob")
#'
#' @importFrom SummarizedExperiment colData assay
#' @importFrom ggplot2 ggplot aes aes_string geom_point theme theme_bw guides guide_legend scale_y_discrete
#' @importFrom ggExtra ggMarginal
#' @export
plot_ggmarginal_sample <- function(se, x_var = "C1QA", y_var = "CD8A", color_var = "BRCAness") {
  # Load TCGA_OV dataset (modify this as per your actual data loading method)
  tcga_ov <- load_TCGA_OV()  # Placeholder function

  # Combine TCGA_OV dataset with the provided se object
  combined_data_tcga <- cbind(colData(tcga_ov), t(assay(tcga_ov)))
  combined_data_se <- rbind(cbind(colData(se), t(assay(se))))

  # Check that the provided variables are in the combined data
  required_vars <- c(x_var, y_var, color_var)
  if (!all(required_vars %in% colnames(combined_data_tcga))) {
    stop("Not all specified variables are present in the combined dataset.")
  }

  # Prepare the plot data
  plot_data_tcga <- as.data.frame(combined_data_tcga[, required_vars])
  plot_data_sample <- as.data.frame(combined_data_se[,c(x_var, y_var)])
  colnames(plot_data_tcga) <- required_vars
  colnames(plot_data_sample) <- c(x_var, y_var)

  # Create the base ggplot
  p <- ggplot(plot_data_tcga, aes_string(x = x_var, y = y_var, color = color_var)) +
    geom_point() +
    theme(legend.position = "left")

  # Highlight the specific sample from se
  p <- p + geom_point(data = plot_data_sample, aes_string(x = x_var, y = y_var),
                      color = "black", fill = "white", shape = 21, size = 3, stroke = 2)

  # Add theme and marginal histograms
  p <- p + theme_bw()
  p <- p + theme(legend.position = "bottom") +
    guides(color=guide_legend(nrow=1, byrow=TRUE)) +
    theme(legend.box.spacing = unit(0, "pt")) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))

  p <- ggMarginal(p, type = "histogram", groupColour = TRUE, groupFill = TRUE)

  for (obj in c("TCGA_OV")) {
    if (exists(obj, envir = .GlobalEnv)) {
      rm(list = obj, envir = .GlobalEnv)
    }
  }
  return(p)
}

#' Plot Deconvolution Data for a Given Sample
#'
#' Creates a bar plot representing the estimated immune cell infiltrate fractions for a specified sample
#' using the first available deconvolution method from a provided list.
#'
#' @param se A `SummarizedExperiment` object containing the deconvolution data in its `colData`.
#' @param sample_id The identifier of the sample for which the plot will be generated.
#' @param deconvolution_methods A vector of strings representing the deconvolution methods to check
#' in the `colData` of the `SummarizedExperiment` object. The function will use the first method that
#' has corresponding data available.
#'
#' @return A ggplot object representing the estimated immune cell infiltrate fractions.
#'
#' @details
#' The function iterates through the provided list of deconvolution methods and checks if there are
#' any columns in the `colData` of the `SummarizedExperiment` object corresponding to these methods.
#' It then uses the data from the first available method to create a bar plot. If no data is found
#' for any of the provided methods, the function stops with an error.
#'
#' @examples
#' # Assuming se is a pre-loaded SummarizedExperiment object and sample_id is known
#' immunedeconv_methods <- c("quantiseq", "timer", "mcp_counter")
#' sample_id <- "sample123"  # Replace with an actual sample ID
#' p <- plot_deconvolution_data(se, sample_id, immunedeconv_methods)
#' print(p)
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw ylab xlab ggtitle scale_y_discrete
#' @importFrom stringr str_wrap
#' @export
plot_deconvolution_data <- function(se, sample_id, deconvolution_methods=c("quantiseq", "timer", "mcp_counter")) {
  # Extract data for the specified sample
  data <- colData(se)[sample_id, ]

  plot_data <- NULL
  selected_method <- NULL

  # Iterate through deconvolution methods
  for (method in deconvolution_methods) {
    method_columns <- grepl(method, names(data), ignore.case = TRUE)

    if (any(method_columns)) {
      selected_method <- method
      method_data <- data[method_columns]

      # Transform 'method_data' into a format suitable for plotting
      plot_data <- as.data.frame(t(method_data))
      colnames(plot_data) <- c("Row", "group_name", "value")
      plot_data$Metric <- unlist(lapply(plot_data$group_name,
                                        function(x) stringr::str_replace_all(x, pattern = paste0("\\|", method), replacement = "")))
      break
    }
  }

  # If no method columns were found
  if (is.null(plot_data)) {
    stop("No deconvolution data found for the specified methods.")
  }

  # Filter out specific categories if needed
  plot_data <- plot_data[plot_data$Metric != "uncharacterized cell", ]

  # Create a ggplot
  p <- ggplot(plot_data, aes(x = value, y = Metric)) +
    geom_bar(stat = "identity", color = "black", fill ="#991915") +
    theme_bw() +
    ylab("") +
    xlab("Cell fraction") +
    ggtitle("") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20))

  return(p)
}


#' Plot Immune Signature Scores for a Single Sample
#'
#' This function generates a bar plot of selected immune signature scores for a specified sample from a `SummarizedExperiment` object. It filters the data to include specific immune-related metrics and presents them with meaningful names.
#'
#' @param se A `SummarizedExperiment` object containing immune-related metrics in its `colData`.
#' @param sample_id A character string specifying the ID of the sample to be plotted.
#'
#' @return A `ggplot` object representing the bar plot of immune signature scores for the specified sample.
#'
#' @details
#' The function extracts data for the given `sample_id` and selects a subset of columns related to immune signatures. These columns are renamed for better readability in the plot. The data is then transformed into a long format suitable for plotting with `ggplot2`. The bar plot visualizes the immune pathway or signature scores, differentiated by the analysis type (Average or GSVA) in separate facets.
#'
#' @examples
#' # se is a pre-loaded SummarizedExperiment object
#' sample_id <- "sample123"
#' plot_immune_signature_one_sample(se, sample_id)
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes geom_bar facet_grid theme_bw ylab xlab ggtitle scale_y_discrete xlim
#' @importFrom stringr str_wrap
#' @export
plot_immune_signature_one_sample <- function(se, sample_id) {
  if (!sample_id %in% rownames(colData(se))) {
    stop("Specified sample_id not found in SummarizedExperiment object.")
  }

  # Extract data for the specified sample
  data <- colData(se)[sample_id, ]
  tcgaStats <- load_tcgaStats()

  # Filter for columns containing "quantiseq"
  sign_columns <- c("IFNG_Ayers", "Inflamed" ,"BRCAness_pos", "BRCAness_neg",
                    "T-cell exhaustion", "T-cell exclusion","Interferon alpha response",
                    "Immunological Constant of Rejection (ICR)")

  nice_name <- c("IFNG", "Inflamed", "BRCAness+","BRCAness-","T cell exhaustion",
                  "T cell exclusion", "Interferon alpha response","Immunological Constant of Rejection")
  data <- data[sign_columns]
  tcgaStats <- tcgaStats[sign_columns,]
  #colnames(data) <- nice_name
  # Perform min-max normalization using lapply
  normalizedData <- as.data.frame(lapply(names(data), function(col) {
    minVal <- tcgaStats[tcgaStats$Feature == col, "Min"]
    maxVal <- tcgaStats[tcgaStats$Feature == col, "Max"]

    new_value <- (data[[col]] - minVal) / (maxVal - minVal)
    if(new_value < 0){
      new_value <- 0
    } else if (new_value > 1){
      new_value <- 1
    }
    return(new_value)
  }))

  # Set the row names to match the original data
  colnames(normalizedData) <- colnames(data)

  # Transform 'quantiseq_data' into a format suitable for plotting
  plot_data <- as.data.frame(t(normalizedData))
  plot_data$group_name <- rownames(plot_data)
  colnames(plot_data) <- c( "value","group_name" )
  plot_data$Metric <- plot_data$group_name
  #plot_data$Analysis <- unlist(lapply(plot_data$Metric, function(x){
  #  ifelse(x %in% c("T cell exhaustion","IFNG"), "Average","GSVA")}))

  # Create a ggplot
  p <- ggplot(plot_data, aes(x = value, y = Metric)) +
    geom_bar(stat = "identity", color = "black", fill ="#175d92") +
    theme_bw()  +
    ggtitle("") +
  ylab("") + xlim(c(0,1)) +
    xlab("Enrichment score*")  + scale_y_discrete(labels = function(x) str_wrap(x, width = 20))

  return(p)
}

