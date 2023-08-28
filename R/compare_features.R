#' graph all features based on hotspot mutation burden
#'
#' @description
#' Creates a scatter plot of numerical features vs hotspot mutation burden and
#' calculates the pearson correlation.
#' Creates a violin plot for categorical features and calculates p-value based
#' on either Wilcox test or Kruskal-Wallis test.
#'
#'
#' @param i a numerical value based on which feature in feature list to pull
#' from
#' @param df a dataframe containing the hotspot mutation count, and values for
#' each feature for all samples in dataset
#' @param feature A list containing all the features.
#'
#' @return graph of comparison of feature vs hotspot mutation burden in all
#' samples
#'
#' @keywords internal
#'
#' @examples
#'
#' plot_features(1, df = sample_df)
#'
#' @noRd
#'

plot_features <- function(i, df, feature) {
  df <- na.omit(df)
  feature_df <- df[, feature]
  if (typeof(feature_df[, i]) %in% c("numeric", "integer", "double")) {
    plot <-
      ggplot2::ggplot(data = df, aes(x = feature_df[, i], y = Count)) +
      ggplot2::geom_point(color = "brown4", fill = "brown4") +
      ggpubr::stat_cor() +
      ggplot2::theme_classic() +
      ggplot2::xlab(feature[[i]]) +
      ggplot2::theme(
        axis.text.x = element_text(
          angle = 45,
          vjust = 0.5,
          hjust = 0.5,
          size = 12
        ),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title = element_text(size = 4),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 12)
      )
  }
  else{
    plot <-
      ggplot2::ggplot(df,
                      aes(
                        x = feature_df[, i],
                        y = Count,
                        group = feature_df[, i],
                        fill = feature_df[, i]
                      )) +
      ggplot2::geom_violin(trim = FALSE) +
      ggpubr::stat_compare_means() +
      ggplot2::scale_fill_manual(values = c(
        "brown4",
                "darksalmon",
                "tan",
                "darkorange",
                "white",
                "grey"
      )) +
      ggplot2::scale_color_manual(values = c(
        "brown4",
                "darksalmon",
                "tan",
                "darkorange",
                "white",
                "grey"
      )) +
      ggplot2::theme_classic() +
      ggplot2::xlab(feature[[i]]) +
      ggplot2::theme(
        axis.text.x = element_text(
          angle = 45,
          vjust = 1,
          hjust = 1,
          size = 12
        ),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title = element_text(size = 4),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none"
      )
  }
  return(plot)
}






#' hotspot comparison by additional features
#'
#' @description
#' This function performs an exploratory data analysis
#' comparing the relationship between user-input features to
#' hotspot mutation burden.
#' @details
#' This function is used to classify the features into
#' sequential features if values are numerical or classifies
#' them into #' categorical features. Sequential features are
#' compared to the mutation count using Pearson correlation.
#' Similarly, in categorical features either Wilcox or
#' Kruskal-Wallis test is used to compare between the groups
#' in the features based on the mutational count. Scatter plot
#' is used to represent the sequential features along with the
#' R and p-value from the pearson correlation. Violin plots are
#' used to plot the groups in the categorical data and Wilcox
#' or Kruskal-Wallis values are shown on the graph.
#'
#'
#' @param data A dataframe containing the clinical features and
#' the mutation count. Dataframe must contain columns with the following names:
#' "Chromosome" <-- Chromosome number where the mutation is located
#' "Position" <-- Genomic position number where the mutation is located
#' "Sample" <-- Unique ID for each sample in dataset
#' "Gene" <-- Name of the gene which mutation is located in (optional)

#' @param feature A list containing all the features.
#' @param regions a dataframe containing the chromosome, start
#' and end base pair position of each region of interest
#' @importFrom ggplot2 ggplot aes theme theme_classic
#' scale_color_manual element_text geom_point geom_smooth xlab annotate
#' ylab geom_violin scale_fill_manual
#' @importFrom ggpubr stat_cor stat_compare_means ggarrange
#' @importFrom utils combn
#' @importFrom stats wilcox.test
#' @importFrom gridExtra grid.arrange
#'
#'
#'
#' @return A grid of all the violin plots for the categorical
#' data and scatter plot for the sequential data.
#'
#'
#' @examples
#'
#' data("compSPOT_example_mutations")
#' data("compSPOT_example_regions")
#' features <- c("AGE", "SEX", "ADJUVANT_TX", "SMOKING_HISTORY",
#' "TUMOR_VOLUME", "KI_67")
#' compare_features(data = compSPOT_example_mutations,
#' regions = compSPOT_example_regions, feature = features)
#'
#' @export
#'

compare_features <- function(data, regions, feature) {
  if (base::isFALSE(is.data.frame(data))) {
    stop("Input data must be in the form of data.frame")
  }
  if (base::isFALSE(is.data.frame(regions))) {
    stop("Input regions must be in the form of data.frame")
  }
  if (base::isFALSE("Chromosome" %in% colnames(data))) {
    stop("input data must contain column 'Chromosome'")
  }
  if (base::isFALSE("Position" %in% colnames(data))) {
    stop("input data must contain column 'Position'")
  }
  if (base::isFALSE("Sample" %in% colnames(data))) {
    stop("input data must contain column 'Sample'")
  }
  if (base::isFALSE("Chromosome" %in% colnames(regions))) {
    stop("input regions must contain column 'Chromosome'")
  }
  if (base::isFALSE("Lowerbound" %in% colnames(regions))) {
    stop("input regions must contain column 'Lowerbound'")
  }
  if (base::isFALSE("Upperbound" %in% colnames(regions))) {
    stop("input regions must contain column 'Upperbound'")
  }
  if (base::isFALSE(is.numeric(data$Position))) {
    stop("column 'Position' must be in the form of numeric")
  }
  if (base::isFALSE(is.numeric(regions$Lowerbound))) {
    stop("column 'Lowerbound' must be in the form of numeric")
  }
  if (base::isFALSE(is.numeric(regions$Upperbound))) {
    stop("column 'Upperbound' must be in the form of numeric")
  }
  count_mutations <- function(x) {
    sub <- subset(data, Sample == x)
    count <- 0
    for (i in seq_len(nrow(regions))) {
      sub <-
        subset(
          sub,
          Chromosome == regions$Chromosome[[i]] &
            Position %in% regions$Lowerbound[[i]]:regions$Upperbound[[i]]
        )
      count <- count + nrow(sub)
    }
    return(count)
  }
  count.ls <- lapply(unique(data$Sample), count_mutations)
  count.df <-
    data.frame("Sample" = unique(data$Sample),
               "Count" = unlist(count.ls))
  clin.df <- data[, c("Sample", feature)]
  clin.df <- unique(clin.df)
  sample_df <- merge(clin.df, count.df, by = "Sample")
  plot_ls <-
    lapply(seq_along(feature),
           plot_features,
           df = sample_df,
           feature = feature)
  plots <-
    ggpubr::ggarrange(plotlist = plot_ls,
                      nrow = 2,
                      ncol = ceiling(length(plot_ls) / 2))
  return(plots)
}

