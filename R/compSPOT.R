#' compSPOT
#'
#' @description
#' It is well known that numerous clones of cells sharing common mutations exist within
#' cancer, precancer, and even clinically normal appearing tissues. The frequency and location of
#' these mutations may aid in the prediction of cancer risk of certain individuals. It has also been
#' well established that certain genomic regions have increased sensitivity to acquiring mutations.
#' Mutation-sensitive genomic regions may therefore be used as markers for prediction of cancer risk.
#' This package contains multiple functions for the establishment of significantly mutated hotspots,
#' comparison of hotspot mutation burden between sub-groups, and exploratory data analysis of the
#' correlation between hotspot mutation burden, and personal risk factors for cancer such as age, gender,
#' and history of carcinogen exposure. This package aims to allow users to identify robust genomic markers
#' which may serve as markers of cancer risk.
#'
#'
#' @importFrom stats ks.test
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 stat_ecdf
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 geom_smooth
#' @importFrom plotly ggplotly
#' @importFrom magrittr %>%
#' @importFrom plotly plot_ly
#' @importFrom plotly add_trace
#' @importFrom plotly layout
#' @importFrom ggpubr stat_cor
#' @importFrom ggpubr ggarrange
#' @importFrom utils combn
#' @importFrom stats wilcox.test
#' @importFrom gridExtra grid.arrange
#' @param data a dataframe containing the chromosome, base pair position, and optionally gene name of each mutation
#' @param regions a dataframe containing the chromosome, start and end base pair position of each region of interest
#' @param pvalue the p-value cutoff for included hotspots
#' @param threshold the cutoff empirical distribution for Kolmogorov-Smirnov test
#' @param include_genes true or false whether gene names are included in regions dataframe
#' @param rank true or false whether regions dataframe is already ranked and includes mutation count of total dataset
#' @param name1 a string containing the name of one group for the comparison
#' @param name2 a string containing the name of the second group for the comparison
#' @param feature A dataframe containing all the features.
#'
#'@examples
#'
#' data("example_mutations")
#' data("example_regions")
#'
#' sig.spots(data = example_mutations, regions = example_regions,
#' pvalue = 0.05, threshold = 0.2, include_genes = TRUE, rank = TRUE)
#'
#' group.spot(data = example_mutations, regions = example_regions, pval = 0.05, threshold = 0.4
#' name1 = "High-Risk", name2 = "Low-Risk", include_genes = TRUE)
#'
#' features <- c("AGE", "SEX", "ADJUVANT_TX", "SMOKING_HISTORY", "TUMOR_VOLUME", "KI_67")
#' feature.spot(data = example_mutations, regions = example_regions, feature = features)
