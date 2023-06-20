#' hotspot group comparison with plot
#'
#' @description
#' Based on a dataframe of hotspots, this function groups the hotspots into two different groups and overlays interactive
#' violin plots based on mutation count for comparison.
#'
#' @details
#' This function compares the mutation frequency per sample for each hotspot in a violin plot.
#'
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
#' @importFrom gridExtra grid.arrange
#' @importFrom magrittr %>%
#' @param plot_data a dataframe containing the hotspot, mutation count, and associated groups
#' @param group1 a string containing the name of one group for the comparison
#' @param group2 a string containing the name of the second group for the comparison
#'
#' @return An interactive violin plot displaying the comparison of mutation counts for each group.
#'
#' @keywords internal
#'
#' @examples
#'
#' data("plot_data")
#' fig <- plot.spot(plot_data = plot_data, group1 = "High-Risk", group2 = "Low-Risk")
#'
#' @noRd


plot.spot <- function(plot_data, group1, group2) {
  hotspot_list <- unique(plot_data$Hotspot)
  for (spot in hotspot_list){
    hotspot <- subset(plot_data, Hotspot == spot)
    group1_hotspot <- subset(hotspot,Group== group1)
    group2_hotspot <- subset(hotspot,Group== group2)
    group1_list_hotspot <- group1_hotspot$Mutation_Count
    group2_list_hotspot <- group2_hotspot$Mutation_Count
    fig <- plot_data %>%
      plotly::plot_ly(type = 'violin')
    fig <- fig %>%
      plotly::add_trace(
        x = ~Hotspot[plot_data$Group == group1],
        y = ~Mutation_Count[plot_data$Group == group1],
        legendgroup = group1,
        scalegroup = group1,
        name = group1,
        box = list(
          visible = TRUE
        ),
        meanline = list(
          visible = TRUE
        ),
        color = I("brown4")
      )
    fig <- fig %>%
      plotly::add_trace(
        x = ~Hotspot[plot_data$Group == group2],
        y = ~Mutation_Count[plot_data$Group == group2],
        legendgroup = group2,
        scalegroup = group2,
        name = group2,
        box = list(
          visible = TRUE
        ),
        meanline = list(
          visible = TRUE
        ),
        color = I("darksalmon")
      )
    fig <- fig %>%
      plotly::layout(
        xaxis = list( title = "Hotspot"),
        yaxis = list( title = "Mutation Count",
                      zeroline = FALSE
        ),
        violinmode = 'group'
      )
    plot2_list <- c()
    for (i in seq_len(length(unique(plot_data$Hotspot)))){
      spot <- unique(plot_data$Hotspot)[[i]]
      plot_data.sub <- subset(plot_data, Hotspot == spot)
      plot2 <- ggplot2::ggplot(data = plot_data.sub, ggplot2::aes(x=Mutation_Count, group = Group, col = Group)) +
        ggplot2::stat_ecdf(geom = "point") +
        ggplot2::theme_classic() +
        ggplot2::scale_color_manual(values = c("brown4", "darksalmon"))
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=0.5, size = 10),
                     axis.title.y = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 4),
                     axis.text.y = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(hjust = 0.5, size = 12))
      plot2_list[[i]] <- plot2
    }
  }
  output <- c()
  output[[1]] <- fig
  output[[2]] <- ggpubr::ggarrange(plotlist = plot2_list,nrow=2, ncol = ceiling(length(plot2_list)/2))
  return(output)
}

#' hotspot comparison by group
#'
#' @description
#' This function compares the mutation frequency of a panel of genomic regions between two sub-groups.
#'
#'
#' @details
#' This function creates a list of mutation frequency per unique sample for each genomic regions separated
#' based on specified sub-groups. The regions with significant differences in mutation distribution are
#' calculated using a Kolmogorov-Smirnov test. The difference in mutation frequency is output in a violin plot.
#'
#' @importFrom stats ks.test
#' @importFrom magrittr %>%
#' @import plotly
#' @param data a dataframe containing the chromosome, base pair position, and optionally gene name of each mutation
#' @param regions a dataframe containing the chromosome, start and end base pair position of each region of interest
#' @param pval a threshold p-value for Kolmogorov-Smirnov test
#' @param threshold the cutoff empirical distribution for Kolmogorov-Smirnov test
#' @param name1 a string containing the name of one group for the comparison
#' @param name2 a string containing the name of the second group for the comparison
#'
#' @return A dataframe containing the genomic regions with significant mutation frequency and summary plots for
#' mutation comparison between groups
#'
#' @examples
#'
#' data("example_mutations")
#' data("example_regions")
#' group.spot(data = example_mutations, regions = example_regions, pval = 0.05, threshold = 0.4,
#' name1 = "High-Risk", name2 = "Low-Risk", include_genes = TRUE)
#'
#' @export
#'


group.spot <- function(data, regions, pvalue, threshold, name1, name2, include_genes){
  if (base::isFALSE(is.data.frame(data))){stop("Input data must be in the form of data.frame")}
  if (base::isFALSE(is.data.frame(regions))){stop("Input regions must be in the form of data.frame")}
  if (base::isFALSE("Chromosome" %in% colnames(data))){stop("input data must contain column 'Chromosome'")}
  if (base::isFALSE("Position" %in% colnames(data))){stop("input data must contain column 'Position'")}
  if (base::isFALSE("Sample" %in% colnames(data))){stop("input data must contain column 'Sample'")}
  if (base::isFALSE("Group" %in% colnames(data))){stop("input data must contain column 'Group'")}
  if (base::isFALSE(length(unique(data$Group)) == 2)){stop("input data must have two unique group labels")}
  if (base::isFALSE("Chromosome" %in% colnames(regions))){stop("input regions must contain column 'Chromosome'")}
  if (base::isFALSE("Lowerbound" %in% colnames(regions))){stop("input regions must contain column 'Lowerbound'")}
  if (base::isFALSE("Upperbound" %in% colnames(regions))){stop("input regions must contain column 'Upperbound'")}
  if (base::isFALSE(is.numeric(data$Position))){stop("column 'Position' must be in the form of numeric")}
  if (base::isFALSE(is.numeric(regions$Lowerbound))){stop("column 'Lowerbound' must be in the form of numeric")}
  if (base::isFALSE(is.numeric(regions$Upperbound))){stop("column 'Upperbound' must be in the form of numeric")}
  if (base::isTRUE(include_genes) & base::isFALSE("Gene" %in% colnames(data)))
  {stop("input data must contain column 'Gene'")}
  if (base::isTRUE(include_genes) & base::isFALSE(is.character(data$Gene))){
    stop("column 'Gene' must be in the form of character")
  }
  if (pvalue > 1){stop("pvalue cannot be greater than 1")}
  if (threshold > 1){stop("threshold cannot be greater than 1")}
  if (pvalue < 0){stop("pvalue cannot be less than 0")}
  if (threshold < 0){stop("threshold cannot be less than 0")}
  plot_data <- data.frame()
  all_diff <- data.frame()
  group1 <- subset(data, Group == name1)
  group2 <- subset(data, Group == name2)
  for (i in 1:nrow(regions)){
    count_1 <- c()
    count_2 <- c()
    for (samp in unique(group1$Sample)){
      group1_sub <- subset(group1, Sample == samp)
      count_1 <- c(count_1, nrow(subset(group1_sub, Chromosome == regions$Chromosome[[i]]
                                        & Position %in% regions$Lowerbound[[i]]:regions$Upperbound[[i]])))
    }
    for (samp in unique(group2$Sample)){
      group2_sub <- subset(group2, Sample == samp)
      count_2 <- c(count_2, nrow(subset(group2_sub, Chromosome == regions$Chromosome[[i]]
                                        & Position %in% regions$Lowerbound[[i]]:regions$Upperbound[[i]])))
    }
    if (include_genes == TRUE){
      spot_data <- data.frame("Hotspot" = paste(regions$Gene[[i]], " ", regions$Chromosome[[i]],
                                                ":", regions$Lowerbound[[i]], "-",
                                                regions$Upperbound[[i]], sep = ""),
                              "Group" = c(rep(name1, length(count_1)),
                                          rep(name2, length(count_2))),
                              "Mutation_Count" = c(count_1, count_2))
    }
    if (include_genes == FALSE){
      spot_data <- data.frame("Hotspot" = paste(regions$Chromosome[[i]],
                                                ":", regions$Lowerbound[[i]], "-",
                                                regions$Upperbound[[i]], sep = ""),
                              "Group" = c(rep(name1, length(count_1)),
                                          rep(name2, length(count_2))),
                              "Mutation_Count" = c(count_1, count_2))

    }
    plot_data <- rbind(plot_data, spot_data)
    if (sum(c(count_1, count_2)) > 0){
      ks <- ks.test(count_1, count_2, alternative = "two.sided")
      if (ks$p.value < pvalue){
        if (sum(count_1) > sum(count_2)){high_group <- "Group 1"}
        else{high_group <- "Group 2"}
        if (include_genes == TRUE){
          diff_hotspots <- data.frame("D" = ks$statistic, "pval" = ks$p.value, "Higher Group" = high_group, "Hotspot" = paste(regions$Gene[[i]], " ", regions$Chromosome[[i]],
                                                                                                                              ":", regions$Lowerbound[[i]], "-",regions$Upperbound[[i]], sep = ""))
        }
        if (include_genes == FALSE){
          diff_hotspots <- data.frame("D" = ks$statistic, "pval" = ks$p.value, "Higher Group" = high_group, "Hotspot" = paste(regions$Chromosome[[i]],
                                                                                                                              ":", regions$Lowerbound[[i]], "-",regions$Upperbound[[i]], sep = ""))
        }
        all_diff <- rbind(all_diff, diff_hotspots)

      }
    }
  }
  plot <- plot.spot(plot_data, group1 = name1, group2 = name2)
  return_items <- c()
  return_items[[1]] <- plot_data
  return_items[[2]] <- plot[[1]]
  return_items[[3]] <- plot[[2]]
  return(return_items)
}

