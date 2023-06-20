#' genomic region ranker
#'
#' @description
#' ranks all genomic regions from highest to lowest mutation frequency in given dataset
#'
#' @param data A dataframe containing the location of each mutation.
#' @param regions a dataframe containing the location of gene region of interest
#' @param include_genes true or false whether gene names are included in regions dataframe
#'
#' @keywords internal
#'
#' @return A dataframe containing the ranked genomic regions from highest to lowest mutation frequency
#'
#'
#' @examples
#'
#' data("mutation_data")
#' spot.rank(data = mutation_data, regions = ea_regions, include_genes = TRUE)
#'
#'@noRd


spot.rank <- function(data, regions, include_genes){
  ranked_regions <- data.frame()
  for (i in seq_len(nrow(regions))){
    data_sub <- subset(data, Chromosome == regions$Chromosome[[i]] &
                         Position %in% regions$Lowerbound[[i]]:regions$Upperbound[[i]])
    if (include_genes == TRUE){
      ranked_row <- data.frame("Chromosome" = regions$Chromosome[[i]], "Lowerbound" = regions$Lowerbound[[i]],
                               "Upperbound" = regions$Upperbound[[i]], "Count" = nrow(data_sub), "Gene" = regions$Gene[[i]])
    }
    if (include_genes == FALSE){
      ranked_row <- data.frame("Chromosome" = regions$Chromosome[[i]], "Lowerbound" = regions$Lowerbound[[i]],
                               "Upperbound" = regions$Upperbound[[i]], "Count" = nrow(data_sub))
    }
    ranked_regions <- rbind(ranked_regions, ranked_row)
  }
  ranked_regions <- ranked_regions[order(-ranked_regions$Count), ]
  return(ranked_regions)
}


#' sig.spot visualization
#'
#' @description
#' Based on a panel of genomic regions, this function visualizes the regions which are found to have
#' significance based on percentage of samples with hotspot mutations.
#'
#' @details
#' This function begins by displaying a plot with a ranked number of hotspot based on number of mutations
#' on the x-axis versus percentage of all samples with mutations in hotspot on the y-axis. Regions of high
#' significance are ranked on the left side of the plot. Based on the type, the points plotted as red
#' dots are the hotspots which were found to have a significant number of mutations whereas the pink were
#' found to not be significant. An interactive label of each point is also available. A second plot is also
#' created using the same data used in the Kolmogorov-Smirnov calculation to show a visualization of the ECDF
#' of hotspots vs non-hotspots.
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
#' @importFrom plotly ggplotly
#' @param df1 a dataframe containing the labeled region number, type of region, percentage of samples with
#' hotspot mutation
#' @param df2 a dataframe containing the data required for calculation of ECDF
#'
#' @return A figure organizing the region number in ranked order versus the percentage of the sample with hotspot
#' mutation which also denotes the label of the data point is displayed when hovering over the item. A second
#' figure providing a visualization of the ECDF for hotspot vs non-hotspots.
#'
#' @keywords internal
#' @examples
#'
#' data("regions")
#' data("ecdf_df")
#' plot <- plot.sigspot(df1 = region, df2 = ecdf_df)
#'
#'@noRd


plot.sigspot <- function(df1, df2){
  plot1 <- ggplot2::ggplot(df1, aes(y = percent, x = number,text = paste("Label",Label))) +
    ggplot2::geom_point(data = subset(df1, type == "Hotspot"),
                        color = "brown4", fill = "brown4", alpha = 1, shape = 21, size = 1.5) +
    ggplot2::geom_point(data = subset(df1, type == "Non-hotspot"),
                        color = "darksalmon", fill = "darksalmon", alpha = 1, shape = 21, size = 1.5) +
    ggplot2::geom_vline(xintercept = (nrow(subset(df1, type == "Hotspot")) + 1), linetype = "longdash",
                        linewidth = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Percentage of Samples \n with Hotspot Mutation") +
    ggplot2::xlab("Ranked Region Number") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 12),
                   axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title = element_text(size = 4),
                   axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 12))
  plot1 <- plotly::ggplotly(plot1, tooltip =  c("text"))
  plot2 <- ggplot2::ggplot(data = df2, aes(x=Count, group = Group, col = Group)) +
    ggplot2::stat_ecdf(geom = "point") +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = c("brown4", "darksalmon"))
  ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 12),
                 axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title = element_text(size = 4),
                 axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 12))
  plot2 <- plotly::ggplotly(plot2, tooltip =  c("text"))
  output <- c()
  output[[1]] <- plot1
  output[[2]] <- plot2
  return(output)
}


#' significant hotspot calculator
#'
#' @description
#' Based on a panel of genomic regions, this function calculates the regions which are found to have
#' significantly higher mutation frequency compared to less mutated regions.
#'
#' @details
#' This function begins by measuring the mutation frequency for each unique sample for each provided genomic region.
#' Beginning with the top ranked hotspot, a Kolmogorov-Smirnov test is preformed on the mutation frequency of the
#' top genomic region compared to the normalized mutation frequency of all the lower-ranked regions. This
#' continues, then running the Kolmogorov-Smirnov test for the normalized mutation frequency of the top 2 genomic
#' regions compared to the normalized mutation frequency of all lower-ranked regions.This process repeats itself,
#' continuously adding an additional genomic regions each time until either the set p-value or empirical distribution
#' threshold is not met. Once this cutoff has been reached, an established list of mutation hotspots is provided.
#'
#' @importFrom stats ks.test
#' @param data a dataframe containing the chromosome, base pair position, and optionally gene name of each mutation
#' @param regions a dataframe containing the chromosome, start and end base pair position of each region of interest
#' @param pvalue the p-value cutoff for included hotspots
#' @param threshold the cutoff empirical distribution for Kolmogorov-Smirnov test
#' @param include_genes true or false whether gene names are included in regions dataframe
#' @param rank true or false whether regions dataframe is already ranked and includes mutation count of total dataset
#'
#' @return A dataframe containing the genomic regions with significant mutation frequency and plots summarizing
#' significant hotspots
#'
#'
#' @examples
#'
#' data("example_mutations")
#' data("example_regions")
#' significant_spots <- sig.spots(data = example_mutations, regions = example_regions,
#' pvalue = 0.05, threshold = 0.2, include_genes = TRUE, rank = TRUE)
#'
#'
#' @export
#'


sig.spots <- function(data, regions, pvalue, threshold, include_genes, rank){
  if (base::isFALSE(is.data.frame(data))){stop("Input data must be in the form of data.frame")}
  if (base::isFALSE(is.data.frame(regions))){stop("Input regions must be in the form of data.frame")}
  if (base::isFALSE("Chromosome" %in% colnames(data))){stop("input data must contain column 'Chromosome'")}
  if (base::isFALSE("Position" %in% colnames(data))){stop("input data must contain column 'Position'")}
  if (base::isFALSE("Sample" %in% colnames(data))){stop("input data must contain column 'Sample'")}
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
  if (base::isFALSE(rank) & base::isFALSE("Count" %in% colnames(data))){
    stop("input data must contain column 'Count', or set rank = TRUE")}
  if (base::isTRUE(rank)){
    regions <- spot.rank(data, regions, include_genes = include_genes)
  }
  pv <- pvalue - 1
  d <- threshold + 1
  i <- 1
  while (pv < pvalue & d > threshold){
    if (i == 2){message(paste((i-1), "Hotspot...", sep = " "))}
    if (i > 2){message(paste((i-1), "Hotspots...", sep = " "))}
    reg1 <- regions[1:i,]
    reg2 <- regions[(i+1):nrow(regions),]
    fun <- function(s, reg_sub){
      reg_list <- c()
      df_sub <- subset(data, Sample == s)
      reg_count <- 0
      for (m in 1:nrow(reg_sub)){
        reg_count <- reg_count + nrow(subset(df_sub, Chromosome == reg_sub$Chromosome[[m]]
                                             & Position %in% reg_sub$Upperbound[[m]]:reg_sub$Lowerbound[[m]]))
      }
      reg_list <- c(reg_list, (reg_count/nrow(reg_sub)))
      return(reg_list)
    }
    reg1_list <- unlist(lapply(unique(data$Sample), fun, reg_sub = reg1))
    reg2_list <- unlist(lapply(unique(data$Sample), fun, reg_sub = reg2))
    ks <- ks.test(reg1_list, reg2_list, alternative = "greater")
    pv <- ks$p.value
    d <- ks$statistic
    i <- i + 1
  }
  sig_spots <- i-2
  ecdf_df <- data.frame("Count" = c(reg1_list, reg2_list),
                        "Group" = c(rep("Hotspot", length(reg1_list)), rep("Non-hotspot", length(reg1_list))))
  regions$number <- 1:nrow(regions)
  regions$type <- c(rep("Hotspot", sig_spots), rep("Non-hotspot", (nrow(regions) - sig_spots)))
  regions$percent <- regions$Count / length(unique(data$Sample)) * 100
  labels <- c()
  if (base::isTRUE(include_genes)){
    for (i in seq_len(nrow(regions))){
      labels <- c(labels, paste(regions$Gene[[i]], " ", regions$Chromosome[[i]], ":", regions$Lowerbound[[i]], "-", regions$Upperbound[[i]], sep = ""))
    }
  }
  if (base::isFALSE(include_genes)){
    for (i in seq_len(nrow(regions))){
      labels <- c(labels, paste(regions$Chromosome[[i]], ":", regions$Lowerbound[[i]], "-", regions$Upperbound[[i]], sep = ""))
    }
  }
  regions$Label <- unlist(labels)
  plots <- plot.sigspot(regions, ecdf_df)
  return_list <- c()
  return_list[[1]] <- regions
  return_list[[2]] <- plots[[1]]
  return_list[[3]] <- plots[[2]]
  return(return_list)

}
