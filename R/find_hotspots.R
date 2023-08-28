#' genomic region ranker
#'
#' @description
#' ranks all genomic regions from highest to lowest
#' mutation frequency in given dataset
#'
#' @param data A dataframe containing the location of
#' each mutation.
#' @param regions a dataframe containing the location
#' of gene region of interest
#' @param include_genes true or false whether gene names
#' are included in regions dataframe
#'
#' @keywords internal
#'
#' @return A dataframe containing the ranked genomic
#' regions from highest to lowest mutation frequency
#'
#'
#' @examples
#'
#' data("mutation_data")
#' spot_rank(data = mutation_data, regions = ea_regions,
#' include_genes = TRUE)
#'
#'@noRd

spot_rank <- function(data, regions, include_genes) {
  ranked_regions <- data.frame()
  for (i in seq_len(nrow(regions))) {
    data_sub <- subset(
      data,
      Chromosome == regions$Chromosome[[i]] &
        Position %in% regions$Lowerbound[[i]]:regions$Upperbound[[i]]
    )
    if (include_genes == TRUE) {
      ranked_row <-
        data.frame(
          "Chromosome" = regions$Chromosome[[i]],
          "Lowerbound" = regions$Lowerbound[[i]],
          "Upperbound" = regions$Upperbound[[i]],
          "Count" = nrow(data_sub),
          "Gene" = regions$Gene[[i]]
        )
    }
    if (include_genes == FALSE) {
      ranked_row <-
        data.frame(
          "Chromosome" = regions$Chromosome[[i]],
          "Lowerbound" = regions$Lowerbound[[i]],
          "Upperbound" = regions$Upperbound[[i]],
          "Count" = nrow(data_sub)
        )
    }
    ranked_regions <- rbind(ranked_regions, ranked_row)
  }
  ranked_regions <- ranked_regions[order(-ranked_regions$Count),]
  return(ranked_regions)
}

#' create labels for final hotspot dataframe
#'
#' @description
#' creates a unique label for each hotspot based on the gene name and genomic
#' location.
#'
#' @param regions a dataframe containing the location
#' of gene region of interest
#' @param include_genes true or false whether gene names
#' are included in regions dataframe
#' #' @param i a numeric value of which row from the region dataframe to make a
#' label for
#'
#' @keywords internal
#'
#' @return A string containing the unique label for hotspot
#'
#'
#' @examples
#'
#' create_labels(i = 1, regions = compSPOT_example_regions,
#' include_genes = TRUE)
#'
#'@noRd
#'


create_labels <- function(i, regions, include_genes) {
  if (base::isTRUE(include_genes)) {
    label <- paste(
      regions$Gene[[i]],
      " ",
      regions$Chromosome[[i]],
      ":",
      regions$Lowerbound[[i]],
      "-",
      regions$Upperbound[[i]],
      sep = ""
    )
  }
  if (base::isFALSE(include_genes)) {
    label <- paste(regions$Chromosome[[i]],
                   ":",
                   regions$Lowerbound[[i]],
                   "-",
                   regions$Upperbound[[i]],
                   sep = "")
  }
  return(label)
}

#' counting mutations per regions
#'
#' @description
#' creates a unique label for each hotspot based on the gene name and genomic
#' location.
#'
#' @param reg_sub subset dataframe of total regions dataframe which contains
#' either top or bottom genomic regions for comparison
#' @param s a string containing ID for unique sample in mutation dataset
#' @param data a dataframe containing the chromosome, base pair position,
#' and optionally gene name of each mutation
#' @keywords internal
#'
#' @return A list of mutation frequencies per genomic region.
#'
#'
#' @examples
#'
#' region_counter(s = "CRUK0062-R1", reg_sub = reg1)
#'
#'@noRd
#'

region_counter <- function(s, reg_sub, data) {
  reg_list <- c()
  df_sub <- subset(data, Sample == s)
  reg_count <- 0
  for (m in seq_len(nrow(reg_sub))) {
    reg_count <-
      reg_count + nrow(
        subset(
          df_sub,
          Chromosome == reg_sub$Chromosome[[m]]
          &
            Position %in% reg_sub$Upperbound[[m]]:reg_sub$Lowerbound[[m]]
        )
      )
  }
  reg_list <- c(reg_list, (reg_count / nrow(reg_sub)))
  return(reg_list)
}





#' sig.spot visualization
#'
#' @description
#' Based on a panel of genomic regions, this function visualizes the
#' regions which are found to have significance based on percentage
#' of samples with hotspot mutations.
#'
#' @details
#' This function begins by displaying a plot with a ranked number of
#' hotspot based on number of mutations on the x-axis versus percentage
#' of all samples with mutations in hotspot on the y-axis. Regions of high
#' significance are ranked on the left side of the plot. Based on the
#' type, the points plotted as red dots are the hotspots which were found
#' to have a significant number of mutations whereas the pink were
#' found to not be significant. An interactive label of each point is also
#' available. A second plot is also created using the same data used in the
#' Kolmogorov-Smirnov calculation to show a visualization of the ECDF
#' of hotspots vs non-hotspots.
#'
#' @importFrom ggplot2 ggplot aes stat_ecdf theme theme_classic
#' scale_color_manual element_text geom_point geom_vline xlab ylab
#' @importFrom plotly ggplotly
#' @param df1 a dataframe containing the labeled region number,
#' type of region, percentage of samples with
#' hotspot mutation
#' @param df2 a dataframe containing the data required for calculation
#' of ECDF
#'
#' @return A figure organizing the region number in ranked order versus the
#' percentage of the sample with hotspot mutation which also denotes the
#' label of the data point is displayed when hovering over the item. A second
#' figure providing a visualization of the ECDF for hotspot vs non-hotspots.
#'
#' @keywords internal
#' @examples
#'
#' data("regions")
#' data("ecdf_df")
#' plot <- plot_sigspot(df1 = region, df2 = ecdf_df)
#'
#'@noRd


plot_sigspot <- function(df1, df2) {
  plot1 <-
    ggplot2::ggplot(df1, ggplot2::aes(
      y = percent,
      x = number,
      text = paste("Label", Label)
    )) +
    ggplot2::geom_point(
      data = subset(df1, type == "Hotspot"),
      color = "brown4",
      fill = "brown4",
      alpha = 1,
      shape = 21,
      size = 1.5
    ) +
    ggplot2::geom_point(
      data = subset(df1, type == "Non-hotspot"),
      color = "darksalmon",
      fill = "darksalmon",
      alpha = 1,
      shape = 21,
      size = 1.5
    ) +
    ggplot2::geom_vline(
      xintercept = (nrow(subset(
        df1, type == "Hotspot"
      )) + 1),
      linetype = "longdash",
      linewidth = 0.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Percentage of Samples \n with Hotspot Mutation") +
    ggplot2::xlab("Ranked Region Number") +
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
  plot1 <- plotly::ggplotly(plot1, tooltip =  c("text"))
  plot2 <-
    ggplot2::ggplot(data = df2, ggplot2::aes(x = Count, group = Group,
                                             col = Group)) +
    ggplot2::stat_ecdf(geom = "step") +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = c("brown4", "darksalmon"))
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      vjust = 0.5,
      hjust = 0.5,
      size = 12
    ),
    axis.title.y = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 4),
    axis.text.y = ggplot2::element_text(size = 12),
    plot.title = ggplot2::element_text(hjust = 0.5, size = 12)
  )
  plot2 <- plotly::ggplotly(plot2, tooltip =  c("text"))
  output <- c()
  output[[1]] <- plot1
  output[[2]] <- plot2
  return(output)
}


#' significant hotspot calculator
#'
#' @description
#' Based on a panel of genomic regions, this function calculates the
#' regions which are found to have significantly higher mutation
#' frequency compared to less mutated regions.
#'
#' @details
#' This function begins by measuring the mutation frequency for each
#' unique sample for each provided genomic region. Beginning with the
#' top ranked hotspot, a Kolmogorov-Smirnov test is preformed on the
#' mutation frequency of the top genomic region compared to the normalized
#' mutation frequency of all the lower-ranked regions. This continues,
#' then running the Kolmogorov-Smirnov test for the normalized mutation
#' frequency of the top 2 genomic regions compared to the normalized
#' mutation frequency of all lower-ranked regions.This process repeats itself,
#' continuously adding an additional genomic regions each time until either
#' the set p-value or empirical distribution threshold is not met. Once this
#' cutoff has been reached, an established list of mutation hotspots is
#' provided.
#'
#' @importFrom stats ks.test
#' @param data a dataframe containing the chromosome, base pair
#' position, and optionally gene name of each mutation. Dataframe must contain
#' columns with the following names:
#' "Chromosome" <-- Chromosome number where the mutation is located
#' "Position" <-- Genomic position number where the mutation is located
#' "Sample" <-- Unique ID for each sample in dataset
#' "Gene" <-- Name of the gene which mutation is located in (optional)
#' @param regions a dataframe containing the chromosome, start and end
#' base pair position of each region of interest
#' @param pvalue the p-value cutoff for included hotspots
#' @param threshold the cutoff empirical distribution for Kolmogorov-Smirnov
#' test
#' @param include_genes true or false whether gene names are included in
#' regions dataframe
#' @param rank true or false whether regions dataframe is already ranked and
#' includes mutation count of total dataset
#'
#' @return A list containing the following:
#' 1. dataframe containing the genomic regions with significant
#' mutation frequency
#' 2. Dotplot showing the percentage of samples with mutations in each ranked
#' genomic region, highlighting significantly mutated hotspots
#' 3. ECDF plot showing the difference in mutation frequency between hotspots
#' and non-hotspots
#'
#'
#' @examples
#'
#' data("compSPOT_example_mutations")
#' data("compSPOT_example_regions")
#' significant_spots <- find_hotspots(data = compSPOT_example_mutations,
#' regions = compSPOT_example_regions,
#' pvalue = 0.05, threshold = 0.2, include_genes = TRUE, rank = TRUE)
#'
#'
#' @export
#'


find_hotspots <-
  function(data,
           regions,
           pvalue,
           threshold,
           include_genes,
           rank) {
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
    if (base::isTRUE(include_genes) &
        base::isFALSE("Gene" %in% colnames(data)))
    {
      stop("input data must contain column 'Gene'")
    }
    if (base::isTRUE(include_genes) &
        base::isFALSE(is.character(data$Gene))) {
      stop("column 'Gene' must be in the form of character")
    }
    if (pvalue > 1) {
      stop("pvalue cannot be greater than 1")
    }
    if (threshold > 1) {
      stop("threshold cannot be greater than 1")
    }
    if (pvalue < 0) {
      stop("pvalue cannot be less than 0")
    }
    if (threshold < 0) {
      stop("threshold cannot be less than 0")
    }
    if (base::isFALSE(rank) &
        base::isFALSE("Count" %in% colnames(data))) {
      stop("input data must contain column 'Count', or set rank = TRUE")
    }
    if (base::isTRUE(rank)) {
      regions <- spot_rank(data, regions, include_genes = include_genes)
    }
    pv <- pvalue - 1
    d <- threshold + 1
    i <- 1
    while (pv < pvalue & d > threshold) {
      reg1 <- regions[seq_len(i),]
      reg2 <- regions[(i + 1):nrow(regions),]
      reg1_list <- unlist(lapply(
        unique(data$Sample),
        region_counter,
        reg_sub = reg1,
        data = data
      ))
      reg2_list <- unlist(lapply(
        unique(data$Sample),
        region_counter,
        reg_sub = reg2,
        data = data
      ))
      ks <- ks.test(reg1_list, reg2_list, alternative = "greater")
      pv <- ks$p.value
      d <- ks$statistic
      i <- i + 1
    }
    sig_spots <- i - 2
    ecdf_df <- data.frame("Count" = c(reg1_list, reg2_list),
                          "Group" = c(rep("Hotspot", length(reg1_list)),
                                      rep("Non-hotspot", length(reg1_list))))
    regions$number <- seq_len(nrow(regions))
    regions$percent <-
      regions$Count / length(unique(data$Sample)) * 100
    perc_cutoff <- regions$percent[[sig_spots]]
    regions$type <-
      c(rep("Hotspot", nrow(subset(
        regions, percent >= perc_cutoff
      ))),
      rep("Non-hotspot", nrow(subset(
        regions, percent < perc_cutoff
      ))))

    labels <-
      lapply(seq_len(nrow(regions)),
             create_labels,
             regions = regions,
             include_genes = include_genes)
    regions$Label <- unlist(labels)
    plots <- plot_sigspot(regions, ecdf_df)
    return_list <- c()
    return_list[[1]] <- regions
    return_list[[2]] <- plots[[1]]
    return_list[[3]] <- plots[[2]]
    return(return_list)
  }
