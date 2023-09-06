#' count hotspot mutations in each group
#'
#' @description
#' Based on sample name from data, count number of mutations
#' in hotspot of interest based on each group
#'
#'
#' @param sample_name a character of sample name of interest
#' @param data a dataframe containing the chromosome, base pair
#' position, and optionally gene name of each mutation.
#' @param regions a dataframe containing the chromosome, start
#' and end base pair position of each region of interest
#' @param spot an integer containing the row number of hotspot
#' from regions dataframe
#' @param include_genes true or false whether gene names are
#' included in regions dataframe
#'
#' @return A dataframe with the hotspot, group, and mutation count from input
#' sample name
#'
#' @keywords internal
#'
#' @examples
#'
#' count <- count_mutations(sample_name = "CRUK0062-R1", data = data,
#' regions = regions, spot = spot, include_genes = TRUE)
#'
#' @noRd
#'
count_mutations <-
  function(sample_name,
           data,
           regions,
           spot,
           include_genes) {
    sub <- subset(data, Sample == sample_name)
    count <-
      nrow(
        subset(
          sub,
          Chromosome == regions$Chromosome[[spot]] &
            Position %in%
            regions$Lowerbound[[spot]]:regions$Upperbound[[spot]]
        )
      )
    if (include_genes == TRUE) {
      spot_data <- data.frame(
        "Hotspot" = paste(
          regions$Gene[[spot]],
          " ",
          regions$Chromosome[[spot]],
          ":",
          regions$Lowerbound[[spot]],
          "-",
          regions$Upperbound[[spot]],
          sep = ""
        ),
        "Group" = unique(sub$Group),
        "Mutation_Count" = count
      )
    }
    if (include_genes == FALSE) {
      spot_data <-
        data.frame(
          "Hotspot" = paste(
            regions$Chromosome[[spot]],
            ":",
            regions$Lowerbound[[spot]],
            "-",
            regions$Upperbound[[spot]],
            sep = ""
          ),
          "Group" = unique(sub$Group),
          "Mutation_Count" = count
        )
    }
    return(spot_data)
  }

#' find differences in mutatino burden per hotspot in different groups
#'
#' @description
#' Based on sample name from data, count number of mutations
#' in hotspot of interest based on each group
#'
#' @param data a dataframe containing the chromosome, base pair
#' position, and optionally gene name of each mutation.
#' @param return a value of either "plot" or "hotspots" based on which type of
#' data will be returned
#' @param regions a dataframe containing the chromosome, start
#' and end base pair position of each region of interest
#' @param spot an integer containing the row number of hotspot
#' from regions dataframe
#' @param include_genes true or false whether gene names are
#' included in regions dataframe
#' @param pvalue a threshold p-value for Kolmogorov-Smirnov test
#' @param threshold the cutoff empirical distribution for
#' Kolmogorov-Smirnov test
#'
#' @return A dataframe with either the mutation count for each group in each
#' hotspots, or a dataframe of the distance and p-value of hotspots with
#' significant differences between the two groups.
#'
#' @keywords internal
#'
#' @examples
#'
#' group <- count_groups(spot = 1, regions = regions, name1 = "High-Risk",
#' name2 = "Low-Risk", include_genes = TRUE, return = "plot")
#'
#' @noRd
#'
count_groups <-
  function(spot,
           data,
           regions,
           name1,
           name2,
           include_genes,
           return,
           pvalue,
           threshold) {
    plot_data <- data.table::rbindlist(
      lapply(
        unique(data$Sample),
        count_mutations,
        data = data,
        regions = regions,
        spot = spot,
        include_genes =
          include_genes
      )
    )
    count_1 <- subset(plot_data, Group == name1)$Mutation_Count
    count_2 <- subset(plot_data, Group == name2)$Mutation_Count
    if (sum(c(count_1, count_2)) > 0) {
      ks <- ks.test(count_1, count_2, alternative = "two.sided")
      if (ks$p.value < pvalue & ks$statistic >= threshold) {
        if (sum(count_1) > sum(count_2)) {
          high_group <- "Group 1"
        }
        else{
          high_group <- "Group 2"
        }
        if (include_genes == TRUE) {
          diff_hotspots <- data.frame(
            "D" = ks$statistic,
            "pval" = ks$p.value,
            "Higher Group" = high_group,
            "Hotspot" = paste(
              regions$Gene[[spot]],
              " ",
              regions$Chromosome[[spot]],
              ":",
              regions$Lowerbound[[spot]],
              "-",
              regions$Upperbound[[spot]],
              sep = ""
            )
          )
        }
        if (include_genes == FALSE) {
          diff_hotspots <-
            data.frame(
              "D" = ks$statistic,
              "pval" = ks$p.value,
              "Higher Group" = high_group,
              "Hotspot" = paste(
                regions$Chromosome[[spot]],
                ":",
                regions$Lowerbound[[spot]],
                "-",
                regions$Upperbound[[spot]],
                sep = ""
              )
            )
        }
      }
      if (ks$p.value >= pvalue |
          ks$statistic < threshold) {
        diff_hotspots <- data.frame()
      }
    }
    if (sum(c(count_1, count_2)) == 0) {
      diff_hotspots <- data.frame()
    }
    if (return == "plot") {
      return(plot_data)
    }
    if (return == "hotspots") {
      return(diff_hotspots)
    }
  }


#' hotspot group comparison with plot
#'
#' @description
#' Based on a dataframe of hotspots, this function groups the
#' hotspots into two different groups and overlays interactive
#' violin plots based on mutation count for comparison.
#'
#'
#' @importFrom ggplot2 ggplot aes stat_ecdf theme theme_classic
#' scale_color_manual element_text geom_point geom_vline xlab ylab
#' geom_violin scale_fill_manual geom_smooth position_dodge geom_dotplot
#' @importFrom plotly ggplotly plot_ly add_trace layout
#' @importFrom magrittr %>%
#' @importFrom gridExtra grid.arrange
#' @param plot_data a dataframe containing the hotspot,
#' mutation count, and associated groups
#' @param group1 a string containing the name of one group for
#' the comparison
#' @param group2 a string containing the name of the second group
#' for the comparison
#'
#' @return A plotly object violin plot displaying the comparison of
#' mutation counts for each group.
#'
#' @keywords internal
#'
#' @examples
#'
#' data("plot_data")
#' fig <- plot_spot(plot_data = plot_data, group1 = "High-Risk",
#' group2 = "Low-Risk")
#'
#' @noRd


plot_spot <- function(plot_data, group1, group2) {
  plot1 <-
    ggplot2::ggplot(data = plot_data,
                    ggplot2::aes(y = Mutation_Count,
                                 x = Hotspot,
                                 fill = Group)) +
    ggplot2::geom_violin(alpha = 0.75,
                         trim = FALSE,
                         position = position_dodge()) +
    ggplot2::scale_fill_manual(values = c("brown4", "darksalmon")) +
    ggplot2::scale_color_manual(values = c("brown4", "darksalmon")) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        vjust = 0.5,
        hjust = 0.5,
        size = 10
      ),
      axis.title.y = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 4),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12)
    )
  plot1 <- plotly::ggplotly(plot1)
  return(plot1)
}




#' ECDF plot of hotspot mutation comparison
#'
#' @description
#' Based on a dataframe of hotspots, this function creates a ECDF plot to
#' compare the difference in hotspot mutation burden between two groups.
#'
#'
#' @importFrom ggplot2 ggplot aes stat_ecdf theme theme_classic
#' scale_color_manual element_text geom_point geom_vline xlab ylab
#' geom_violin scale_fill_manual geom_smooth
#' @param spot a string containing the name of a unique hotspot in present in
#' plot_data
#' @param plot_data a dataframe containing the hotspot,
#' mutation count, and associated groups
#'
#' @return A ggplot object ECDF plot displaying the comparison of
#' mutation counts for each group.
#'
#' @keywords internal
#'
#' @examples
#'
#' data("plot_data")
#' fig <- plot_ecdf(spot = "TP53 17:7578429-7578528",
#' plot_data = plot_data,
#' group1 = "High-Risk",
#' group2 = "Low-Risk")
#'
#' @noRd


plot_ecdf <- function(spot, plot_data) {
  plot_data.sub <- subset(plot_data, Hotspot == spot)
  plot2 <- ggplot2::ggplot(data = plot_data.sub,
                           ggplot2::aes(x = Mutation_Count, group = Group,
                                        col = Group)) +
    ggplot2::stat_ecdf(geom = "step") +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = c("brown4", "darksalmon")) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 0.5,
        size = 10
      ),
      axis.title.y = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 4),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12)
    )
  return(plot2)
}

#' hotspot comparison by group
#'
#' @description
#' This function compares the mutation frequency of a panel of
#' genomic regions between two sub-groups.
#'
#'
#' @details
#' This function creates a list of mutation frequency per unique
#' sample for each genomic regions separated based on specified
#' sub-groups. The regions with significant differences in mutation
#' distribution are calculated using a Kolmogorov-Smirnov test. The
#' difference in mutation frequency is output in a violin plot.
#'
#' @importFrom stats ks.test
#' @importFrom magrittr %>%
#' @importFrom data.table rbindlist
#' @param data a dataframe containing the chromosome, base pair
#' position, and optionally gene name of each mutation. Dataframe must contain
#' columns with the following names:
#' "Chromosome" <-- Chromosome number where the mutation is located
#' "Position" <-- Genomic position number where the mutation is located
#' "Sample" <-- Unique ID for each sample in dataset
#' "Gene" <-- Name of the gene which mutation is located in (optional)
#' "Group" <-- Group classification ID (group.spot only)
#' @param regions a dataframe containing the chromosome, start and
#' end base pair position of each region of interest
#' @param pvalue a threshold p-value for Kolmogorov-Smirnov test
#' @param threshold the cutoff empirical distribution for
#' Kolmogorov-Smirnov test
#' @param name1 a string containing the name of one group for the
#' comparison
#' @param name2 a string containing the name of the second group for
#' the comparison
#' @param include_genes true or false whether gene names are included
#' in regions dataframe
#'
#' @return a list containing the following:
#' 1. A dataframe with the hotspot, group, and mutation count from input
#' sample name
#' 2. A plotly object violin plot comparing the mutation frequency per sample in groups as
#' given by "name1" and "name2" variables
#' 3. An array of ECDF plots comparing the mutation frequency per sample in
#' groups as given by "name1" and "name2" variables
#'
#' @examples
#'
#' data("compSPOT_example_mutations")
#' data("compSPOT_example_regions")
#' compare_groups(data = compSPOT_example_mutations,
#' regions = compSPOT_example_regions, pvalue = 0.05, threshold = 0.4,
#' name1 = "High-Risk", name2 = "Low-Risk", include_genes = TRUE)
#'
#' @export
#'


compare_groups <-
  function(data,
           regions,
           pvalue,
           threshold,
           name1,
           name2,
           include_genes) {
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
    if (base::isFALSE("Group" %in% colnames(data))) {
      stop("input data must contain column 'Group'")
    }
    if (base::isFALSE(length(unique(data$Group)) == 2)) {
      stop("input data must have two unique group labels")
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
    plot_data <-
      data.table::rbindlist(
        lapply(
          seq_len(nrow(regions)),
          count_groups,
          data = data,
          regions = regions,
          name1 = "High-Risk",
          name2 = "Low-Risk",
          include_genes = include_genes,
          return = "plot",
          pvalue = pvalue,
          threshold = threshold
        )
      )
    all_diff <- data.table::rbindlist(
      lapply(
        seq_len(nrow(regions)),
        count_groups,
        data = data,
        regions = regions,
        name1 = "High-Risk",
        name2 = "Low-Risk",
        include_genes = include_genes,
        return = "hotspots",
        pvalue = pvalue,
        threshold = threshold
      )
    )
    plot2_list <- lapply(unique(plot_data$Hotspot), plot_ecdf,
                         plot_data = plot_data)
    plot <- plot_spot(plot_data, group1 = name1,
                      group2 = name2)
    return_items <-
      list(all_diff,
           plot,
           ggpubr::ggarrange(
             plotlist = plot2_list,
             nrow = 2,
             ncol = ceiling(length(plot2_list)
                            / 2)
           ))
    return(return_items)
  }
