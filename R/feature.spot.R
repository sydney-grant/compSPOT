#' hotspot comparison by additional features
#'
#' @description
#' This function classifying the features into categorical and sequential features
#' and return a grid of violin plots for categorical features and scatter plot for sequential features.
#' @details
#' This function is used to classify the features into sequential features if values are integers or classifies them into categorical #' features.
#' Sequential features are compared to the mutation count using Pearson correlation.
#' Similarly, in categorical features wilcox test is used to compare between the groups in the features based on their mutational
#' count.
#' Scatter plot is used to represent the sequential features along with the R and p-value from the pearson correlation.
#' Similarly violin plots are used to plot the groups in the categorical data and wilcox or kruska-wallis values are shown on the
#' graph.
#'
#'
#' @param data A dataframe containing the clinical features and the mutation count.
#' @param features A dataframe containing all the features.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom plotly ggplotly
#' @importFrom ggpubr stat_cor
#' @importFrom ggpubr stat_compare_means
#' @importFrom ggpubr ggarrange
#' @importFrom utils combn
#' @importFrom stats wilcox.test
#' @importFrom gridExtra grid.arrange
#'
#'
#'
#' @return A grid of all the violin plots for the categorical data and scatter plot for the sequential data.
#'
#'
#' @examples
#'
#' data("example_mutations")
#' data("example_regions")
#' features <- c("AGE", "SEX", "ADJUVANT_TX", "SMOKING_HISTORY", "TUMOR_VOLUME", "KI_67")
#' feature.spot(data = example_mutations, regions = example_regions, feature = features)
#'
#' @export
#'

feature.spot <-function(data, regions, feature){
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
  count.ls <- c()
  for (sample in unique(data$Sample)){
    sub <- subset(data, Sample == sample)
    count <- 0
    for (i in seq_len(nrow(regions))){
      sub <- subset(sub, Chromosome == regions$Chromosome[[i]] & Position %in% regions$Lowerbound[[i]]:regions$Upperbound[[i]])
      count <- count + nrow(sub)
    }
    count.ls <- c(count.ls, count)
  }
  count.df <- data.frame("Sample" = unique(data$Sample), "Count" = unlist(count.ls))
  clin.df <- data[,c("Sample", feature)]
  clin.df <- unique(clin.df)
  sample_df <- merge(clin.df, count.df, by = "Sample")
  feature_list <- feature
  plot_ls <- list()
  Sequential <-function(f){
    result <- eval(parse(text=paste0("sample_df$",f)))
    plot <-ggplot2::ggplot(data=sample_df,aes(x=result,y=Count))+
      ggplot2::geom_point(color = "brown4", fill = "brown4")+
      ggpubr::stat_cor()+
      ggplot2::geom_smooth(method = 'lm', color = "darksalmon")+
      ggplot2::theme_classic() +
      ggplot2::xlab(f) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 12),
                     axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title = element_text(size = 4),
                     axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 12))
    return(plot)
  }
  Categorical <-function(f)
  {
    result <- eval(parse(text=paste0("sample_df$",f)))
    cat_list <-unique(result)
    comb <- combn(cat_list,2)
    for (cat in 1:(length(comb)-1)){
      cat1 <- comb[cat]
      cat2 <- comb[cat+1]
      message(paste("Wilcoxon test for Feature",f,"-", cat1,"vs", cat2))
      subset1 <- subset(sample_df,result==cat1)$Count
      subset2 <- subset(sample_df,result==cat2)$Count
      wilcox_result <- wilcox.test(subset1,subset2)
      print(wilcox_result)
    }
    plotc <- ggplot2::ggplot(sample_df,aes(x=result,y=Count, fill =result))+
      ggplot2::geom_violin(trim = FALSE)+
      ggpubr::stat_compare_means( label.y = 1.15)+
      ggplot2::scale_fill_manual(values = c("brown4", "darksalmon", "tan", "darkorange", "white", "grey")) +
      ggplot2::scale_color_manual(values = c("brown4", "darksalmon", "tan", "darkorange", "white", "grey")) +
      ggplot2::theme_classic() +
      ggplot2::xlab(f) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
                     axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title = element_text(size = 4),
                     axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 12),
                     legend.position = "none")

    return(plotc)
  }
  for (i in seq_len(length(feature_list))){
    f <- feature_list[[i]]

    if(sapply(sample_df[f],class)=="numeric" | sapply(sample_df[f],class)=="integer")
    {
      plot_s <- Sequential(f)
      plot_ls[[i]] <- plot_s
    }
    else
    {
      plot_c <- Categorical(f)
      plot_ls[[i]] <- plot_c
    }
  }
  plots <- ggpubr::ggarrange(plotlist=plot_ls,nrow=2, ncol = ceiling(length(plot_ls)/2))
  return(plots)

}
