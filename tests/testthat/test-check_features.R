feature <- c("AGE", "SEX", "SMOKING_HISTORY")

data("compSPOT_example_mutations")
data("compSPOT_example_regions")

count_mutations_features <- function(x, regions, data) {
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




test_that("correct number of plots", {
  count.ls <-
    lapply(
      unique(compSPOT_example_mutations$Sample),
      count_mutations_features,
      regions = compSPOT_example_regions,
      data = compSPOT_example_mutations
    )
  count.df <-
    data.frame(
      "Sample" = unique(compSPOT_example_mutations$Sample),
      "Count" = unlist(count.ls)
    )
  clin.df <- compSPOT_example_mutations[, c("Sample", feature)]
  clin.df <- unique(clin.df)
  sample_df <- merge(clin.df, count.df, by = "Sample")
  
  plot_ls <-
    lapply(seq_along(feature),
           plot_features,
           df = sample_df,
           feature = feature)
  expect_equal(length(feature), length(plot_ls))
})
