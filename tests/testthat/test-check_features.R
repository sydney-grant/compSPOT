feature <- c("AGE", "SEX", "SMOKING_HISTORY",
              "TUMOR_VOLUME", "KI_67")

data("compSPOT_example_mutations")
data("compSPOT_example_regions")


count.ls <- lapply(unique(data$Sample), count_mutations, regions = compSPOT_example_regions, data = compSPOT_example_mutations)
count.df <-
  data.frame("Sample" = unique(compSPOT_example_mutations$Sample),
             "Count" = unlist(count.ls))
clin.df <- data[, c("Sample", feature)]
clin.df <- unique(clin.df)
sample_df <- merge(clin.df, count.df, by = "Sample")

plot_ls <- lapply(seq_along(feature), plot_features, df = sample_df, feature = feature)


test_that("correct number of plots", {
  expect_equal(length(feature), length(plot_ls))
})
