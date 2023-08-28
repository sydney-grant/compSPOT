feature <- c("AGE", "SEX", "ADJUVANT_TX", "SMOKING_HISTORY",
              "TUMOR_VOLUME", "KI_67")

plot_ls <- lapply(seq_along(feature), plot_features, df = sample_df, feature = feature)


test_that("correct number of plots", {
  expect_equal(length(feature), length(plot_ls))
})
