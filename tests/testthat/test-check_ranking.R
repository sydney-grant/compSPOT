data("compSPOT_example_mutations")
data("compSPOT_example_regions")

shuffled_regions <-
  compSPOT_example_regions[sample(1:nrow(compSPOT_example_regions)),]
ranked_regions <- spot_rank(data = compSPOT_example_mutations,
                            regions = shuffled_regions, include_genes = TRUE)

comp_values <- function(i, reg){
  v1 <- reg$Count[[i]]
  v2 <- reg$Count[[i + 1]]
  test_that("correct ranking", {
    expect_gte(v1, v2)
  })
}

lapply(seq_len(nrow(ranked_regions))[[-length(seq_len(nrow(ranked_regions)))]],
       comp_values, reg = ranked_regions)


