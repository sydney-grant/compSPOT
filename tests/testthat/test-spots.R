data("compSPOT_example_mutations")
data("compSPOT_example_regions")

pval <- 0.05
thres <- 0.2

test_that("spots are valid", {
  spots <-
    find_hotspots(
      data = compSPOT_example_mutations,
      regions = compSPOT_example_regions,
      pvalue = pval,
      threshold = thres,
      include_genes = TRUE,
      rank = TRUE
    )
  hotspots <- subset(spots[[1]], type == "Hotspot")
  nonhotspots <- subset(spots[[1]], type == "Non-hotspot")
  cutoff <- hotspots$percent[[nrow(hotspots)]]
  if (max(nonhotspots$percent) < min(hotspots$percent)) {
    check = TRUE
  }
  else{
    check = FALSE
  }
  
  expect_true(check)
})
