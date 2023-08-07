data("example_mutations")
data("example_regions")
spots <- sig.spots(data = example_mutations, regions = example_regions, pvalue = 0.05, threshold = 0.2,
                   include_genes = TRUE, rank = TRUE)


test_that("spots are valid", {
  hotspots <- subset(spots[[1]], type == "Hotspot")
  nonhotspots <- subset(spots[[1]], type == "Non-hotspot")
  h_count <- hotspots$Count / length(unique(data$Sample))
  nh_count <- nonhotspots$Count / length(unique(data$Sample))

  suppressWarnings({ks <- ks.test(h_count, nh_count, alternative = "greater")})
  pv <- ks$p.value
  d <- ks$statistic

  if (pv < pvalue & d > threshold){check = "yes"}
  else{check = "no"}

  expect_equal(check, "yes")
})
