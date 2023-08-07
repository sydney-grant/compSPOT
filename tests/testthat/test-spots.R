data("example_mutations")
data("example_regions")

pval <- 0.05
thres <- 0.2
spots <- sig.spots(data = example_mutations, regions = example_regions, pvalue = pval, threshold = thres,
                   include_genes = TRUE, rank = TRUE)


test_that("spots are valid", {
  hotspots <- subset(spots[[1]], type == "Hotspot")
  nonhotspots <- subset(spots[[1]], type == "Non-hotspot")
  h_count <- hotspots$Count / length(unique(example_mutations$Sample))
  nh_count <- nonhotspots$Count / length(unique(example_mutations$Sample))

  suppressWarnings({ks <- ks.test(h_count, nh_count, alternative = "greater")})
  pv <- ks$p.value
  d <- ks$statistic

  if (pv < pval & d > thres){check = "yes"}
  else{check = "no"}

  expect_equal(check, "yes")
})
