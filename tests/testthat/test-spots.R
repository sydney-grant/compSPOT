data("compSPOT_example_mutations")
data("compSPOT_example_regions")

pval <- 0.05
thres <- 0.2
spots <- sig.spots(data = example_mutations, regions = example_regions,
                   pvalue = pval, threshold = thres,
                   include_genes = TRUE, rank = TRUE)


test_that("spots are valid", {
  hotspots <- subset(spots[[1]], type == "Hotspot")
  nonhotspots <- subset(spots[[1]], type == "Non-hotspot")
  count_h <- function(x){
    count <- 0
    sub <- subset(example_mutations, Sample == x)
    for (i in seq_len(nrow(hotspots))){
      sub2 <- subset(sub, Chromosome == hotspots$Chromosome[[i]] & Position %in%
                       hotspots$Lowerbound[[i]]:hotspots$Upperbound[[i]])
      count <- count + nrow(sub2)
    }
    return(count/nrow(hotspots))
  }
  count_nh <- function(x){
    count <- 0
    sub <- subset(example_mutations, Sample == x)
    for (i in seq_len(nrow(nonhotspots))){
      sub2 <- subset(sub, Chromosome == nonhotspots$Chromosome[[i]] & Position
                     %in%
                       nonhotspots$Lowerbound[[i]]:nonhotspots$Upperbound[[i]])
      count <- count + nrow(sub2)
    }
    return(count/nrow(nonhotspots))
  }
  h_count <- lapply(unique(example_mutations$Sample), count_h)
  nh_count <- lapply(unique(example_mutations$Sample), count_nh)

  suppressWarnings({ks <- ks.test(unlist(h_count), unlist(nh_count),
                                  alternative = "greater")})
  pv <- ks$p.value
  d <- ks$statistic

  if (pv < pval & d > thres){check = "yes"}
  else{check = "no"}

  expect_equal(check, "yes")
})
