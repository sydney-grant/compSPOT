## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("compSPOT")

## -----------------------------------------------------------------------------
library(compSPOT)

## ----load mutations-----------------------------------------------------------
data("example_mutations")
head(example_mutations)

## ----load regions-------------------------------------------------------------
data("example_regions")
head(example_regions)

## ----sig.spots----------------------------------------------------------------
significant_spots <- sig.spots(data = example_mutations, regions = example_regions, pvalue = 0.05, threshold = 0.2, include_genes = TRUE, rank = TRUE)

## ----table 1------------------------------------------------------------------
head(significant_spots[[1]])

## ----figure 1a----------------------------------------------------------------
significant_spots[[2]]

## ----figure 1b----------------------------------------------------------------
significant_spots[[3]]

## ----group.spot---------------------------------------------------------------
hotspots <- subset(significant_spots[[1]], type == "Hotspot")

group_comp <- group.spot(data = example_mutations, regions = hotspots, pval = 0.05, threshold = 0.4, 
                         name1 = "High-Risk", name2 = "Low-Risk", include_genes = TRUE)

## ----table 2------------------------------------------------------------------
head(group_comp[[1]])

## ----figure 2a----------------------------------------------------------------
group_comp[[2]]

