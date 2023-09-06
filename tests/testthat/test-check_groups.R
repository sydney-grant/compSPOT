data("compSPOT_example_mutations")
data("compSPOT_example_regions")


test_that("check number significant spots", {
    significant_spots <-
        find_hotspots(
            data = compSPOT_example_mutations,
            regions = compSPOT_example_regions,
            pvalue = 0.05,
            threshold = 0.2,
            include_genes = TRUE,
            rank = TRUE
        )
    hotspots <- subset(significant_spots[[1]], type == "Hotspot")
    comp_test <- compare_groups(
        data = compSPOT_example_mutations,
        regions = hotspots,
        pvalue = 0.05,
        threshold = 0.4,
        name1 = "High-Risk",
        name2 = "Low-Risk",
        include_genes = TRUE
    )
    expect_true(nrow(comp_test[[1]]) <= length(unique(hotspots$Label)))
})
