#' example genomic regions
#' @name example_regions
#' @docType data
#' @aliases example_regions
#' @title Genomic Coordinates of Regions of Interest
#' @usage example_regions
#' @keywords datasets
#'
#' @return A dataframe containing the chromosome number, base pair location, and gene names of 200 genomic regions 
#' highly mutated in Lung Squamous Cell Carcinoma identified using seq.hotSPOT
#'
#' @format ## `example_regions`
#' a dataframe with 2 columns and 200 rows:
#' \describe{
#'  \item{Lowerbound}{Base pair position of the start of the region}
#'  \item{Upperbound}{Base pair position of the end of the region}
#'  \item{Chromosome}{Chromosome which the mutation is located on}
#'  \item{Gene}{Name of gene affected by mutation}
#'  
#' }
#'
#' @description A dataframe containing the chromosome number, lowerbound and upperbound base pair locations of each region of 
#' interest along with the name of the gene where the region is located. Each row indicates a unique region. Regions were
#' identified using the seq.hotSPOT package based on Lung Squamous Cell Carcinoma highly mutated regions.
#'
#' @references
#' Grant SR et al;
#' HotSPOT: A Computational Tool to Design Targeted Sequencing Panels to Assess Early Photocarcinogenesis. 
#' Cancers (Basel). 2023 Mar 5;15(5):1612. doi: 10.3390/cancers15051612. PMID: 36900402; PMCID: PMC10001346.
#' 
#' Grant S, Wei L, Paragh G (2023). seq.hotSPOT: Targeted sequencing panel design based on mutation hotspots. 
#' R package version 1.0.0, https://github.com/sydney-grant/seq.hotSPOT.
#'
NULL