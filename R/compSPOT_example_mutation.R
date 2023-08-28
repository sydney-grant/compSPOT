#' compSPOT example mutation data
#' @name compSPOT_example_mutations
#' @docType data
#' @aliases compSPOT_example_mutations
#' @title Single Nucleotide Variants and Patient Features in Lung Cancer
#' Patients
#' @usage compSPOT_example_mutations
#' @keywords datasets
#'
#' @return A dataframe containing the chromosome number, base pair location,
#' sample ID, gene name, patient features including:
#' age, sex, adjuvant therapy treatment, smoking history, tumor volume, ki67
#' quantification, and risk-classification for
#' cancer progression.
#'
#' @format ## `example_mutations`
#' a dataframe with 11 columns and 22947 rows:
#' \describe{
#'  \item{Sample}{ID assigned to indicate each unique sample}
#'  \item{Gene}{Name of gene affected by mutation}
#'  \item{Chromosome}{Chromosome which the mutation is located on}
#'  \item{Position}{Base pair position of mutation}
#'  \item{AGE}{Age of the patient}
#'  \item{SEX}{Sex of the patient}
#'  \item{ADJUVANT_TX}{Statement of whether or not patient recieved adjuvant
#'  therapy}
#'  \item{SMOKING_HISTORY}{Patient's history of smoking}
#'  \item{TUMOR_VOLUME}{Measured volumne of patient's lung tumor}
#'  \item{KI_67}{Quantification of ki67 markers observed in each patient}
#'  \item{Group}{Risk classification of patients based on observed survival}
#' }
#'
#' @description A dataframe containing the chromosome number, base pair
#' location, sample ID, gene name,
#' patient features including: age, sex, adjuvant therapy treatment, smoking
#' history, tumor volume,
#' ki67 quantification, and risk-classification for cancer progression. Data
#' curated from cBioPortal dataset:
#' Non-Small Cell Lung Cancer (TRACERx, NEJM & Nature 2017)
#'
#' @references
#' Abbosh C et al.;
#' TRACERx consortium; PEACE consortium;
#' Swanton C. Phylogenetic ctDNA analysis depicts early-stage lung cancer
#' evolution.
#' Nature. 2017 Apr 26;545(7655):446-451. doi: 10.1038/nature22364. Erratum in:
#' Nature. 2017 Dec 20;:
#' PMID: 28445469; PMCID: PMC5812436.
#'
#' Jamal-Hanjani M et al.;
#' TRACERx Consortium. Tracking the Evolution of Non-Small-Cell Lung Cancer.
#' N Engl J Med. 2017 Jun 1;376(22):2109-2121. doi: 10.1056/NEJMoa1616288.
#' Epub 2017 Apr 26. PMID: 28445112.
#'
#'
NULL
