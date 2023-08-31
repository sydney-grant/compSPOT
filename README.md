
## Table of Contents
1. [Introduction](#introduction)
2. [Installation & Setup](#load_package)
3. [Dataset Formatting](#format_data)
  1. [Mutation Data](#mutation_data)
  2. [Genomic Regions](#region_data)
4. [compSPOT Functions](#functions)
  1. [Identifying Mutation Hotspots](#sig.spot)
  2. [Comparison of Mutation Hotspots Between Group](#group.spot)
  3. [Exploratory Data Analysis of Mutated Regions and Personal Risk Factors](#feature.spot)


## Introduction to compSPOT <a name="introduction"/>

Clonal cell groups share common mutations within cancer, precancer, and even 
clinically normal appearing tissues. The frequency and location of these 
mutations may predict prognosis and cancer risk. It has also been well 
established that certain genomic regions have increased sensitivity to acquiring 
mutations. Mutation-sensitive genomic regions may therefore serve as markers 
for predicting cancer risk. This package contains multiple functions to 
establish significantly mutated hotspots, compare hotspot mutation burden 
between samples, and perform exploratory data analysis of the correlation 
between hotspot mutation burden and personal risk factors for cancer, such as 
age, gender, and history of carcinogen exposure. This package allows users to 
identify robust genomic markers to help establish cancer risk.

Currently, minimal resources exist which enable researchers to design their own 
targeted sequencing panels based on specific biological questions and tissues 
of interest. `compSPOT` has been designed to work sequentially with Bioconductor 
package `seq.hotSPOT`. Highly mutated genomic regions identified by `seq.hotSPOT` 
may be used for discovery of significant mutation hotspots with `compSPOT`. 
`compSPOT` may also be used to discover differences in hotspot mutation burden 
between different groups of interest, and the association of mutation burden with 
clinical features. `compSPOT` may be used in combination with the Bioconductor 
package `RTCGA.mutations`, which can be used to pull mutation datasets from the 
TCGA database to be used as input data in various cancer types. Additionally, 
the package `RTCGA.clinical` may be also used to identify highly mutated regions 
in subsets of patients with specific clinical features of interest.


## Installation & Setup <a name="load_package"/>

``` {r install package}
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")                                                         
BiocManager::install("compSPOT")
```

``` {r load library}
library(compSPOT)
```

## Formatting of Input Data <a name = "format_data"/>

### Mutation Data <a name = "mutation_data"/>

The mutation dataset should include the following columns:\
"Chromosome" <-- Chromosome number where the mutation is located\
"Position" <-- Genomic position number where the mutation is located\
"Sample" <-- Unique ID for each sample in dataset\
"Gene" <-- Name of the gene which mutation is located in (optional)\
"Group" <-- Group classification ID (for compare_groups only)\
Clinical Parameters <-- (for compare_features only)


### Genomic Regions <a name = "region_data"/>

The regions dataset should include the following columns:\
"Chromosome" <-- Chromosome number where the region is located\
"Lowerbound" <-- Genomic position number where the region begins\
"Upperbound" <-- Genomic position number where the region ends\
"Gene" <-- Name of the gene which mutation is located in (optional)\
"Count" <-- Number of mutations in mutation dataset which are found within the
region (optional)



## compSPOT Functions <a name = "functions"/>

The compSPOT package contains three main functions for (1) selection of mutation hotspots
(2) comparison of hotspot mutation burden between groups, and (3) comparison of mutation hotspot
enrichment based on clinical and personal risk factors. All functions return both numerical outputs
based on analysis summary and data visualization components for quick and easy interpretation of results.

### Identifying Mutation Hotspots with find_hotspots <a name = "sig.spot"/>

Our previously published Bioconductor package `seq.hotSPOT` 
(doi: 10.3390/cancers15051612) identifies highly mutated genomic regions based 
on SNV datasets. While this tool can identify long lists of mutated regions, 
we sought to establish a method for identifying which of these genomic regions 
have significantly higher mutation frequency compared to others and may be used 
as markers of carcinogenic progression.


Methods: This function begins by measuring the mutation frequency for each 
unique sample for each provided genomic region. Beginning with the top-ranked 
hotspot, a Kolmogorov-Smirnov test is performed on the mutation frequency of 
the top genomic region compared to the normalized mutation frequency of all the 
lower-ranked regions. This continues, then running the Kolmogorov-Smirnov test 
for the normalized mutation frequency of the top 2 genomic regions compared to 
the normalized mutation frequency of all lower-ranked regions. This process 
repeats itself, continuously adding an additional genomic regions each time 
until either the set p-value or empirical distribution threshold is not met. 
Once this cutoff has been reached, an established list of mutation hotspots is 
provided.


### Comparison Mutation Hotspot Burden with compare_groups <a name = "group.spot"/>

Previously, we have shown mutation hotspots identified using seq.hotSPOT may be 
used to differentiate between samples with history of frequent vs infrequent 
carcinogen exposure (doi: 10.3390/cancers15051612, doi: 10.3390/ijms24097852). 
compare_groups provides an automated approach for statistical and visual comparison 
between mutation enrichment of different groups of interest.


Methods: This function creates a list of mutation frequency per unique sample 
for each genomic region separated based on specified sub-groups. The regions 
with significant differences in mutation distribution are calculated using a 
Kolmogorov-Smirnov test. The difference in mutation frequency is output in a 
violin plot.

For this example dataset, the sig.spot function identified 6 hotspots. We will 
use these 6 hotspots to compared the mutation burden between Lung Cancer 
patients with high- and low-risk of disease progression.


### Exploratory Data Analysis of Mutation Hotspot Burden and Personal Risk Factors with compare_features <a name = "feature.spot"/>

Mutation enrichment in cancer mutation hotspots has been shown to relate to 
personal cancer risk factors such as age, gender, and carcinogen exposure 
history and may be used in combination to create predictive models of cancer 
risk (doi: 10.3390/ijms24097852). feature.spot provides a baseline analysis of 
any set of clinical features to identify trends in the enrichment of mutations 
and personal risk factors.

Methods: This function first classifies the features into sequential or 
categorical features. Sequential features are compared to the mutation count 
using Pearson Correlation. Similarly, in categorical features Wilcox Rank Sum 
and Kruska-Wallis Tests are used to compare groups within the features based on 
their mutational count. 




