# Comp-gLASSO-JASA

## Introduction

This paper analyzes both simulated and real data.

The simulated data are generated with random number seed.

The first real dataset is the OTU abundance data downloadable from the Tara Oceans Project data repository (http://taraoceans.sb-roscoff.fr/EukDiv/), including:

Database_W5_OTU_occurences.tsv: OTU abundance data (zipped due to storage limit on GitHub)
taxa_mapping.csv: taxonomic hierarchy for each OTU
TARA_truth.csv: literature validated genus interactions

The code aggregates the OTU abudance data (Database_W5_OTU_occurences.tsv) into genus abundance data based on the taxonomic information (taxa_mapping.csv), estimates the genus interaction networks, and compares with the literature validated genus interactions (TARA_truth.csv).

The second real dataset is the OTU abundance data from the OSU Zebrafish Project, including:

asv.tab: OTU abundance data
tax.tab: taxonomic hierarchy for each OTU
metadata.tab: meta data for samples (zebrafish)

The code aggregates the OTU abudance data (asv.tab) into genus abundance data based on the taxonomic information (tax.tab) and estimates the genus interaction networks separately for uninfected and infected zebrafish based on their parasite burden.

## Instructions

Set the R working directory to the main folder and run the R wrapper file "Wrapper.R" to reproduce all figures and tables in the manuscript. Below are additional details.

SIMULATION

To reproduce Figures 1-2 in Section 3.2, run the following R scripts:

setwd("./Simulation")
source("Simulation_Dense.R")

To reproduce Figures S1-S2 in the supplementary materials, run the following R scripts:

setwd("./Simulation")
source("Simulation_Sparse.R")

REAL DATA/TARA

Run the following R scripts:

setwd("./TARA")
source("data_preprocessing.R")
source("TARA_path.R")

To produce Figure 3(a) in Section 4.1, run the following R scripts:

source("TARA_table_plot_path.R")

To produce Figure 3(b) in Section 4.1 and Figure S3 in the supplementary materials, as well as Table 1 in Section 4.1 and Table S1 in the supplementary materials, run the following R scripts:

source("TARA_StARS.R")
source("TARA_table_plot_StARS.R")

REAL DATA/Zebrafish

Run the following R scripts:

setwd("./Zebrafish")
source("data_preprocessing.R")
source("fish_StARS.R")

To produce Figures 4-5 in Section 4.2 and Table S2 in the supplementary materials, run the following R scripts:

source("fish_table_plot_StARS.R")
