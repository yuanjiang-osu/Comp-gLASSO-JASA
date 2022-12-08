# Comp-gLASSO-JASA

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
