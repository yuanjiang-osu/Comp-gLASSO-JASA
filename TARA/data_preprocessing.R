library(tidyverse)

# Extract pair information from literature
truth <- read_csv("TARA_truth.csv")
truth <- truth %>% rename(interaction_type = 'Ecological interaction type')
truth_filter <- truth %>% filter(!is.na(genus_org1) & !is.na(genus_org2))
unique_pairs <- distinct(truth_filter, genus_org1, genus_org2) %>% arrange(genus_org1, genus_org2)
truth_filter

tara_all <- read_tsv("Database_W5_OTU_occurences.tsv")
taxa_mapping <- read_csv("taxa_mapping.csv")
tara_full <- inner_join(tara_all, taxa_mapping, by = c("cid" = "OTU_id"))
unique_pairs_data = unique_pairs %>% filter(genus_org1 %in% unique(tara_full$genus) & genus_org2 %in% unique(tara_full$genus))
unique_pairs_data

OTU_in_table <- tara_full %>% filter(genus %in% unique(c(unique_pairs_data$genus_org1, unique_pairs_data$genus_org2)))
num_only <- data.frame(OTU_in_table %>% dplyr::select(starts_with("TARA")))
num_gen <- data.frame(OTU_in_table %>% dplyr::select(starts_with("TARA"), genus))

# Aggregate to genus level
## Remove samples of read less than 100
sum(colSums(num_only) > 100)
num_only_sample_filtered <- num_only[, colSums(num_only) > 100]
num_gen_sample_filtered <- num_gen[, colSums(num_only) > 100]

num_gen_sample_filtered_genus_level = num_gen_sample_filtered %>% group_by(genus) %>% summarize_all(sum)
gen_names = num_gen_sample_filtered_genus_level$genus
num_only_sample_filtered_genus_level = num_gen_sample_filtered_genus_level %>% select(-genus)

# Find a reference genus from the OTU's that were not involved in the literature
tara_rest <- tara_full %>% filter(!cid %in% OTU_in_table$cid) # Use OTU's that has taxomic info as reference candidate
tara_rest_num_gen <- data.frame(tara_rest %>% select(starts_with("TARA"), genus))

tara_rest_num_gen_genus_level = tara_rest_num_gen %>% group_by(genus) %>% summarize_all(sum)
tara_rest_num_genus_level = tara_rest_num_gen_genus_level %>% select(-genus)
# Compositions per sample
compo <- prop.table(as.matrix(tara_rest_num_genus_level), 2)
# Average relative abundance per OTU
rel_abd <- apply(compo, 1, mean)
sd_compo <- apply(compo, 1, sd)
stab <- rel_abd / sd_compo
# Sparsity
spar <- apply(tara_rest_num_genus_level == 0, 1, sum) / dim(tara_rest_num_genus_level)[1]

# Select reference OTU based on high stability and low sparsity
#o <- order(spar,decreasing = F, na.last = T)
#o <- order(-stab,decreasing = F, na.last = T)
o <- order(rel_abd, decreasing = T)
# o <- order(sd_compo,decreasing = F, na.last = T)
ref_with_info <- tara_rest_num_gen_genus_level[o,]
ref <- tara_rest_num_genus_level[o,]
which(o == 1)
stab[2]
spar[2]
sd_compo[2]
sd_compo
rel_abd[2]
round(compo[2, ], 5)
# Seems like Acrosphaera is the winner

tara_rest_num_genus_level[1,]
length(tara_rest_num_genus_level[1,])
# Combine column of reference with reference genus
ref_1 <- ref[1,] %>% select(one_of(names(num_only_sample_filtered_genus_level))) # Select the filtered samples
ref_comb <- rbind(num_only_sample_filtered_genus_level, ref_1)
x = t(ref_comb)
saveRDS(x, "TARA_reference.rds")
saveRDS(gen_names, "TARA_genus_names.rds")
saveRDS(unique_pairs_data, "TARA_pairs_literature.rds")
