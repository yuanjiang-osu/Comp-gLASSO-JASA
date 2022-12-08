rm(list = ls())

library(tidyverse)

asv.data <- read.table("asv.tab", header = TRUE)
names(asv.data) <- gsub("[.]", "-", names(asv.data))
tax <- read.table("tax.tab", header = TRUE)

genera <- unique(as.character(tax$Genus))
genus.data <- matrix(NA, nrow = nrow(asv.data), ncol = length(genera))

for (i in 1 : length(genera)) 
{
  print(i)
  genus <- genera[i]
  if(!is.na(genus)) {
    asv.genus <- rownames(tax)[tax$Genus == genus & !is.na(tax$Genus)]
  } else {
    asv.genus <- rownames(tax)[is.na(tax$Genus)]
  }
  # print(asv.genus)
  print(dim(data.frame(asv.data[, c(asv.genus)])))
  genus.data[, i] <- apply(data.frame(asv.data[, c(asv.genus)]), 1, sum)
}

genus.data <- data.frame(genus.data)
names(genus.data) <- genera
names(genus.data)[is.na(names(genus.data))] <- "NONE"
rownames(genus.data) <- rownames(asv.data)

# Verify row sums of genus data are identical to asv data
sum(apply(genus.data, 1, sum) == apply(asv.data, 1, sum))

####################################

fish <- genus.data
fish <- rownames_to_column(fish, var = "ID")
info <- read.csv("metadata.tab")
info <- rownames_to_column(info, var = "ID")

fish_full <- inner_join(fish, info, by = c("ID" = "ID"))

# Remove fish with 0 DaysPE as they are not independent samples of the rest
fish_full <- subset(fish_full, DaysPE != 0)
col_to_drop = colnames(info)
fish_full_num_only = fish_full %>% select(-one_of(col_to_drop))

# Remove genus that exits in less than 5% of the samples
sum(colSums(fish_full_num_only > 0) > dim(fish_full_num_only)[1] * 0.05)
fish_full_filtered <- fish_full_num_only[, colSums(fish_full_num_only > 0) > dim(fish_full_num_only)[1] * 0.05]
fish_infected = fish_full_filtered[fish_full$Total>0 & !is.na(fish_full$Total), ]
fish_not_infected = fish_full_filtered[fish_full$Total==0 | is.na(fish_full$Total), ]

# Swap reference to the last column
library(dplyr)
x_not_infected = relocate(fish_not_infected, "NONE", .after=last_col())
x_infected = relocate(fish_infected, "NONE", .after=last_col())
saveRDS(colnames(x_infected), "genus_names.rds")
rownames(x_not_infected)=NULL
colnames(x_not_infected)=NULL
rownames(x_infected)=NULL
colnames(x_infected)=NULL
saveRDS(as.matrix(x_infected), "infected.rds")
saveRDS(as.matrix(x_not_infected), "not_infected.rds")
