library(tidyverse)

x <- readRDS("TARA_reference.rds")

n <- dim(x)[1]
K <- dim(x)[2] - 1

genera <- readRDS("TARA_genus_names.rds")

results <- readRDS("TARA_path.rds")

pairs_lit <- readRDS("TARA_pairs_literature.rds")

genus_1 <- matrix(rep(genera, times = K), nrow = K, ncol = K)
genus_2 <- matrix(rep(genera, each = K), nrow = K, ncol = K)
row_index <- matrix(rep(seq(1, K), times = K), nrow = K, ncol = K)
col_index <- matrix(rep(seq(1, K), each = K), nrow = K, ncol = K)
genus_1.unique <- genus_1[row_index > col_index]
genus_2.unique <- genus_2[row_index > col_index]

pairs_lit_sym <- rbind(cbind(pairs_lit$genus_org1, pairs_lit$genus_org2), cbind(pairs_lit$genus_org2, pairs_lit$genus_org1))
pairs_lit_sym_zip = rep(c("0"), dim(pairs_lit_sym)[1])
for (i in 1:dim(pairs_lit_sym)[1]) {
  pairs_lit_sym_zip[i] = paste(pairs_lit_sym[i, 1], pairs_lit_sym[i, 2])
}
Omega_lit <- data.frame(genus_1 = as.character(genus_1.unique), genus_2 = as.character(genus_2.unique))
in_lit <- rep(0, dim(Omega_lit)[1])
for (i in 1:dim(Omega_lit)[1]) {
  in_lit[i] = paste(as.character(Omega_lit$genus_1[i]), as.character(Omega_lit$genus_2[i])) %in% pairs_lit_sym_zip
}
sum(in_lit)
Omega_lit <- data.frame(genus_1 = as.character(genus_1.unique), genus_2 = as.character(genus_2.unique), lit = in_lit)

# Compo-glasso
n1.rho <- results$length_rholist
selected.1 <- rep(NA, n1.rho)
selected.lit.1 <- rep(NA, n1.rho)
for (i.rho.1 in 1:n1.rho) {
  Omega.1.ksi <- results$Omegas.1.orig[, , i.rho.1]
  Omega.1.ksi.unique <- Omega.1.ksi[row_index > col_index]
  Omega.1.unique <- (Omega.1.ksi != 0)[row_index > col_index]
  Omega.1.ksi.all <- data.frame(instability = Omega.1.ksi.unique, selected = Omega.1.unique, genus_1 = genus_1.unique, genus_2 = genus_2.unique, lit = in_lit)
  Omega.1.StARS.selected <- Omega.1.ksi.all %>% filter(selected == TRUE)
  selected.1[i.rho.1] <- dim(Omega.1.StARS.selected)[1]
  Omega.1.StARS.selected.lit <- Omega.1.ksi.all %>% filter(selected == TRUE & lit == 1)
  selected.lit.1[i.rho.1] <- dim(Omega.1.StARS.selected.lit)[1]
}

# glasso
n2.rho <- results$length_rholist
selected.2 <- rep(NA, n2.rho)
selected.lit.2 <- rep(NA, n2.rho)
for (i.rho.2 in 1:n2.rho) {
  Omega.2.ksi <- results$Omegas.2.orig[, , i.rho.2]
  Omega.2.ksi.unique <- Omega.2.ksi[row_index > col_index]
  Omega.2.unique <- (Omega.2.ksi != 0)[row_index > col_index]
  Omega.2.ksi.all <- data.frame(instability = Omega.2.ksi.unique, selected = Omega.2.unique, genus_1 = genus_1.unique, genus_2 = genus_2.unique, lit = in_lit)
  # Omega.1.ksi.selected <- Omega.1.ksi.all[Omega.1.ksi.all$selected, ]
  Omega.2.StARS.selected <- Omega.2.ksi.all %>% filter(selected == TRUE)
  selected.2[i.rho.2] <- dim(Omega.2.StARS.selected)[1]
  Omega.2.StARS.selected.lit <- Omega.2.ksi.all %>% filter(selected == TRUE & lit == 1)
  selected.lit.2[i.rho.2] <- dim(Omega.2.StARS.selected.lit)[1]
}

# MB
n3.rho <- results$length_rholist
selected.3 <- rep(NA, n3.rho)
selected.lit.3 <- rep(NA, n3.rho)
for (i.rho.3 in 1:n3.rho) {
  Omega.3.ksi <- results$Omegas.3.orig[, , i.rho.3]
  Omega.3.ksi.unique <- Omega.3.ksi[row_index > col_index]
  Omega.3.unique <- (Omega.3.ksi != 0)[row_index > col_index]
  Omega.3.ksi.all <- data.frame(instability = Omega.3.ksi.unique, selected = Omega.3.unique, genus_1 = genus_1.unique, genus_2 = genus_2.unique, lit = in_lit)
  Omega.3.StARS.selected <- Omega.3.ksi.all %>% filter(selected == TRUE)
  selected.3[i.rho.3] <- dim(Omega.3.StARS.selected)[1]
  Omega.3.StARS.selected.lit <- Omega.3.ksi.all %>% filter(selected == TRUE & lit == 1)
  selected.lit.3[i.rho.3] <- dim(Omega.3.StARS.selected.lit)[1]
}

library(ggplot2)
method = rep(c("Compo-glasso", "glasso", "mb"), each = results$length_rholist)
selected = c(selected.1, selected.2, selected.3)
selected.lit = c(selected.lit.1, selected.lit.2, selected.lit.3)
d = data.frame(method, selected, selected.lit)

pdf(file = "TARA_path_half.pdf", width = 6, height = 6)
p <- ggplot(data = d, aes(x = selected, y = selected.lit)) +
  geom_step(aes(color = method, linetype = method)) + 
  scale_linetype_manual(values=c("solid", "dashed", "dotted"), name = "Method", labels = c('Compo-glasso', "glasso", "MB")) +
  labs(x = "Number of Selected Edges", y = "Number of Selected Edges Reported in the Literature") +
  scale_color_discrete(name = "Method", labels = c('Compo-glasso', "glasso", "MB")) + 
  scale_y_continuous(breaks = seq(0, 15, by = 2), limits = c(0, 14)) + 
  xlim(0, 200) +
  theme_bw() +
  theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()
