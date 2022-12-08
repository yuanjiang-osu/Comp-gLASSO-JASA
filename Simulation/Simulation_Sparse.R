# Simulation: ROC for sparse data (Section 3.3)
rm(list = ls())
library(snow)

source("Functions_ROC.R")
start_time <- proc.time()

n_rep <- 100  # number of replicates
n_cpu <- 10  # number of cores for parallel computing
length_rholist <- 70
type <- 1

for(seq_depth in c("1", "2", "3", "4"))
{
  for(z_var in c("l", "h"))
  {
    seed_list <- sample(1 : n_rep, n_cpu, replace = FALSE)
    cluster <- makeCluster(n_cpu, type = "SOCK", outfile = "")
    res_list <- clusterApply(cluster, seed_list, Simulation_ROC, 
                             n = 100, K = 200, length_rholist = length_rholist, 
                             type = type, seq_depth = seq_depth, z_var = z_var, 
                             n_rep = n_rep, n_cpu = n_cpu)
    stopCluster(cluster)
    saveRDS(res_list, paste("T", type, "_n_100_K_200_", seq_depth, 
                            "_", z_var, ".rds", sep = ""))
  }
}

end_time <- proc.time()
print(end_time - start_time)

# Plot: ROC for sparse data (Section 3.3)

library(tidyverse)

T <- list()
T[[1]] <- T[[2]] <- list()

for(z_var in c("h", "l"))
{
  i <- which(z_var == c("h", "l"))
  for(seq_depth in c("1", "2", "3", "4"))
  {
    j <- which(seq_depth == c("1", "2", "3", "4"))
    results <- read_rds(paste("T", type, "_n_100_K_200_", seq_depth, "_", z_var, ".rds", sep = ""))
    T[[i]][[j]] <- Plot_ROC(n_rep = n_rep, n_cpu = n_cpu, length_rholist = length_rholist) +
      labs(caption = paste("Sparsity Level ", seq_depth, sep = "")) + 
      theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14, face = 'bold'),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    results <- NULL
    gc()
  }
}

z_var.legend <- c("High Compositional Variation", "Low Compositional Variation")
plot <- list()
for(z_var in c("h", "l"))
{
  i <- which(z_var == c("h", "l"))
  panel <- plot_grid(T[[i]][[1]], T[[i]][[2]], T[[i]][[3]], T[[i]][[4]], nrow = 1)
  title <- ggdraw() + draw_label(z_var.legend[i], size = 20, fontface = 'bold')
  plot[[i]] <- plot_grid(title, panel, ncol = 1, rel_heights = c(0.2, 1.5)) 
}

pdf(file = "ROC_sparse.pdf", width = 10.5, height = 6)
p <- plot_grid(plot[[1]], plot[[2]], nrow = 2, rel_heights = c(2, 2))
print(p)
dev.off()

# Simulation: Bar graph for sparse data (Section 3.3)

rm(list = ls())

source("Functions_Bar_Graph.R")
start_time <- proc.time()

n_rep <- 50  # number of replicates
n_cpu <- 10  # number of cores for parallel computing
type <- 1
results = list()

for(seq_depth in c("1", "2", "3", "4"))
{
  for(z_var in c("l", "h"))
  {
    for (i in 1:n_rep) {
      cat(i, "-th replicate", "\n")
      results[[i]] <- Simulation_Bar_Graph(seed = i, n = 100, K = 200, 
                                           type = type, seq_depth = seq_depth, z_var = z_var, 
                                           length_rholist = 20, n_rep = n_rep, n_cpu = n_cpu)
    }
    saveRDS(results, paste("T", type, "_n_100_K_200_", seq_depth, 
                           "_", z_var, "_StARS.rds", sep = ""))
  }
}

end_time <- proc.time()
print(end_time - start_time)

# Plot: Bar graph for sparse data (Section 3.3)

methods <- c("Comp-gLASSO", "gLASSO", "MB")
types <- 1
seq_depths <- c("1", "2", "3", "4")
z_vars <- c("h", "l")

Summary_Bar_Graph(methods, types, seq_depths, z_vars, n_rep)
results <- readRDS("Results_StARS.rds")
Plot_Bar_Graph_Sparse(results, z_vars)

