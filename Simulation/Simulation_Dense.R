# Simulation: ROC for dense data (Section 3.2)
rm(list = ls())
library(snow)

source("Functions_ROC.R")
start_time <- proc.time()

n_rep <- 100  # number of replicates
n_cpu <- 10   # number of cores for parallel computing
length_rholist <- 70

for(type in 1 : 3)
{
  for(seq_depth in c("l", "h"))
  {
    for(z_var in c("l", "h"))
    {
      seed_list <- sample(1:n_rep, n_cpu, replace = FALSE)
      cluster <- makeCluster(n_cpu, type = "SOCK", outfile = "")
      res_list <- clusterApply(cluster, seed_list, Simulation_ROC, 
                               n = 100, K = 200, length_rholist = length_rholist, 
                               type = type, seq_depth = seq_depth, 
                               z_var = z_var, n_rep = n_rep, n_cpu = n_cpu)
      stopCluster(cluster)
      saveRDS(res_list, paste("T", type, "_n_100_K_200_", seq_depth, 
                              "_", z_var, ".rds", sep = ""))
    }
  }
}

end_time <- proc.time()
print(end_time - start_time)

# Plot: ROC for dense data (Section 3.2)

library(tidyverse)

T <- list()
T[[1]] <- T[[2]] <- T[[3]] <- list()

for(type in 1 : 3)
{
  j <- 1
  for(seq_depth in c("h", "l"))
  {
    for(z_var in c("h", "l"))
    {
      results <- read_rds(paste("T", type, "_n_100_K_200_", seq_depth, "_", z_var, ".rds", sep = ""))
      T[[type]][[j]] <- Plot_ROC(n_rep = n_rep, n_cpu = n_cpu, length_rholist = length_rholist) +
        labs(caption = paste(seq_depth, ", ", z_var, sep = "")) + 
        theme(legend.position="none", plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      results <- NULL
      gc()
      j <- j + 1
    }
  }
}

type.legend <- c("Chain", "Random", "Hub")
plot <- list()
for(type in 1 : 3)
{
  panel <- plot_grid(T[[type]][[1]], T[[type]][[2]], T[[type]][[3]], T[[type]][[4]], nrow = 1)
  title <- ggdraw() + draw_label(type.legend[type], size = 20, fontface = 'bold')
  plot[[type]] <- plot_grid(title, panel, ncol = 1, rel_heights = c(0.2, 1.5)) 
}

pdf(file = "ROC_dense.pdf", width = 14, height = 12)
p <- plot_grid(plot[[1]], plot[[2]], plot[[3]], nrow = 3, rel_heights = c(2, 2, 2))
print(p)
dev.off()

# Simulation: Bar graph for dense data (Section 3.2)
rm(list = ls())

source("Functions_Bar_Graph.R")
start_time <- proc.time()

n_rep <- 50  # number of replicates
n_cpu <- 10  # number of cores for parallel computing
results = list()

for(type in 1 : 3)
{
  for(seq_depth in c("l", "h"))
  {
    for(z_var in c("l", "h"))
    {
      for (i in 1 : n_rep) {
        cat(i, "-th replicate", "\n")
        results[[i]] <- Simulation_Bar_Graph(seed = i, n = 100, K = 200, 
                                             type = type, seq_depth = seq_depth, z_var = z_var, 
                                             length_rholist = 20, n_rep = n_rep, n_cpu = n_cpu)
      }
      saveRDS(results, paste("T", type, "_n_100_K_200_", seq_depth, 
                             "_", z_var, "_StARS.rds", sep = ""))
    }
  }
}

end_time <- proc.time()
print(end_time - start_time)

# Plot: Bar graph for dense data (Section 3.2)

library(tidyverse)

methods <- c("Comp-gLASSO", "gLASSO", "MB")
types <- 1 : 3
seq_depths <- c("h", "l")
z_vars <- c("h", "l")

Summary_Bar_Graph(methods, types, seq_depths, z_vars, n_rep)
results <- readRDS("Results_StARS.rds")
Plot_Bar_Graph(results, types)



