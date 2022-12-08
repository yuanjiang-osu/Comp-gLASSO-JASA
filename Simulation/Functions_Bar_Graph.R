Simulation_Bar_Graph <- function(seed, n = 100, K = 50, length_rholist = 20, type = 1, seq_depth = "l", z_var = "l", 
                                 n_rep = 50, n_cpu = 10)  
{
  gc()
  library(MASS)
  library(glasso)
  library(huge)
  library(propagate)
  library(foreach)
  library(doParallel)
  source("../CompoGlasso.R")
  
  # n <- 100  # number of samples
  # K <- 50   # number of OTUs (reference OTU not counted)
  # K <- 200
  
  # Type of inverse covariance matrix
  # type <- 1 Chain
  # type <- 2 Random
  # type <- 3 Cluster
  
  # Size of sequencing depth
  # seq_depth <- 'l'
  # seq_depth <- 'h'
  # seq_depth <- '1'
  # seq_depth <- '2'
  # seq_depth <- '3'
  # seq_depth <- '4'
  
  # Size of multivariate normal variance
  # z_var <- 'l'
  # z_var <- 'h'
  
  offset <- K + 1   # Laplace smoothing
  option <- 2
  
  num_cores_NR <- 1 # Parallel processes in Newton's method
  # n_cpu = 10 Parallel processes in StARS
  # n_rep = 50 Number of StARS subsample
  n_sub = floor(7 * sqrt(n)) # size of StARS subsample
  beta <- 0.05 # StARS parameter
  
  results <- list() ##Setting list to save results
  
  ## data generation
  ###Setting for loop:
  temp <- generate_cov(K, type)
  Sigma.True <- temp$Sigma.True
  Omega.True <- temp$Omega.True
  
  avg <- 0
  mu <- rep(avg, K)
  if (z_var == 'l') err <- 1 else if (z_var == 'h') err <- 5
  
  z <- matrix(mvrnorm(n, mu, err * Sigma.True), nrow = n, ncol = K)
  
  p <- matrix(0, nrow = n, ncol = K + 1)
  p[, -(K + 1)] <- exp(z)/((apply(exp(z), 1, sum) + 1) %*% matrix(1, ncol = K))
  p[, (K + 1)] <- 1/(apply(exp(z), 1, sum) + 1)
  cat("Sum of Multinomial variance:", sum(apply(p, MARGIN = 1, FUN = var)), "\n")
  
  # Setting sequencing depth
  if (seq_depth == 'l') {
    M <- runif(n, 20 * K, 40 * K) 
  } else if (seq_depth == 'h') {
    M <- runif(n, 100 * K, 200 * K)
  } else if (seq_depth == '4') {
    M <- runif(n, 0.5 * K, 1 * K)
  } else if (seq_depth == '3') {
    M <- runif(n, 1 * K, 2 * K) 
  } else if (seq_depth == '2') {
    M <- runif(n, 2 * K, 4 * K)
  } else if (seq_depth == '1') {
    M <- runif(n, 4 * K, 8 * K)
  }
  cat("M's simulated \n")
  
  x <- matrix(0, n, K + 1)
  for(i in 1 : n)
  {
    x[i, ] <- rmultinom(1, size = M[i], prob = p[i, ])
  }
  cat("X's initiated \n")
  
  z.hat <- z_hat_offset(x, offset, option)
  Sigma.2 <- cov(z.hat)
  cat("Sigma 2 comupted \n")
  
  # Glasso
  rho.list.2 <- exp(seq(log(max(Sigma.2) * 2), log(max(Sigma.2) / 100), length = length_rholist))
  cat("rho.list.2:", "\n")
  print(rho.list.2)
  
  n2.rho <- length(rho.list.2)
  TP.2_s <- rep(NA, n2.rho)
  FP.2_s <- rep(NA, n2.rho)
  ER.2_s <- rep(NA, n2.rho)
  PR.2_s <- rep(NA, n2.rho)
  
  path.2 <- huge(Sigma.2, method = "glasso", lambda = rho.list.2)
  Omegas.2 <- path.2$icov
  Omegas.2.orig <- array(0, dim = c(K, K, n2.rho))
  Omegas.2.prob <- array(0, dim = c(K, K, n2.rho))
  Omegas.2.ksi <- array(0, dim = c(K, K, n2.rho))
  D_b.2 <- array(0, n2.rho)
  sup_D_b.2 <- array(0, n2.rho)
  
  for(i.rho.2 in 1 : n2.rho)
  {
    cat("i.rho.2 = :", i.rho.2, "\n")
    
    # StARS selection
    lmda = rho.list.2[i.rho.2]
    
    # Parallel the StARS selection subsampling:
    cl<-makeCluster(n_cpu)
    registerDoParallel(cl)
    Omegas.2.sub <- foreach (i_rep = 1:n_rep, combine = cbind, .export = "obj", .packages = c("huge", "propagate")) %dopar%
      {
        source("../CompoGlasso.R")
        x_sub <- x[sample(1:n, size = n_sub), ]
        z.hat.sub <- z_hat_offset(x_sub, offset = K + 1, option = 2)
        Sigma.2.sub <- bigcor(z.hat.sub, fun = "cov", verbose = FALSE)
        Sigma.2.sub <- Sigma.2.sub[1:nrow(Sigma.2.sub), 1:ncol(Sigma.2.sub)]
        Omega.2.sub <- as.matrix(huge(Sigma.2.sub, method = "glasso", lambda = lmda)$icov[[1]])
        Omega.2.sub
      }
    stopCluster(cl)
    Omegas.2.sub <- array(as.numeric(unlist(Omegas.2.sub)), dim=c(K, K, n_rep))
    
    Omega.2.prob <- rowMeans((Omegas.2.sub > 0) * 1, dims = 2) 
    Omega.2.ksi <- 2 * Omega.2.prob * (1 - Omega.2.prob)
    Omegas.2.prob[, , i.rho.2] <- Omega.2.prob
    Omegas.2.ksi[, , i.rho.2] <- Omega.2.ksi
    D_b.2[i.rho.2] <- (sum(Omega.2.ksi) / 2) / (K * (K - 1) / 2)
    sup_D_b.2[i.rho.2] <- max(D_b.2)
    Omegas.2.orig[, , i.rho.2] <- as.matrix(Omegas.2[[i.rho.2]])
    
    # Updated 04/08/2022 by Yuan Jiang to avoid ERROR
    if (is.na(sup_D_b.2[i.rho.2])) {
      break
    } else {
      if (sup_D_b.2[i.rho.2] > beta) break
    }
  }
  rho.2 <- rho.list.2[i.rho.2 - 1]
  Omega.2 <- Omegas.2.orig[, , i.rho.2 - 1]
  non_empty_2 <- ((sum(Omega.2 != 0) - K)/2) / (K * (K - 1) / 2)
  non_empty_2
  
  TP.2 <- sum((Omega.2 - diag(diag(Omega.2))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.True != 0) - K) #Also recall
  FP.2 <- sum(Omega.2 != 0 & Omega.True == 0)/sum(Omega.True == 0) #False positive rates
  PR.2 <- sum((Omega.2 - diag(diag(Omega.2))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.2 != 0) - K) # Precision
  TP.2
  FP.2
  PR.2
  F1.2 <- ifelse(TP.2 == 0 | PR.2 == 0, 0, 2 * TP.2 * PR.2/(TP.2 + PR.2))
  F1.2
  
  # Neighbourhood Selection
  rho.list.3 <- exp(seq(log(max(Sigma.2) * 2), log(max(Sigma.2) / 100), length = length_rholist))
  cat("rho.list.3:", "\n")
  print(rho.list.3)
  
  n3.rho <- length(rho.list.3)
  TP.3_s <- rep(NA, n3.rho)
  FP.3_s <- rep(NA, n3.rho)
  ER.3_s <- rep(NA, n3.rho)
  PR.3_s <- rep(NA, n3.rho)
  
  Omegas.3.orig <- array(0, dim = c(K, K, n2.rho))
  Omegas.3.prob <- array(0, dim = c(K, K, n2.rho))
  Omegas.3.ksi <- array(0, dim = c(K, K, n2.rho))
  D_b.3 <- array(0, n3.rho)
  sup_D_b.3 <- array(0, n3.rho)
  
  for(i.rho.3 in 1 : n3.rho)
  {
    cat("i.rho.3 = :", i.rho.3, "\n")
    
    # StaARS selection
    lmda = rho.list.3[i.rho.3]
    
    # Parallel the StARS selection subsampling:
    cl<-makeCluster(n_cpu)
    registerDoParallel(cl)
    Omegas.3.sub <- foreach (i_rep = 1:n_rep, combine = cbind, .export = "obj", .packages = c("huge", "propagate")) %dopar%
      {
        source("../CompoGlasso.R")
        x_sub <- x[sample(1:n, size = n_sub), ]
        z.hat.sub <- z_hat_offset(x_sub, offset = K + 1, option = 2)
        Sigma.3.sub <- bigcor(z.hat.sub, fun = "cov", verbose = FALSE)
        Sigma.3.sub <- Sigma.3.sub[1:nrow(Sigma.3.sub), 1:ncol(Sigma.3.sub)]
        Omega.3.sub <- as.matrix(huge(Sigma.3.sub, method = "mb", lambda = lmda)$path[[1]]) + diag(K)
        Omega.3.sub
      }
    stopCluster(cl)
    Omegas.3.sub <- array(as.numeric(unlist(Omegas.3.sub)), dim=c(K, K, n_rep))
    
    Omega.3.prob <- rowMeans((Omegas.3.sub > 0) * 1, dims = 2) 
    Omega.3.ksi <- 2 * Omega.3.prob * (1 - Omega.3.prob)
    Omegas.3.prob[, , i.rho.3] <- Omega.3.prob
    Omegas.3.ksi[, , i.rho.3] <- Omega.3.ksi
    D_b.3[i.rho.3] <- (sum(Omega.3.ksi) / 2) / (K * (K - 1) / 2)
    sup_D_b.3[i.rho.3] <- max(D_b.3)
    Omegas.3.orig[, , i.rho.3] <- as.matrix(huge(Sigma.2, method = "mb", lambda = lmda)$path[[1]]) + diag(K)
    
    # Updated 04/08/2022 by Yuan Jiang to avoid ERROR
    if (is.na(sup_D_b.3[i.rho.3])) {
      break
    } else {
      if (sup_D_b.3[i.rho.3] > beta) break
    }
  }
  rho.3 <- rho.list.3[i.rho.3 - 1]
  Omega.3 <- Omegas.3.orig[, , i.rho.3 - 1]
  non_empty_3 <- ((sum(Omega.3 != 0) - K)/2) / (K * (K - 1) / 2)
  non_empty_3
  
  TP.3 <- sum((Omega.3 - diag(diag(Omega.3))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.True != 0) - K) #Also recall
  FP.3 <- sum(Omega.3 != 0 & Omega.True == 0)/sum(Omega.True == 0) #False positive rates
  PR.3 <- sum((Omega.3 - diag(diag(Omega.3))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.3 != 0) - K) # Precision
  TP.3
  FP.3
  PR.3
  
  F1.3 <- ifelse(TP.3 == 0 | PR.3 == 0, 0, 2 * TP.3 * PR.3/(TP.3 + PR.3))
  F1.3
  
  
  # Compo-glasso
  rho.list.1 <- exp(seq(log(max(Sigma.2) * 2), log(max(Sigma.2) / 100), length = length_rholist))
  
  n1.rho <- length(rho.list.1)
  TP.1_s <- rep(NA, n1.rho)
  FP.1_s <- rep(NA, n1.rho)
  ER.1_s <- rep(NA, n1.rho)
  PR.1_s <- rep(NA, n1.rho)
  Omegas.1 <- array(0, dim = c(K, K, n1.rho))
  
  Omegas.1.orig <- Compo_glasso(x, rho.list = rho.list.1)$Omegas.1
  Omegas.1.prob <- array(0, dim = c(K, K, n1.rho))
  Omegas.1.ksi <- array(0, dim = c(K, K, n1.rho))
  D_b.1 <- array(0, n1.rho)
  sup_D_b.1 <- array(0, n1.rho)
  
  for(i.rho.1 in 1 : n1.rho)    
  {
    cat("i.rho.1 =", i.rho.1, "\n")
    rho <- rho.list.1[i.rho.1]
    
    # n_sub = floor(n / 2)
    # n_rep <- 100
    # Omegas.1.sub <- array(0, dim = c(K, K, n_rep))
    cl<-makeCluster(n_cpu)
    registerDoParallel(cl)
    Omegas.1.sub <- foreach (i_rep = 1:n_rep, combine = cbind, .export = "obj", .packages = c("huge", "propagate")) %dopar%
      {
        cat("StARS selection subsampling number =", i_rep, "\n")
        source("../CompoGlasso.R")
        Omega.0 <- Omegas.1.orig[, , i.rho.1]
        x_sub <- x[sample(1:n, size = n_sub), ]
        Omega.1.sub <- Compo_glasso(x_sub, rho.list = rho)$Omega.1
        Omega.1.sub
      }
    stopCluster(cl)
    Omegas.1.sub <- array(as.numeric(unlist(Omegas.1.sub)), dim=c(K, K, n_rep))
    
    Omega.1.prob <- rowMeans((Omegas.1.sub > 0) * 1, dims = 2) 
    Omegas.1.prob[, , i.rho.1] <- Omega.1.prob
    Omega.1.ksi <- 2 * Omega.1.prob * (1 - Omega.1.prob)
    Omegas.1.ksi[, , i.rho.1] <- Omega.1.ksi
    D_b.1[i.rho.1] <- (sum(Omega.1.ksi) / 2) / (K * (K - 1) / 2)
    sup_D_b.1[i.rho.1] <- max(D_b.1)
    
    # Updated 04/08/2022 by Yuan Jiang to avoid ERROR
    if (is.na(sup_D_b.1[i.rho.1])) {
      break
    } else {
      if (sup_D_b.1[i.rho.1] > beta) break
    }
  }
  
  rho.1 <- rho.list.1[i.rho.1 - 1]
  Omega.1 <- Omegas.1.orig[, , i.rho.1 - 1]
  non_empty_1 <- ((sum(Omega.1 != 0) - K) / 2) / (K * (K - 1) / 2)
  non_empty_1
  
  TP.1 <- sum((Omega.1 - diag(diag(Omega.1))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.True != 0) - K)
  FP.1 <- sum(Omega.1 != 0 & Omega.True == 0)/sum(Omega.True == 0) #False positive rates
  PR.1 <- sum((Omega.1 - diag(diag(Omega.1))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.1 != 0) - K)
  TP.1
  FP.1
  PR.1
  
  F1.1 <- ifelse(TP.1 == 0 | PR.1 == 0, 0, 2 * TP.1 * PR.1/(TP.1 + PR.1))
  F1.1
  
  results <- list(Omega.True, Omega.1, Omegas.1.orig, i.rho.1 - 1, non_empty_1, D_b.1, sup_D_b.1, Omegas.1.prob, Omegas.1.ksi, 
                  Omega.2, Omegas.2.orig, i.rho.2 - 1, non_empty_2, D_b.2, sup_D_b.2, Omegas.2.prob, Omegas.2.ksi,
                  Omega.3, Omegas.3.orig, i.rho.3 - 1, non_empty_3, D_b.3, sup_D_b.3, Omegas.3.prob, Omegas.3.ksi,
                  TP.1, FP.1, PR.1, F1.1, TP.2, FP.2, PR.2, F1.2, TP.3, FP.3, PR.3, F1.3)
  names(results) <- c("Omega.True", "Omega.1", "Omegas.1.orig", "i.rho.1", "non_empty_1", "D_b.1", "sup_D_b.1", "Omegas.1.prob", "Omegas.1.ksi", 
                      "Omega.2", "Omegas.2.orig", "i.rho.2", "non_empty_2", "D_b.2", "sup_D_b.2", "Omegas.2.prob", "Omegas.2.ksi",
                      "Omega.3", "Omegas.3.orig", "i.rho.3", "non_empty_3", "D_b.3", "sup_D_b.3", "Omegas.3.prob", "Omegas.3.ksi",
                      "TP.1", "FP.1", "PR.1", "F1.1", "TP.2", "FP.2", "PR.2", "F1.2", "TP.3", "FP.3", "PR.3", "F1.3")
  return(results)
}


Summary_Bar_Graph <- function(methods, types, seq_depths, z_vars, n_rep)
{
  library(tidyverse)
  library(reshape2)
  
  # Count the number of settings/method combination
  n_method <- length(methods)
  n_type <- length(types)
  n_seq_depth <- length(seq_depths)
  n_var <- length(z_vars)
  n_comb <- n_method * n_var * n_seq_depth * n_type
  
  # Initialization
  TP <- rep(0, n_comb * n_rep)
  FP <- rep(0, n_comb * n_rep)
  PR <- rep(0, n_comb * n_rep)
  F1 <- rep(0, n_comb * n_rep)
  method.vec <- rep(0, n_comb * n_rep)
  type.vec <- rep(0, n_comb * n_rep)
  seq_depth.vec <- rep(0, n_comb * n_rep)
  z_var.vec <- rep(0, n_comb * n_rep)
  
  counter <- 0
  # T1, n = 100, K = 200, m, m
  for(type in types)
  {
    for(seq_depth in seq_depths)
    {
      for(z_var in z_vars)
      {
        results <- read_rds(paste("T", type, "_n_100_K_200_", seq_depth, "_", z_var, "_StARS.rds", sep = ""))
        for (i in 1:n_rep) {
          TP[i + counter * n_rep] <- results[[i]]$TP.1
          FP[i + counter * n_rep] <- results[[i]]$FP.1
          PR[i + counter * n_rep] <- results[[i]]$PR.1
          F1[i + counter * n_rep] <- results[[i]]$F1.1
          method.vec[i + counter * n_rep] <- "Comp-gLASSO"
          type.vec[i + counter * n_rep] <- type
          seq_depth.vec[i + counter * n_rep] <- seq_depth
          z_var.vec[i + counter * n_rep] <- z_var
        }
        counter <- counter + 1
        
        for (i in 1:n_rep) {
          TP[i + counter * n_rep] <- results[[i]]$TP.2
          FP[i + counter * n_rep] <- results[[i]]$FP.2
          PR[i + counter * n_rep] <- results[[i]]$PR.2
          F1[i + counter * n_rep] <- results[[i]]$F1.2
          method.vec[i + counter * n_rep] <- "gLASSO"
          type.vec[i + counter * n_rep] <- type
          seq_depth.vec[i + counter * n_rep] <- seq_depth
          z_var.vec[i + counter * n_rep] <- z_var
        }
        counter <- counter + 1
        
        for (i in 1:n_rep) {
          TP[i + counter * n_rep] <- results[[i]]$TP.3
          FP[i + counter * n_rep] <- results[[i]]$FP.3
          PR[i + counter * n_rep] <- results[[i]]$PR.3
          F1[i + counter * n_rep] <- results[[i]]$F1.3
          method.vec[i + counter * n_rep] <- "MB"
          type.vec[i + counter * n_rep] <- type
          seq_depth.vec[i + counter * n_rep] <- seq_depth
          z_var.vec[i + counter * n_rep] <- z_var
        }
        counter <- counter + 1
      }
    }
  }
  
  result_collection <- data.frame(type = type.vec, z_var = z_var.vec, seq_depth = seq_depth.vec,
                                  method = method.vec, TP = TP, FP = FP, PR = PR, F1 = F1)
  results <- result_collection %>% unite(var_read, c(seq_depth, z_var), sep = ", ", remove = FALSE)
  saveRDS(results, file = "Results_StARS.rds")
}


Plot_Bar_Graph <- function(results, types)
{
  library(tidyverse)
  library(reshape2)
  library(cowplot)
  
  type.legend <- c("Chain", "Random", "Hub")
  T <- list()
  T[[1]] <- T[[2]] <- T[[3]] <- list()
  
  for(type in types)
  {
    R <- results %>% filter(type == type) %>% group_by(method, var_read) %>% summarise(TP_mean = mean(TP), 
                                                                                       TP_sd = sd(TP), 
                                                                                       l_bar = max(0, TP_mean-TP_sd), 
                                                                                       h_bar = TP_mean+TP_sd)
    T_R <- ggplot(R, aes(x = var_read, y = TP_mean, fill = method, width = 0.8)) +
      geom_bar(stat = 'identity', position = 'dodge') + 
      geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                    position=position_dodge()) + 
      labs(x = "", y = "", fill = "", title = "Recall") +
      ylim(0, 1) +
      theme_light() + 
      theme(legend.position = "None", plot.title = element_text(hjust = 0.5, size = 16), 
            axis.text=element_text(size = 14))
    
    P <- results %>% filter(type == type) %>% group_by(method, var_read) %>% summarise(PR_mean = mean(PR), 
                                                                                       PR_sd = sd(PR), 
                                                                                       l_bar = max(0, PR_mean-PR_sd), 
                                                                                       h_bar = PR_mean+PR_sd)
    T_P <- ggplot(P, aes(x = var_read, y = PR_mean, fill = method, width = 0.8)) +
      geom_bar(stat = 'identity', position = 'dodge') + 
      geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                    position=position_dodge()) + 
      labs(x = "", y = "", fill = "", title = "Precision") +
      ylim(0, 1) +
      theme_light() + 
      theme(legend.position = "None", plot.title = element_text(hjust = 0.5, size = 16), 
            axis.text=element_text(size = 14))
    
    F <- results %>% filter(type == type) %>% group_by(method, var_read) %>% summarise(F1_mean = mean(F1), 
                                                                                       F1_sd = sd(F1), 
                                                                                       l_bar = max(0, F1_mean-F1_sd), 
                                                                                       h_bar = F1_mean+F1_sd)
    T_F <- ggplot(F, aes(x = var_read, y = F1_mean, fill = method, width = 0.8)) +
      geom_bar(stat = 'identity', position = 'dodge') + 
      geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                    position=position_dodge()) + 
      labs(x = "", y = "", fill = "", title = "F1") +
      ylim(0, 1) +
      theme_light() + 
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16),
            axis.text=element_text(size = 14))
    
    panel <- plot_grid(T_R, T_P, T_F, nrow = 1)
    title <- ggdraw() + draw_label(type.legend[type], size = 20, fontface = 'bold')
    T[[type]] <- plot_grid(title, panel, ncol = 1, rel_heights = c(0.2, 1.5)) 
  }
  pdf(file = "bar_graph_dense.pdf", width = 12, height = 12)
  p <- plot_grid(T[[1]], T[[2]], T[[3]], nrow = 3, rel_heights = c(2, 2, 2))
  print(p)
  dev.off()
}


Plot_Bar_Graph_Sparse <- function(results, z_vars)
{
  library(tidyverse)
  library(reshape2)
  library(cowplot)
  
  z_var.legend <- c("High Compositional Variation", "Low Compositional Variation")
  T <- list()
  T[[1]] <- T[[2]] <- list()
  
  for(z_var in z_vars)
  {
    i <- which(z_var == z_vars)
    R <- results %>% filter(z_var == z_var) %>% group_by(method, z_var, seq_depth) %>% summarise(TP_mean = mean(TP), 
                                                                                                         TP_sd = sd(TP), 
                                                                                                         l_bar = max(0, TP_mean-TP_sd), 
                                                                                                         h_bar = TP_mean+TP_sd)
    T_R <- ggplot(R, aes(x = seq_depth, y = TP_mean, fill = method, width = 0.8)) +
      geom_bar(stat = 'identity', position = 'dodge') + 
      geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                    position=position_dodge()) + 
      labs(x = "", y = "", fill = "", title = "Recall") +
      ylim(0, 1) +
      theme_light() + 
      theme(legend.position = "None", plot.title = element_text(hjust = 0.5, size = 16), 
            axis.text=element_text(size = 14))
    
    P <- results %>% filter(z_var == z_var) %>% group_by(method, z_var, seq_depth) %>% summarise(PR_mean = mean(PR), 
                                                                                                         PR_sd = sd(PR), 
                                                                                                         l_bar = max(0, PR_mean-PR_sd), 
                                                                                                         h_bar = PR_mean+PR_sd)
    T_P <- ggplot(P, aes(x = seq_depth, y = PR_mean, fill = method, width = 0.8)) +
      geom_bar(stat = 'identity', position = 'dodge') + 
      geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                    position=position_dodge()) + 
      labs(x = "", y = "", fill = "", title = "Precision") +
      ylim(0, 1) +
      theme_light() + 
      theme(legend.position = "None", plot.title = element_text(hjust = 0.5, size = 16), 
            axis.text=element_text(size = 14))
    
    F <- results %>% filter(z_var == z_var) %>% group_by(method, z_var, seq_depth) %>% summarise(F1_mean = mean(F1), 
                                                                                                         F1_sd = sd(F1), 
                                                                                                         l_bar = max(0, F1_mean-F1_sd), 
                                                                                                         h_bar = F1_mean+F1_sd)
    T_F <- ggplot(F, aes(x = seq_depth, y = F1_mean, fill = method, width = 0.8)) +
      geom_bar(stat = 'identity', position = 'dodge') + 
      geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                    position=position_dodge()) + 
      labs(x = "", y = "", fill = "", title = "F1") +
      ylim(0, 1) +
      theme_light() + 
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16),
            axis.text=element_text(size = 14))
    
    panel <- plot_grid(T_R, T_P, T_F, nrow = 1)
    title <- ggdraw() + draw_label(z_var.legend[i], size = 20, fontface = 'bold')
    T[[i]] <- plot_grid(title, panel, ncol = 1, rel_heights = c(0.2, 1.5)) 
  }
  pdf(file = "bar_graph_sparse.pdf", width = 16, height = 10)
  p <- plot_grid(T[[1]], T[[2]], nrow = 2, rel_heights = c(2, 2))
  print(p)
  dev.off()
}

