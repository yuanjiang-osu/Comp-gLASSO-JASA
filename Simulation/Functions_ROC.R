Simulation_ROC <- function(seed, n = 100, K = 200, length_rholist = 70, type = 1, 
                           seq_depth = "h", z_var = "h", n_rep = 100, n_cpu = 10)  
{
  gc()
  # n_rep <- 10
  # n_cpu <- 5
  library(MASS)
  library(glasso)
  library(huge)
  source("../CompoGlasso.R")
  
  # n <- 100 # number of samples
  # K <- 200 # number of OTUs (reference OTU not counted)
  # length_rholist <- 70 # number of penalty parameters
  
  # Type of inverse covariance matrix
  # type <- 1 # chain
  # type <- 3 # random
  # type <- 4 # cluster
  
  # Size of sequencing depth
  # seq_depth <- 'h' 
  # seq_depth <- 'l'
  # seq_depth <- '1'
  # seq_depth <- '2'
  # seq_depth <- '3'
  # seq_depth <- '4'
  
  # Size of multivariate normal variance
  # z_var <- 'h'
  # z_var <- 'l'
  
  offset <- K + 1   ##20170406:To see if this influence a lot
  option <- 2 # Doesn't account for the cases when reference OTU = 0?

  results <- list() ##Setting list to save results
  
  ## data generation
  ###Setting for loop:
  for (n_pa in 1:(n_rep / n_cpu)) {
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
      M <- runif(n, 1 * K, 2 * K)
    } else if (seq_depth == '3') {
      M <- runif(n, 2 * K, 4 * K) 
    } else if (seq_depth == '2') {
      M <- runif(n, 4 * K, 8 * K)
    } else if (seq_depth == '1') {
      M <- runif(n, 8 * K, 16 * K)
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
    path.2 <- huge(Sigma.2, method = "glasso")
    cat("Glasso path created \n")
    
    rho.list.2 <- exp(seq(log(max(Sigma.2) * 2), log(max(Sigma.2) / 1000), length = length_rholist))
    
    # SPEIC-EASI gl
    cat("rho.list.2:", "\n")
    print(rho.list.2)
    
    n2.rho <- length(rho.list.2)
    TP.2_s <- rep(NA, n2.rho)
    FP.2_s <- rep(NA, n2.rho)
    ER.2_s <- rep(NA, n2.rho)
    RC.2_s <- rep(NA, n2.rho)
    
    path.2 <- huge(Sigma.2, method = "glasso", lambda = rho.list.2)
    Omegas.2 <- path.2$icov
    for(i.rho.2 in 1 : n2.rho)
    {
      cat("i.rho.2 = :", i.rho.2, "\n")
      Omega.2 <- as.matrix(Omegas.2[[i.rho.2]])
      TP.2_s[i.rho.2] <- sum((Omega.2 - diag(diag(Omega.2))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.True != 0) - K)
      FP.2_s[i.rho.2] <- sum(Omega.2 != 0 & Omega.True == 0)/sum(Omega.True == 0) #False positive rates
      ER.2_s[i.rho.2] <- sum((Omega.2 - Omega.True)^2)  #Estimation errors
      RC.2_s[i.rho.2] <- sum((Omega.2 - diag(diag(Omega.2))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.2 != 0) - K) #Precision
    }
    
    # SPEIC-EASI mb
    rho.list.3 <- exp(seq(log(max(Sigma.2) * 2), log(max(Sigma.2) / 4000), length = length_rholist))
    
    cat("rho.list.3:", "\n")
    print(rho.list.3)
    
    n3.rho <- length(rho.list.3)
    TP.3_s <- rep(NA, n3.rho)
    FP.3_s <- rep(NA, n3.rho)
    ER.3_s <- rep(NA, n3.rho)
    RC.3_s <- rep(NA, n3.rho)
    
    Omegas.3 <- array(0, dim = c(K, K, n3.rho))
    path.3 <- huge(Sigma.2, method = "mb", lambda = rho.list.3)
    Omegas.adj.3 <- path.3$path
    for(i.rho.3 in 1 : n3.rho)
    {
      cat("i.rho.3 = :", i.rho.3, "\n")
      Omega.3 <- as.matrix(Omegas.adj.3[[i.rho.3]]) + diag(K)
      Omegas.3[, , i.rho.3] = Omega.3
      TP.3_s[i.rho.3] <- sum((Omega.3 - diag(diag(Omega.3))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.True != 0) - K)
      FP.3_s[i.rho.3] <- sum(Omega.3 != 0 & Omega.True == 0)/sum(Omega.True == 0) #False positive rates
      ER.3_s[i.rho.3] <- sum((Omega.3 - Omega.True)^2)  #Estimation errors
      RC.3_s[i.rho.3] <- sum((Omega.3 - diag(diag(Omega.3))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.3 != 0) - K) #Precision
    }
    
    # Compo-glasso
    rho.list.1 = rho.list.2
    n1.rho <- length(rho.list.1)
    TP.1_s <- rep(NA, n1.rho)
    FP.1_s <- rep(NA, n1.rho)
    ER.1_s <- rep(NA, n1.rho)
    RC.1_s <- rep(NA, n1.rho)
    
    Omegas.1 <- array(0, dim = c(K, K, n1.rho))
    flag.1 <- rep(0, n1.rho)
    Sigmas_rho <- list()
    Omegas_rho <- list()
    for(i.rho.1 in 1 : n1.rho)    
    {
      cat("i.rho.1 =", i.rho.1, "\n")
      rho <- rho.list.1[i.rho.1]
      Sigmas_iter <- list()
      Omegas_iter <- list()
      
      Sigma.0 <- cov(z.hat)
      if (i.rho.1 == 1) {
        Omega.0 <- as.matrix(Omegas.2[[length_rholist]])
      }
      else if (i.rho.1 > 1) {
        Omega.0 <- as.matrix(Omegas.2[[1]])
      }
      Omega.1 <- matrix(0, K, K)
      iter <- 0
      z.0 <- z.hat
      z.1 <- matrix(0, n, K)
      mu.0 <- apply(z.0, 2, mean)
      z.start <- z.hat
      z.end <- matrix(0, n, K)
      while (mean((Omega.0 - Omega.1) ^ 2) > 0.0000001 || mean((z.start - z.end) ^ 2) > 0.1)
      {
        cat("iter = ", iter + 1, "mean((Omega.0 - Omega.1) ^ 2) = ", mean((Omega.0 - Omega.1) ^ 2), "\n")
        if (iter != 0) {
          Omega.0 <- Omega.1
        }
        
        if (mean((z.start - z.end) ^ 2) > 0.1) {
          if (iter != 0) {
            z.start <- z.end
          }
          iter <- iter + 1
          # cat("iter = ", iter, "\n")
          z.0 <- z.start
          mu.0 <- apply(z.0, 2, mean)
          z.1 <- matrix(0, n, K)
          z.iter <- 0
          while (mean((z.0 - z.1) ^ 2) > 0.0000001){
            if (z.iter != 0) {
              z.0 <- z.1
            }
            z.iter <- z.iter + 1
            for (j in 1:n) {
              dipi <- M[j] * exp(z.0[j,]) / (t(rep(1, K)) %*% exp(z.0[j,]) + 1) - x[j, 1:K] +  as.vector(Omega.0 %*% (z.0[j,] - mu.0))
              tripi <- M[j] * diag(exp(z.0[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) - 
                M[j] * (exp(z.0[j,])) %*% t(exp(z.0[j,])) / (as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1)) ^ 2 + Omega.0
              z.1[j,] <- z.0[j,] - solve(tripi) %*% dipi
            }
            cat("z iteration = ", z.iter, "mean square difference in z's = ", mean((z.0 - z.1) ^ 2), "\n")
          }
          z.end <- z.1
          Sigma.1 <- cov(z.end)
          
          mod <- huge(x = Sigma.1, lambda = rho, method = "glasso")
          
          Omega.1 = as.matrix(mod$icov[[1]])
          
          Sigmas_iter[[iter]] <- Sigma.1
          Omegas_iter[[iter]] <- Omega.1
        }
      }
      
      Sigmas_rho[[i.rho.1]] <- Sigmas_iter
      Omegas_rho[[i.rho.1]] <- Omegas_iter
      
      Omegas.1[, , i.rho.1] <- Omega.1
      TP.1_s[i.rho.1] <- sum((Omega.1 - diag(diag(Omega.1))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.True != 0) - K)
      FP.1_s[i.rho.1] <- sum(Omega.1 != 0 & Omega.True == 0)/sum(Omega.True == 0) #False positive rates
      ER.1_s[i.rho.1] <- sum((Omega.1 - Omega.True)^2)  #Estimation errors
      RC.1_s[i.rho.1] <- sum((Omega.1 - diag(diag(Omega.1))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.1 != 0) - K)
    }
    
    results[[n_pa]] <- list(Omega.True, Omegas.1, Omegas.2, Omegas.3,
                            TP.1_s, TP.2_s, TP.3_s,  
                            FP.1_s, FP.2_s, FP.3_s, 
                            ER.1_s, ER.2_s, ER.3_s, 
                            RC.1_s, RC.2_s, RC.3_s,
                            rho.list.1, rho.list.2, rho.list.3)
    names(results[[n_pa]]) <- c("Omega.True", "Omegas.1", "Omegas.2", "Omegas.3",
                                "TP.1_s", "TP.2_s", "TP.3_s", 
                                "FP.1_s", "FP.2_s", "FP.3_s", 
                                "ER.1_s", "ER.2_s", "ER.3_s", 
                                "RC.1_s", "RC.2_s", "RC.3_s", 
                                "rho.list.1", "rho.list.2", "rho.list.3")
  }
  return(results)
}


Plot_ROC <- function(n_rep, n_cpu, length_rholist)
{
  library(cowplot)
  # n_rep <- 100
  # n_cpu <- 10
  # length_rholist <- 70
  
  # ROC curve of TP/FP from all rhos, in one simulation
  
  TP.1_array <- array(0, c(n_cpu, n_rep / n_cpu, length_rholist))
  TP.2_array <- array(0, c(n_cpu, n_rep / n_cpu, length_rholist))
  TP.3_array <- array(0, c(n_cpu, n_rep / n_cpu, length_rholist))

  FP.1_array <- array(0, c(n_cpu, n_rep / n_cpu, length_rholist))
  FP.2_array <- array(0, c(n_cpu, n_rep / n_cpu, length_rholist))
  FP.3_array <- array(0, c(n_cpu, n_rep / n_cpu, length_rholist))
  
  for (i in 1:n_cpu) {
    for (j in 1:(n_rep / n_cpu)) {
      for (k in 1:length_rholist) {
        TP.1 <- results[[i]][[j]]$TP.1_s[k]
        TP.2 <- results[[i]][[j]]$TP.2_s[k]
        TP.3 <- results[[i]][[j]]$TP.3_s[k]

        FP.1 <- results[[i]][[j]]$FP.1_s[k]
        FP.2 <- results[[i]][[j]]$FP.2_s[k]
        FP.3 <- results[[i]][[j]]$FP.3_s[k]
        
        TP.1_array[i, j, k] <- TP.1
        TP.2_array[i, j, k] <- TP.2
        TP.3_array[i, j, k] <- TP.3

        FP.1_array[i, j, k] <- FP.1
        FP.2_array[i, j, k] <- FP.2
        FP.3_array[i, j, k] <- FP.3
      }
    }
  }
  
  avg_TP.1 <- apply(TP.1_array, 3, mean, na.rm = TRUE) 
  avg_TP.2 <- apply(TP.2_array, 3, mean, na.rm = TRUE) 
  avg_TP.3 <- apply(TP.3_array, 3, mean, na.rm = TRUE) 

  avg_FP.1 <- apply(FP.1_array, 3, mean, na.rm = TRUE) 
  avg_FP.2 <- apply(FP.2_array, 3, mean, na.rm = TRUE)
  avg_FP.3 <- apply(FP.3_array, 3, mean, na.rm = TRUE) 
  
  roc <- data.frame(TPR = c(avg_TP.1, avg_TP.2, avg_TP.3), FPR = c(avg_FP.1, avg_FP.2, avg_FP.3), 
                    method = c(rep("Comp-gLASSO", length_rholist), rep("gLASSO", length_rholist), rep("MB", length_rholist)))
  roc_curve <- ggplot(data = roc, aes(x = FPR, y = TPR, group = method)) +
    geom_path(aes(colour = method, linetype = method)) +
    scale_color_manual(name = "Method", breaks = c("Comp-gLASSO", "gLASSO", "MB"), values =c("red", "green", "blue")) +
    scale_linetype_manual(name = "Method", breaks= c("Comp-gLASSO", "gLASSO", "MB"), values=c("solid", "longdash", "dotted")) + 
    xlim(0, 1) + 
    ylim(0, 1.05)
  
  return(roc_curve)
}