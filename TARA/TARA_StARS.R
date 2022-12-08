# Aim 2: real data 
# rm(list = ls())
# Load necessary libraries
gc()
library(MASS)
library(glasso)
library(huge)
library(propagate)
library(foreach)
library(doParallel)
library(space)
source("../CompoGlasso.R")

# Load screened and pre-processed data
x <- readRDS("TARA_reference.rds")

start_time <- proc.time()
n <- dim(x)[1] #number of samples
K <- dim(x)[2] - 1  #number of OTUs (reference OTU not counted)
M <- as.numeric(apply(x, 1, sum))

num_cores_StARS <- 1
num_cores_NR <- 1
iter_sub <- 100
n_sub = floor(10 * sqrt(n))
beta <- 0.05
length_rholist <- 30

# Perform log-ratio transformation
p.hat <- colMeans(x/M)
offset = K + 1
x.adj <- t(t(x) + p.hat * offset)
z.hat <- log(x.adj[, -(K + 1)]/x.adj[, (K + 1)])

# Glasso with stability selection
Sigma.2 <- bigcor(z.hat, fun = "cov", verbose = FALSE)
Sigma.2 <- Sigma.2[1:nrow(Sigma.2), 1:ncol(Sigma.2)]
rho.list.2 <- exp(seq(log(max(Sigma.2)), log(max(Sigma.2) / 100), length = length_rholist))
cat("rho.list.2:", "\n")
print(rho.list.2)

n2.rho <- length(rho.list.2)
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
  
  # StaARS selection
  lmda = rho.list.2[i.rho.2]
  
  # Parallel the StARS selection subsampling:
  cl<-makeCluster(num_cores_StARS)
  registerDoParallel(cl)
  Omegas.2.sub <- foreach (i_iter_sub = 1:iter_sub, combine = cbind, .export = "obj", .packages = c("huge", "propagate")) %dopar%
  {
    z.hat.sub <- z.hat[sample(1:n, size = n_sub), ]
    Sigma.2.sub <- bigcor(z.hat.sub, fun = "cov", verbose = FALSE)
    Sigma.2.sub <- Sigma.2.sub[1:nrow(Sigma.2.sub), 1:ncol(Sigma.2.sub)]
    Omega.2.sub <- as.matrix(huge(Sigma.2.sub, method = "glasso", lambda = lmda)$icov[[1]])
    Omega.2.sub
  }
  stopCluster(cl)
  
  Omegas.2.sub.tmp = array(dim = c(K, K, iter_sub))
  for (i in 1:iter_sub){
    Omegas.2.sub.tmp[, , i] = Omegas.2.sub[[i]]
  }
  Omegas.2.sub = Omegas.2.sub.tmp
  
  
  Omega.2.prob <- rowMeans((Omegas.2.sub > 0) * 1, dims = 2, na.rm = TRUE) 
  Omega.2.ksi <- 2 * Omega.2.prob * (1 - Omega.2.prob)
  Omegas.2.prob[, , i.rho.2] <- Omega.2.prob
  Omegas.2.ksi[, , i.rho.2] <- Omega.2.ksi
  D_b.2[i.rho.2] <- (sum(Omega.2.ksi) / 2) / (K * (K - 1) / 2)
  sup_D_b.2[i.rho.2] <- max(D_b.2)
  Omegas.2.orig[, , i.rho.2] <- as.matrix(Omegas.2[[i.rho.2]])
  if (sup_D_b.2[i.rho.2] > beta) break
}
rho.2 <- rho.list.2[i.rho.2 - 1]
Omega.2 <- Omegas.2.orig[, , i.rho.2 - 1]
non_empty_2 <- ((sum(Omega.2 != 0) - K) / 2) / (K * (K - 1) / 2)
non_empty_2

# MB with stability selection
Sigma.3 <- bigcor(z.hat, fun = "cov", verbose = FALSE)
Sigma.3 <- Sigma.3[1:nrow(Sigma.3), 1:ncol(Sigma.3)]
rho.list.3 <- exp(seq(log(max(Sigma.3) * 2), log(max(Sigma.3) / 1000), length = length_rholist))
cat("rho.list.3:", "\n")
print(rho.list.3)

n3.rho <- length(rho.list.3)
Omegas.3.orig <- array(0, dim = c(K, K, n3.rho))
Omegas.3.prob <- array(0, dim = c(K, K, n3.rho))
Omegas.3.ksi <- array(0, dim = c(K, K, n3.rho))
D_b.3 <- array(0, n3.rho)
sup_D_b.3 <- array(0, n3.rho)

for(i.rho.3 in 1 : n3.rho)
{
  cat("i.rho.3 = :", i.rho.3, "\n")
  
  # StaARS selection
  lmda = rho.list.3[i.rho.3]
  
  # Parallel the StARS selection subsampling:
  cl<-makeCluster(num_cores_StARS)
  registerDoParallel(cl)
  Omegas.3.sub <- foreach (i_iter_sub = 1:iter_sub, combine = cbind, .export = "obj", .packages = c("huge", "propagate", "space")) %dopar%
    {
      z.hat.sub <- z.hat[sample(1:n, size = n_sub), ]
      Sigma.3.sub <- bigcor(z.hat.sub, fun = "cov", verbose = FALSE)
      Sigma.3.sub <- Sigma.3.sub[1:nrow(Sigma.3.sub), 1:ncol(Sigma.3.sub)]
      Omega.3.sub <- space.neighbor(z.hat.sub, lam1 = lmda)$ParCor
      Omega.3.sub
    }
  stopCluster(cl)
  Omegas.3.sub.tmp = array(dim = c(K, K, iter_sub))
  for (i in 1:iter_sub){
    Omegas.3.sub.tmp[, , i] = Omegas.3.sub[[i]]
  }
  Omegas.3.sub = Omegas.3.sub.tmp
  
  Omega.3.prob <- rowMeans((Omegas.3.sub > 0) * 1, dims = 2) 
  Omega.3.ksi <- 2 * Omega.3.prob * (1 - Omega.3.prob)
  Omegas.3.prob[, , i.rho.3] <- Omega.3.prob
  Omegas.3.ksi[, , i.rho.3] <- Omega.3.ksi
  D_b.3[i.rho.3] <- (sum(Omega.3.ksi) / 2) / (K * (K - 1) / 2)
  sup_D_b.3[i.rho.3] <- max(D_b.3)
  Omegas.3.orig[, , i.rho.3] <- space.neighbor(z.hat, lam1 = lmda)$ParCor
  if (sup_D_b.3[i.rho.3] > beta) break
}
rho.3 <- rho.list.3[i.rho.3 - 1]
Omega.3 <- Omegas.3.orig[, , i.rho.3 - 1]
non_empty_3 <- ((sum(Omega.3 != 0) - K) / 2) / (K * (K - 1) / 2)
non_empty_3

# Compo-glasso with StARS selection
rho.list.1 <- exp(seq(log(max(Sigma.2) * 1), log(max(Sigma.2) / 300), length = length_rholist))
n1.rho <- length(rho.list.1)
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
  
  cl<-makeCluster(num_cores_StARS)
  registerDoParallel(cl)
  Omegas.1.sub <- foreach (i_iter_sub = 1:iter_sub, combine = cbind, .export = "obj", .packages = c("huge", "propagate")) %dopar%
  {
    source("../CompoGlasso.R")
    cat("StARS selection subsampling number =", i_iter_sub, "\n")
    Omega.0 <- Omegas.1.orig[, , i.rho.1]
    x_sub <- x.adj[sample(1:n, size = n_sub), ]
    Omega.1.sub <- Compo_glasso(x_sub, option = 0, rho.list = rho, O_ratio = 100, z_ratio = 100, max_iter = 100)$Omega.1
    Omega.1.sub
  }
  stopCluster(cl)
  Omegas.1.sub.tmp = array(dim = c(K, K, iter_sub))
  for (i in 1:iter_sub){
    Omegas.1.sub.tmp[, , i] = Omegas.1.sub[[i]]
  }
  Omegas.1.sub = Omegas.1.sub.tmp
  
  Omega.1.prob <- rowMeans((Omegas.1.sub > 0) * 1, dims = 2) 
  Omegas.1.prob[, , i.rho.1] <- Omega.1.prob
  Omega.1.ksi <- 2 * Omega.1.prob * (1 - Omega.1.prob)
  Omegas.1.ksi[, , i.rho.1] <- Omega.1.ksi
  D_b.1[i.rho.1] <- (sum(Omega.1.ksi) / 2) / (K * (K - 1) / 2)
  sup_D_b.1[i.rho.1] <- max(D_b.1)

  if (sup_D_b.1[i.rho.1] > beta) break
  
}

rho.1 <- rho.list.1[i.rho.1 - 1]
Omega.1 <- Omegas.1.orig[, , i.rho.1 - 1]
non_empty_1 <- ((sum(Omega.1 != 0) - K) / 2) / (K * (K - 1) / 2)
non_empty_1

results <- list(Omega.1, Omegas.1.orig, i.rho.1 - 1, rho.list.1, non_empty_1, D_b.1, sup_D_b.1, Omegas.1.prob, Omegas.1.ksi, 
                Omega.2, Omegas.2.orig, i.rho.2 - 1, rho.list.2, non_empty_2, D_b.2, sup_D_b.2, Omegas.2.prob, Omegas.2.ksi,
                Omega.3, Omegas.3.orig, i.rho.3 - 1, rho.list.3, non_empty_3, D_b.3, sup_D_b.3, Omegas.3.prob, Omegas.3.ksi)
names(results) <- c("Omega.1", "Omegas.1.orig", "i.rho.1", "rho.list.1", "non_empty_1", "D_b.1", "sup_D_b.1", "Omegas.1.prob", "Omegas.1.ksi",
                    "Omega.2", "Omegas.2.orig", "i.rho.2", "rho.list.2", "non_empty_2", "D_b.2", "sup_D_b.2", "Omegas.2.prob", "Omegas.2.ksi",
                    "Omega.3", "Omegas.3.orig", "i.rho.3", "rho.list.3", "non_empty_3", "D_b.3", "sup_D_b.3", "Omegas.3.prob", "Omegas.3.ksi")
saveRDS(results, "TARA_StARS.rds") 
end_time <- proc.time()
print(end_time - start_time)