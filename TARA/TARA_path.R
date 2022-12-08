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
n <- dim(x)[1] # number of samples
K <- dim(x)[2] - 1  # number of OTUs (reference OTU not counted)
M <- as.numeric(apply(x, 1, sum))

num_cores_StARS <- 1
num_cores_NR <- 1
iter_sub <- 100
n_sub = floor(10 * sqrt(n))
beta <- 0.05
length_rholist <- 70

# Perform log-ratio transformation
z.hat <- z_hat_offset(x, offset = K + 1, option = 2)

# Glasso
Sigma.2 <- bigcor(z.hat, fun = "cov", verbose = FALSE)
Sigma.2 <- Sigma.2[1:nrow(Sigma.2), 1:ncol(Sigma.2)]
rho.list.2 <- exp(seq(log(max(Sigma.2)), log(14.0952), length = length_rholist))
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
  
  lmda = rho.list.2[i.rho.2]
  
  Omegas.2.orig[, , i.rho.2] <- as.matrix(Omegas.2[[i.rho.2]])
}
Omega.2 <- Omegas.2.orig[, , i.rho.2]
non_empty_2 <- ((sum(Omega.2 != 0) - K) / 2) / (K * (K - 1) / 2)
non_empty_2
((sum(Omega.2 != 0) - K) / 2)

# Neighbourhood selection
Sigma.3 <- bigcor(z.hat, fun = "cov", verbose = FALSE)
Sigma.3 <- Sigma.3[1:nrow(Sigma.3), 1:ncol(Sigma.3)]
rho.list.3 <- exp(seq(log(max(Sigma.3) / 10), log(0.2387954), length = length_rholist)) # 0.16:300
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
  
  lmda = rho.list.3[i.rho.3]
  
  Omegas.3.orig[, , i.rho.3] <- space.neighbor(z.hat, lam1 = lmda)$ParCor
}
Omega.3 <- Omegas.3.orig[, , i.rho.3]
non_empty_3 <- ((sum(Omega.3 != 0) - K) / 2) / (K * (K - 1) / 2)
non_empty_3
((sum(Omega.3 != 0) - K) / 2)

# Compo-glasso
rho.list.1 <- exp(seq(log(max(Sigma.2) * 1), log(13.19391), length = length_rholist))
n1.rho <- length(rho.list.1)
Omegas.1 <- array(0, dim = c(K, K, n1.rho))

Omegas.1.orig <- Compo_glasso(x, rho.list = rho.list.1, O_ratio = 100, z_ratio = 100)$Omegas.1

i.rho.1 <- n1.rho
Omega.1 <- Omegas.1.orig[, , i.rho.1]
non_empty_1 <- ((sum(Omega.1 != 0) - K) / 2) / (K * (K - 1) / 2)
non_empty_1
((sum(Omega.1 != 0) - K) / 2)

results <- list(Omegas.1.orig, Omegas.2.orig, Omegas.3.orig, length_rholist,
                rho.list.1, rho.list.2, rho.list.3)
names(results) <- c("Omegas.1.orig", "Omegas.2.orig", "Omegas.3.orig","length_rholist",
                    "rho.list.1", "rho.list.2", "rho.list.3")
saveRDS(results, "TARA_path.rds") 
end_time <- proc.time()
print(end_time - start_time)
