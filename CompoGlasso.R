# Function to generate network (inverse covariance matrix) in the simulations
generate_cov <- function(K, type)
{
  library(huge)
  if(type == 1){
    # Chain network
    Omega.True <- matrix(0, nrow = K, ncol = K)
    for(i in 1 : K) {
      Omega.True[i, i] = 1.5
      for(j in 1 : K) {
        if (abs(i - j) == 1) Omega.True[i, j] = 0.5
      }
    }
    Sigma.True <- solve(Omega.True)
    Sigma.True[abs(Sigma.True) < 1e-10] <- 0
  }
  else if(type == 2){
    # Random network
    if (K == 200) {
    Omega.True <- as.matrix(read.table("T2_K200_3_over_K.txt"))
    }
    Omega.True = matrix(Omega.True, ncol = K,byrow=TRUE)
    Sigma.True <- solve(Omega.True)
  }
  else if(type == 3){
    # Hub network
    L <- huge.generator(graph = "hub", vis = FALSE, d = K)
    Sigma.True <- L$sigma
    Omega.True <- L$omega
    Omega.True[abs(Omega.True) < 1e-10] <- 0
  }
  return(list(Sigma.True = Sigma.True, Omega.True = Omega.True))
}

# Function for additive log-ratio transformation
# Option give different options to handle zero's in the numerator/denominator in the ratio
# Only option = 0 or 2 is used in our work
z_hat_offset <- function(x, offset, option)
{
  K <- dim(x)[2] - 1
  M <- apply(x, 1, sum)
  N <- dim(x)[1]
  
  # no offset
  if(option == 0){
    z.hat <- log(x[, -(K + 1)]/x[, (K + 1)])
  }
  # x offset by a constant
  if(option == 1){
    p.hat <- rep(1/(K + 1), K + 1)
    x.adj <- t(t(x) + p.hat * offset)
    z.hat <- log(x.adj[, -(K + 1)]/x.adj[, (K + 1)])
  }
  # x offset proportionally
  else if(option == 2){
    p.hat <- colMeans(x/M)
    x.adj <- t(t(x) + p.hat * offset)
    z.hat <- log(x.adj[, -(K + 1)]/x.adj[, (K + 1)])
  }
  # ratio offset proportionally only for zero x_j or x_{K+1}
  else if(option == 3){
    p.hat <- colMeans(x/M)
    z.hat <- matrix(NA, nrow = nrow(x), ncol = K)
    for(j in 1 : K){
      zero.ind <- (x[, j] == 0) | (x[, K + 1] == 0)
      z.hat[!zero.ind, j] <- log(x[!zero.ind, j]/x[!zero.ind, (K + 1)])
      z.hat[zero.ind, j] <- log(p.hat[j]/p.hat[K + 1])
    }
  }
  # ratio offset proportionally for all x's
  else if(option == 4){
    p.hat <- colMeans(x/M)
    z.hat <- matrix(NA, nrow = nrow(x), ncol = K)
    for(j in 1 : K){
      zero.ind <- (x[, j] == 0) | (x[, K + 1] == 0)
      z.hat[!zero.ind, j] <- log(x[!zero.ind, j]/x[!zero.ind, (K + 1)] + p.hat[j]/p.hat[K + 1])
      z.hat[zero.ind, j] <- log(p.hat[j]/p.hat[K + 1])
    }
  }
  # ratio offset proportionally only for zero x_j or x_{K+1}
  else if(option == 5){
    p.hat <- matrix(0, N, K + 1)
    z.hat <- matrix(NA, nrow = nrow(x), ncol = K)
    for (i in 1:N) {
        zeros <- which(x[i, ] == 0)
        nzeros <- which(x[i, ] != 0)
        p.hat[i, zeros] <- (x[i, zeros] + offset)/(M[i] + offset * length(zeros))
        p.hat[i, nzeros] <- (x[i, nzeros])/(M[i] + offset * length(zeros))
        z.hat[i, ] <- log(p.hat[i, -(K + 1)]/p.hat[i, K + 1])
      }
  }
  return(z.hat)
}

# Function to calculate the objective function in equation (7) of the paper
obj <- function(x, z, Omega, K) {
  M <- sum(x)
  mu = mean(z)
  f = 1 / 2 * t(z - mu) %*% Omega %*% (z - mu) - (t(x) %*% z - M * log(as.numeric(t(rep(1, K)) %*% exp(z) + 1))) 
  return(as.numeric(f))
}


# Newton-Raphson procedure to optimize the objective function in equation (7) of the paper (with parallelization)
NR_para <- function(x, z.0, Omega.0, alpha_0 = 1, delta = 5 , epsilon = 0.01, threshold = 0.0001, num_cores = 1)
{
  # Initialization
  n = dim(z.0)[1]
  K = dim(z.0)[2]
  M <- as.numeric(apply(x, 1, sum))
  mu.0 = apply(z.0, 2, mean)
  z.1 <- matrix(0, n, K)
  
  # Parallelization begins
  library(foreach)
  library(doParallel)
  cl<-makeCluster(num_cores)
  registerDoParallel(cl)
  z.new <- foreach (j = 1:n, .combine = rbind, .export = "obj") %dopar% {
    z.iter <- 0
    alpha <- alpha_0
    # Loop to update z
    while (mean((z.0[j,] - z.1[j,]) ^ 2) > threshold) {
      if (z.iter != 0 && (obj(x[j, 1:K], z.1[j,], Omega.0, K) <= (obj(x[j, 1:K], z.0[j,], Omega.0, K) + epsilon * alpha * h_0))) {
        z.0[j,] <- z.1[j,]
      }
      z.iter <- z.iter + 1
      # Gradient of the objective function in equation (7) with respect to z
      dipi <- M[j] * exp(z.0[j,]) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) - x[j, 1:K] +  as.vector(Omega.0 %*% (z.0[j,] - mu.0))
      # Hessian matrix of the objective function in equation (7) with respect to z
      tripi <- M[j] * diag(exp(z.0[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) - 
        M[j] * (exp(z.0[j,])) %*% t(exp(z.0[j,])) / (as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1)) ^ 2 + Omega.0
      # Newton-Raphson's updating rule
      z.1[j,] <- z.0[j,] - alpha * solve(tripi) %*% dipi
      # Shrink the step size using Armijo's Rule
      dk = (-1) * solve(tripi) %*% dipi
      h_0 = t(z.0[j,] - mu.0) %*% Omega.0 %*% dk - x[j, 1:K] %*% dk + 
        M[j] * as.numeric(t(dk) %*% exp(z.1[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.1[j,]) + 1)
      if (obj(x[j, 1:K], z.1[j,], Omega.0, K) > (obj(x[j, 1:K], z.0[j,], Omega.0, K) + epsilon * alpha * h_0)) {
        alpha <- alpha / delta
      }
    }
    z.1[j,]
  }
  stopCluster(cl)
  
  return(unname(z.new, force = TRUE))
}

# Newton-Raphson procedure to optimize the objective function in equation (7) of the paper
NR <- function(x, z.0, Omega.0, alpha_0 = 1, delta = 5 , epsilon = 0.01, threshold = 0.0001)
{
  # Initialization
  n = dim(z.0)[1]
  K = dim(z.0)[2]
  M <- as.numeric(apply(x, 1, sum))
  mu.0 = apply(z.0, 2, mean)
  z.1 <- matrix(0, n, K)
  
  for (j in 1:n) 
  {
    z.iter <- 0
    alpha <- alpha_0
    # Loop to update z
    while (mean((z.0[j,] - z.1[j,]) ^ 2) > threshold) {
      if (z.iter != 0 && (obj(x[j, 1:K], z.1[j,], Omega.0, K) <= (obj(x[j, 1:K], z.0[j,], Omega.0, K) + epsilon * alpha * h_0))) {
        z.0[j,] <- z.1[j,]
      }
      z.iter <- z.iter + 1
      # Gradient of the objective function in equation (7) with respect to z
      dipi <- M[j] * exp(z.0[j,]) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) - x[j, 1:K] +  as.vector(Omega.0 %*% (z.0[j,] - mu.0))
      # Hessian matrix of the objective function in equation (7) with respect to z
      tripi <- M[j] * diag(exp(z.0[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) - 
        M[j] * (exp(z.0[j,])) %*% t(exp(z.0[j,])) / (as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1)) ^ 2 + Omega.0
      # Newton-Raphson's updating rule
      z.1[j,] <- z.0[j,] - alpha * solve(tripi) %*% dipi
      # Shrink the step size using Armijo's Rule
      dk = (-1) * solve(tripi) %*% dipi
      h_0 = t(z.0[j,] - mu.0) %*% Omega.0 %*% dk - x[j, 1:K] %*% dk + 
        M[j] * as.numeric(t(dk) %*% exp(z.0[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1)
      if (obj(x[j, 1:K], z.1[j,], Omega.0, K) > (obj(x[j, 1:K], z.0[j,], Omega.0, K) + epsilon * alpha * h_0)) {
        alpha <- alpha / delta
      }
    }
  }
  return(z.1)
}

# Main function for compositional graphical lasso
Compo_glasso <- function(x, rho.list, option = 2, offset = (K + 1), para_NR = FALSE, num_cores_NR = 1, z_ratio = 1000, O_ratio = 1000, max_iter = 50)  
{
  gc()
  library(MASS)
  library(glasso)
  library(huge)
  library(propagate)
  
  n <- dim(x)[1] # number of samples
  K <- dim(x)[2] - 1 # number of OTUs (reference OTU not counted)
  M <- as.numeric(apply(x, 1, sum)) # Sequencing depths
  
  results <- list() # List to save results
  
  # Additive log-ratio transformation
  z.hat <- z_hat_offset(x, offset, option)
  colnames(z.hat) = NULL
  rownames(z.hat) = NULL
  cat("z_hat comupted \n")
  
  rho.list.1  <- rho.list
  n1.rho <- length(rho.list.1)
  Sigma.2 <- bigcor(z.hat, fun = "cov", verbose = FALSE)
  Sigma.2 <- Sigma.2[1:nrow(Sigma.2), 1:ncol(Sigma.2)]
  
  # Run graphical lasso on the additive log-ratio transformed data
  path.2 <- huge(Sigma.2, method = "glasso", lambda = max(rho.list))
  Omegas.2 <- path.2$icov
  Sigma.0 <- Sigma.2
  Omegas.1 <- array(0, dim = c(K, K, n1.rho))
  
  for(i.rho.1 in 1 : n1.rho)
  {
    if (n1.rho != 1) cat("i.rho.1 =", i.rho.1, "\n")
    rho <- rho.list.1[i.rho.1]
    
    # Initialization
    if (i.rho.1 == 1) {
      Omega.0 <- as.matrix(Omegas.2[[1]]) 
    }
    else if (i.rho.1 > 1) {
      Omega.0 <- Omegas.1[, , (i.rho.1 - 1)]
    }
    
    # Iteration between Newton-Raphson and graphical lasso
    Omega.1 <- matrix(0, K, K)
    iter <- 0
    z.start <- z.hat
    z.end <- matrix(0, n, K)
    O_thr <- mean((Omega.0 - Omega.1) ^ 2) / O_ratio # convergence threshold
    z_thr <- mean((z.start - z.end) ^ 2) / z_ratio # convergence threshold
    while ((mean((Omega.0 - Omega.1) ^ 2) > O_thr || (mean((z.start - z.end) ^ 2) > z_thr)) && iter <= max_iter)
    {
      cat("iter = ", iter + 1, "mean((Omega.0 - Omega.1) ^ 2) = ", mean((Omega.0 - Omega.1) ^ 2), "\n")
      if (iter != 0) {
        Omega.0 <- Omega.1
        z.start <- z.end
      }
      iter <- iter + 1
      
      # Update z: Newton-Raphson
      if (para_NR == TRUE) z.end <- NR_para(x, z.start, Omega.0, num_cores = num_cores_NR) # with parallelization
      else z.end <- NR(x, z.start, Omega.0) # without parallelization
      cat("mean square z.end - z.start =", mean((z.start - z.end) ^ 2), "\n")
      
      # Update Omega: Graphical Lasso
      Sigma.1 <- bigcor(z.end, fun = "cov", verbose = FALSE)
      Sigma.1 <- Sigma.1[1:nrow(Sigma.1), 1:ncol(Sigma.1)]
      mod <- huge(x = Sigma.1, lambda = rho, method = "glasso")
      Omega.1 = as.matrix(mod$icov[[1]])
      
      # Added 04/08/2022 by Yuan Jiang to avoid ERROR
      if (is.na(mean((Omega.0 - Omega.1) ^ 2)) || is.na(mean((z.start - z.end) ^ 2)) ) 
        break
    }
    
    # Updated 04/08/2022 by Yuan Jiang to avoid ERROR
    # Omegas.1[, , i.rho.1] <- Omega.1
    if (any(is.na(Omega.1))) {
      Omegas.1[, , i.rho.1] <- Omega.0
    } else {
      Omegas.1[, , i.rho.1] <- Omega.1
    }
    
  }
  
  results <- list(Omega.1, Omegas.1, rho.list.1)
  names(results) <- c("Omega.1", "Omegas.1", "rho.list.1")
  return(results)
}