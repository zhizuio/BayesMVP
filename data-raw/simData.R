## code to prepare `simData` and `simDataTest` datasets goes here

requireNamespace("MASS", quietly = TRUE)
requireNamespace("BDgraph", quietly = TRUE)

# Parameters
n <- 100
p <- 200
q <- 20
snr <- 10 # may tune snr to 3 and 5

sim.mvp <- function(n, q, p, seed = 123) {
  # BayesMVP parameters
  covariancePrior <- "IW"
  gammaPrior <- "hierarchical"

  threshold_Gamma <- 0.5

  # X's correlation with two blocks
  Prime <- list(
    c(1:10), c(11:15)
  )

  Sigma <- diag(p)
  for (i in Prime) Sigma[i, i] <- 0.1

  ### Construction of G (dependencies about response variables)
  G <- matrix(0, q, q)

  Prime2 <- list(
    c(1:5),
    c(6:12),
    c(13:15),
    c(16:20)
  )
  Res <- Prime2

  for (i in Prime2) G[i, i] <- 1

  G[1:8, 19:20] <- 1
  G[15:20, 1:3] <- 1

  # X with additional correlations corresponding to residual graph
  for (i in 1:length(Prime2)) {
    Sigma[1:3 + (i + 3) * 10, 1:3 + (i + 3) * 10] <- 0.1
    Sigma[1:4 + (i + 1) * 30, 1:4 + (i + 1) * 30] <- 0.1
  }

  diag(Sigma) <- 1

  #### X with more correlations in irrelevant X's
  Sigma0 <- Sigma
  Sigma0[146:150, 146:150] <- 0.1
  diag(Sigma0) <- 1

  ############## Simulation of X
  mu <- c(rep(0, p))
  set.seed(1234 + seed)
  X <- MASS::mvrnorm(n, mu = mu, Sigma = Sigma0)

  ### Construction of B
  set.seed(12345)
  B <- matrix(rnorm(p * q, mean = 0, sd = 1), nrow = p)

  ### Construction of Gamma
  Gamma <- matrix(0, p, q)

  # Correlation present in Sigma
  Gamma[1:10, c(1:2, 8:10)] <- 1
  Gamma[11:15, 6:10] <- 1

  # Response variable that works together (G)
  for (i in 1:length(Prime2)) {
    Gamma[1:3 + (i + 3) * 10, Prime2[[i]]] <- 1
    Gamma[1:4 + (i + 1) * 30, Prime2[[i]]] <- 1
  }

  # Additional indep. covariates to be selected (random here), wide need to be 1 (1 x selected only)
  Gamma[183, 5:8] <- 1
  Gamma[173, 8:10] <- 1
  Gamma[193, 2:4] <- 1
  
  ###############################
  # Construction of MRF graph
  ###############################
  
  mrfG <- diag(q * p)
  
  ## use the adjacency matrix of the precision matrix of correlated X's
  Omega <- matrix(as.numeric(solve(Sigma) != 0), ncol = p)
  
  ## block X1:10 & Y1:2, Y8:10
  Omega1_10 <- Omega
  Omega1_10[-(1:10), ] <- 0  
  Omega1_10[, -(1:10)] <- 0 
  A.1_10 <- matrix(0, nrow=q, ncol=q)
  A.1_10[1:2, 1:2] <- 1
  A.1_10[8:10, 8:10] <- 1
  mrfG <- mrfG+kronecker(A.1_10,Omega1_10)
  
  
  #block X11:15 & Y6:10
  Omega11_15<-Omega
  Omega11_15[-(11:15), ] <- 0  
  Omega11_15[, -(11:15)] <- 0 
  A.11_15 <- matrix(0, nrow=q, ncol=q)
  A.11_15[6:10, 6:10] <- 1
  diag(A.11_15) <- 0
  mrfG<-mrfG+kronecker(A.11_15,Omega11_15)
  
  
  ##block X41-43 & Y1-5
  Omega41_43<-Omega
  Omega41_43[-(41:43), ] <- 0  
  Omega41_43[, -(41:43)] <- 0 
  A.41_43 <- matrix(0, nrow=q, ncol=q)
  A.41_43[1:5, 1:5] <- 1
  diag(A.41_43) <- 0
  mrfG<-mrfG+kronecker(A.41_43,Omega41_43)
  
  ##block X51:53 & Y6:12
  Omega51_53<-Omega
  Omega51_53[-(51:53), ] <- 0  
  Omega51_53[, -(51:53)] <- 0 
  A.51_53 <- matrix(0, nrow=q, ncol=q)
  A.51_53[6:12, 6:12] <- 1
  diag(A.51_53) <- 0
  mrfG<-mrfG+kronecker(A.51_53,Omega51_53)
  
  ##block X61-64 & Y1-5
  Omega61_64<-Omega
  Omega61_64[-(61:64), ] <- 0  
  Omega61_64[, -(61:64)] <- 0 
  A.61_64 <- matrix(0, nrow=q, ncol=q)
  A.61_64[1:5, 1:5] <- 1
  diag(A.61_64) <- 0
  mrfG<-mrfG+kronecker(A.61_64,Omega61_64)
  
  ##block X61-63 & Y13-15
  Omega61_63<-Omega
  Omega61_63[-(61:63), ] <- 0  
  Omega61_63[, -(61:63)] <- 0 
  A.61_63 <- matrix(0, nrow=q, ncol=q)
  A.61_63[13:15, 13:15] <- 1
  diag(A.61_63) <- 0
  mrfG<-mrfG+kronecker(A.61_63,Omega61_63)
  
  ##block X71:73 & Y16:20
  Omega71_73<-Omega
  Omega71_73[-(71:73), ] <- 0  
  Omega71_73[, -(71:73)] <- 0 
  A.71_73 <- matrix(0, nrow=q, ncol=q)
  A.71_73[16:20, 16:20] <- 1
  diag(A.71_73) <- 0
  mrfG<-mrfG+kronecker(A.71_73,Omega71_73)
  
  ##block X91:94 & Y6:12
  Omega91_94<-Omega
  Omega91_94[-(91:94), ] <- 0  
  Omega91_94[, -(91:94)] <- 0 
  A.91_94 <- matrix(0, nrow=q, ncol=q)
  A.91_94[6:12, 6:12] <- 1
  diag(A.91_94) <- 0
  mrfG<-mrfG+kronecker(A.91_94,Omega91_94)
  
  ##block X121:124 & Y13:15
  Omega121_124<-Omega
  Omega121_124[-(121:124), ] <- 0  
  Omega121_124[, -(121:124)] <- 0 
  A.121_124 <- matrix(0, nrow=q, ncol=q)
  A.121_124[13:15, 13:15] <- 1
  diag(A.121_124) <- 0
  mrfG<-mrfG+kronecker(A.121_124,Omega121_124)
  
  
  ## block X151:154 & Y16:20
  Omega151_154<-Omega
  Omega151_154[-(151:154), ] <- 0  
  Omega151_154[, -(151:154)] <- 0 
  A.151_154 <- matrix(0, nrow=q, ncol=q)
  A.151_154[16:20, 16:20] <- 1
  diag(A.151_154) <- 0
  mrfG<-mrfG+kronecker(A.151_154,Omega151_154)
  
  diag(mrfG) <- 0
  mrfG[upper.tri(mrfG)] <- 0
  
  # Summarize the true edges into a list
  mrfG <- which(mrfG==1, arr.ind=TRUE)
  ############################################################
  
  # We select the B coefficients (one selected) with Gamma
  B[Gamma == 0] <- 0

  ### XB
  XB <- X %*% B

  ### Finally, Y

  ### Simulate residuals
  set.seed(123456 + seed)
  U_tilde <- matrix(rnorm(n * q, mean = 0, sd = 1), nrow = n, ncol = q)

  # M
  v_r <- mean(diag(var(XB))) / snr

  M <- matrix(0.5, q, q)
  diag(M) <- rep(1, q)

  # set.seed(1234567 + seed)
  P <- BDgraph::rgwish(n = 1, adj = G, b = 3, D = v_r * M)
  varE <- solve(P)

  ### control signal-noise ratio

  factor <- 1.5
  factor_min <- 0.001
  factor_max <- 1000
  count <- 0
  maxit <- 10000

  factor_prev <- 1

  repeat{
    ### Sample the error
    varE <- varE / factor * factor_prev
    cVar <- chol(as.matrix(varE))
    U <- U_tilde %*% cVar
    S <- diag(diag(t(cVar)))
    sigma <- S * S

    ### Sample X0
    t0 <- 4
    X0 <- sample(1:t0, n, replace = T, prob = rep(1 / t0, t0))
    X0 <- model.matrix(~ factor(X0) + 0)[, ]
    set.seed(12 + seed)
    B0 <- matrix(rnorm(t0 * q, 0, 2), nrow = t0, ncol = q)

    # X0XB <- B0[X0] + XB
    X0XB <- X0 %*% B0 + XB

    ### Sample Y
    Z <- X0XB + U

    S <- diag(diag(cVar))
    sigma <- S * S
    ### S/N Ratio
    emp_snr <- mean(diag(var(X0XB) %*% solve(sigma)))
    ##############

    if (abs(emp_snr - snr) < (snr / 2) | count > maxit) {
      break
    } else {
      if (emp_snr < snr) { # increase factor
        factor_min <- factor
      } else { # decrease factor
        factor_max <- factor
      }
      factor_prev <- factor
      factor <- (factor_min + factor_max) / 2
    }
    count <- count + 1
  }
  cat("count = ", count, "\n")
  cat("emp_snr = ", emp_snr, "\n")

  colnames(X0) <- paste0("X0", seq_len(ncol(X0)))
  Y <- matrix(as.numeric(Z > 0), nrow = n, ncol = q)
  colnames(X) <- paste0("X", 1:p)
  colnames(Y) <- paste0("Y", 1:q)

  return(list(Y = Y, Z = Z, X = X, B = B, X0 = X0, B0 = B0, Gamma = Gamma, G = G, mrfG = mrfG))
}

simData <- sim.mvp(n, q, p, seed = 123)
simDataTest <- sim.mvp(n, q, p, seed = 1234)
# save(simData, file = "simData.rda")
# save(simDataTest, file = "simDataTest.rda")

usethis::use_data(simData, overwrite = TRUE)
usethis::use_data(simDataTest, overwrite = TRUE)
