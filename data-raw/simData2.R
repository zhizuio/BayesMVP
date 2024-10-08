## code to prepare `simData` and `simDataTest` datasets goes here

requireNamespace("MASS", quietly = TRUE)
requireNamespace("BDgraph", quietly = TRUE)

# Parameters
n <- 100
p <- 200
q <- 20
snr <- 10 # may tune snr to 3 and 5

sim.mvp2 <- function(n, q, p, seed = 123){
  
  threshold_Gamma<-0.5
  
  # X's correlation with two blocks
  Prime <- list(
    c(1:10),c(11:15))
  
  Sigma <- diag(p)
  for (i in Prime) {
    Sigma[i, i] <- 0.1
  }
  
  
  ###Construction of G (dependencies about response variables)
  G <- matrix(0, q, q)
  
  Prime2 <- list(
    c(1:5),
    c(6:12),
    c(13:15),
    c(16:20)
  )
  Res <- Prime2
  
  for (i in Prime2) {
    G[i, i] <- 1
  }
  
  G[1:8, 19:20] <- 1
  G[15:20, 1:3] <- 1
  
  # X with additional correlations corresponding to residual graph
  
  for (i in 1:length(Prime2)) {
    Sigma[1:3 +(i+3)*10,1:3 +(i+3)*10] <- 0.1
    Sigma[1:4 +(i+1)*30,1:4 +(i+1)*30] <- 0.1
  }
  
  diag(Sigma) <- 1
  
  
  #### X more correlations without being modelled
  Sigma0 <- Sigma
  Sigma0[146:150,146:150] <- 0.1
  diag(Sigma0) <- 1
  
  
  ####On simule alors X avec Sigma 0 mais le reste reste avce Sigma (mrfG par ex). Prior knowledge based on Sigma/Sigma0
  ###Ensuite, on peut supprimer des infos de Sigma (prior knowledge we don't know yet)
  
  
  ##############Simulation of X
  mu<-c(rep(0,p))
  
  #X
  set.seed(1234 + seed)
  X<-MASS::mvrnorm(n, mu=mu,Sigma=Sigma0)
  #X<-mvrnorm(n, mu=mu,Sigma=Sigma)
  
  ###Construction of B
  set.seed(12345)
  B<-matrix(rnorm(p*q, mean = 0,sd = 1), nrow = p)
  
  ### Random Gamma
  Gamma <- matrix(rbinom(p*q, 1, 0.1), p, q)
  
  #We select the B coefficients (one selected) with Gamma
  B[Gamma == 0] <- 0
  
  ###XB
  XB<-X %*% B
  
  ###Finally, Y
  
  ###Simulate residuals
  set.seed(123456 + seed)
  U <- matrix(rnorm(n*q, mean=0, sd=1), nrow=n, ncol=q)
  
  ### Sample X0
  t0 <- 4
  X0 <- sample(1:t0, n, replace=T, prob=rep(1/t0,t0))
  X0 <- model.matrix(~ factor(X0) + 0)[,]
  set.seed(12 + seed)
  B0 <- matrix(rnorm(t0*q, 0, 2), nrow = t0, ncol = q)
  
  #X0XB <- B0[X0] + XB
  X0XB <- X0 %*% B0 + XB
  
  ### Sample Y
  Z <- X0XB + U
  
  #S<-diag(diag(cVar))
  #sigma<-S*S
  ### S/N Ratio
  #emp_snr <- mean(diag(var(X0XB) %*% solve(sigma)))
  ####
  
  
  colnames(X0) <- paste0("X0", seq_len(ncol(X0)))
  Y <- matrix(as.numeric(Z > 0), nrow = n, ncol = q)
  colnames(X) <- paste0("X", 1:p)
  colnames(Y) <- paste0("Y", 1:q)
  
  return(list(Y = Y, Z = Z, X0 = X0, B0 = B0, X = X, B = B, Gamma = Gamma))
}

simData2 <- sim.mvp2(n, q, p, seed = 123)
simDataTest2 <- sim.mvp2(n, q, p, seed = 1234)
# save(simData2, file = "simData2.rda")
# save(simDataTest2, file = "simDataTest2.rda")

usethis::use_data(simData2, overwrite = TRUE)
usethis::use_data(simDataTest2, overwrite = TRUE)
