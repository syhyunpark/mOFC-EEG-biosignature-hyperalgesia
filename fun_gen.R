## functions for data generation
u.v.gen <- function(seed, K_T, K_F, R){
  set.seed(seed)
  U <-  randortho(K_T)[, 1:R]
  V <-  randortho(K_F)[, 1:R]
  return(list(U = U, V = V))
}


delta.gen <- function(seed, J, R, p){
  
  set.seed(seed)
  Delta <- array(NA, dim = c(J, R, p))
  for(j in 1:J){
    for(r in 1:R){
      Delta[j,r, ] <- runif(p, -3, 3)
    }
  }
  return(Delta)
}


x.gen <- function(seed, n, p){
  set.seed(seed)
  #return(matrix(rep(1, n), ncol = 1))
  
  return(matrix(runif(n*p), ncol = p))
  
}


dat.gen <- function(seed, n, J, TT, FF,  K_T, K_F, U, V, Delta, X, Sigma2s){
  
  set.seed(seed) 
  
  K <- K_T * K_F
  tt <- seq(0, 1, length.out = TT)  # time points 
  ff <- seq(0, 1, length.out = FF) # frequency points
  
  # generate the basis evaluations
  bsMat_tt <- ns(tt,  knots = seq(0, 1, length.out = K_T)[2:(K_T-1)], intercept = TRUE)
  bsMat_ff <- ns(ff,  knots = seq(0, 1, length.out = K_F)[2:(K_F-1)], intercept = TRUE)
  
  O <- kronecker(bsMat_ff, bsMat_tt)
  #tmp <- svd(O)  # 1st way to orthogonalize
  #O_tilde <- tmp$u %*% diag(tmp$d)
  # O_tilde  <- O%*% tmp$v  # 2nd way to orthogonalize
  
  C_t <- svd(bsMat_tt)$v
  C_f <- svd(bsMat_ff)$v
  
  O_tilde <- kronecker(bsMat_ff%*%C_f, bsMat_tt%*%C_t) 
  #O%*%kronecker(C_f, C_t)  # same result as the previous line
  
   
  # generate the basis coefficients 
  A <- array(0, dim = c(n, J, K_T, K_F))
  B <- array(NA, dim = c(n, K_T, K_F))
  C <- array(NA, dim = c(n, J,  K_T, K_F))
  E <- array(NA, dim = c(n, J,  K_T, K_F))
  R <- dim(Delta)[2]  
  
  for(i in 1:n){
    for(j in 1:J){
      for(r in 1:R){
        A[i, j,,] <-   A[i, j,,] + as.numeric(t(Delta[j,r,])%*%t(matrix(X[i,], ncol = ncol(X))))* U[, r]%*%t(V[,r])
      }
    }
  }

  for(i in 1:n){
    B[i,,] <- matrix(rnorm(K, 0, sd = sqrt(Sigma2s$Sigma2_gamma)), K_T, K_F)
  }
  
  for(i in 1:n){
    for(j in 1:J){
    C[i,j,, ] <- matrix(rnorm(K, 0, sd = sqrt(Sigma2s$Sigma2_omega)), K_T, K_F)
    }
  }
  
  for(i in 1:n){
    for(j in 1:J){
      E[i,j,,] <- A[i,j,,] + B[i,,] +  C[i,j,,]
    }
  }
  
  # generate the responses 
  Y <- Y_clean <- array(NA, dim = c(n, J, TT*FF))
  for(i in 1:n){
    for(j in 1:J){
      Y_clean[i,j,] <-  O_tilde %*% c(E[i,j,,])
      Y[i,j,] <- Y_clean[i,j,] + rnorm(TT*FF, 0, sd =  sqrt(Sigma2s$Sigma2_epsilon))
    }
  }
  
    
  dat <- list(tt = tt, ff = ff, bsMat_tt= bsMat_tt, bsMat_ff = bsMat_ff, C_t =  C_t, C_f = C_f,  O = O, O_tilde = O_tilde, Y = Y, Y_clean = Y_clean, A = A, B = B, C = C, E= E, X = X, U = U, V = V, Delta = Delta)
  
  return(dat)
}

