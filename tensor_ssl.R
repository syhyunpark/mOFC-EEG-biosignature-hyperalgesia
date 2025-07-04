## revise the code to allow missing conditions
## first revise the sampling of eta 

tensor.decomp.ssl <- function(seed, Y, X, JJ,  TT, FF, R, K_T, K_F, n_sample, n_burn, params = NULL, inits = NULL, save_all = FALSE, type_simple = TRUE, threshold = 0.0001){
  
  # control the randomness of sampling 
  set.seed(seed) 

  # extract dimensions
  n <- dim(Y)[1]
  J <- dim(Y)[2] # maximum number of conditions
  p <- dim(X)[2]
  K <- K_T * K_F
  
  
  ## if no J_missing, then do the fast algorithm 
  idx_no_missing <- 1:n
  if(sum(JJ) == n*J){
    J_missing <- FALSE
  }else{
    J_missing <- TRUE
    idx_missing <- which(apply(JJ, 1, sum) <J)
    idx_no_missing <- setdiff(1:n, idx_missing)
  }
  
  J_n <- sum(JJ)
  

  R_num <- 1 # start from the simplest model
  
  K_T_tilde <- K_T - R_num + 1
  K_F_tilde <- K_F - R_num + 1
  
  # preprocess the data (assume evenly spaced in the time and frequency domain)
  tmp_prep <- prepocess(TT, FF, K_T, K_F)
  O <- tmp_prep$O
  O_tilde <- tmp_prep$O_tilde
  d_ls <- diag(t(O_tilde)%*%O_tilde) # a vector of length K
  C_t <- tmp_prep$C_t
  C_f <- tmp_prep$C_f
  
  # initialize objects to save samples 
  n_total <- n_sample+n_burn # total number of iterations
  
  # number of saved iterations 
  if(save_all){
    n_save <- n_total
  }else{
    n_save <- n_sample
  }
  
  # --- save posterior samples
  omega_s <- array(NA, dim = c(n_save, K, J, n)) # for each l = 1,..., K, one value per (i, j) combination
  gamma_s <- array(NA, dim = c(n_save, K, n)) # for each l = 1,..., K, one value per i 
  alpha_s <- array(NA, dim = c(n_save, K, J, n)) # for each l = 1,..., K, one value per j 
  u_s <- array(NA, dim = c(n_save, R, K_T))
  v_s <- array(NA, dim = c(n_save, R, K_F))
  delta_s <- array(NA, dim = c(n_save,J, R, p))
  Sigma2_epsilon_s <- rep(NA, n_save)
  Sigma2_gamma_s <- rep(NA, n_save)
  Sigma2_omega_s <- matrix(NA, n_save, n)
  tau_s <- matrix(0, n_total, R)
  
  
  # --- initialization 
  tmp_init <- init_params(inits,  K_T, K_F, J, R, n, p) 
  Sigma2_epsilon_current <- tmp_init$Sigma2_epsilon
  Sigma2_gamma_current <- tmp_init$Sigma2_gamma
  Sigma2_omega_current <- tmp_init$Sigma2_omega
  u_current <- tmp_init$u
  v_current <- tmp_init$v
  delta_current <- tmp_init$delta
  pi_current <- 0.5
  tau_current <- rep(NA, R)
  tau_current[1] <- 1
  
  # --- extract hyper-parameters
  tmp_params <- set_params(params) 
  Sigma2_theta <- tmp_params$Sigma2_theta
  Sigma2_eta <- tmp_params$Sigma2_eta
  Sigma2_delta <- tmp_params$Sigma2_delta
  aa_gamma <- tmp_params$aa_gamma # for the Inverse Gamma distribution for the variance of gamma 
  bb_gamma <- tmp_params$bb_gamma  # for the Inverse Gamma distribution for variance of gamma 
  aa_omega <- tmp_params$aa_omega # for the Inverse Gamma distribution for the variance of omega
  bb_omega <- tmp_params$bb_omega  # for the Inverse Gamma distribution for variance of omega
  aa_beta <- tmp_params$aa_beta # for the beta distribution for pi
  bb_beta <- tmp_params$bb_beta # for the beta distribution for pi
  h0 <- tmp_params$h0 # for the mixture Laplace  distribution
  h1 <- tmp_params$h1 # for the mixture Laplace  distribution
  
  
  # --- orthogonalize basis evaluations
  Y_tilde <- array(NA, dim = c(n, J, K))  # this is of dimension n, J, K (corresponding to \tilde{y}_{i,j}'s)
  O_tilde_tmp <- diag(1/d_ls) %*% t(O_tilde)
  
  Y_tilde <- aperm(apply(Y, c(1, 2), function(vec) O_tilde_tmp %*% vec), c(2, 3, 1))
  
  #for(i in 1:n){
  #  for(j in 1:J){
  #    Y_tilde[i,j,] <- O_tilde_tmp%*%Y[i,j,]
  #  }
  #}
  
  Z <- generate_Z1(n, J) # an (nJ) \times n matrix 
  
  UV_patch_vec <- matrix(NA, R, K) # need to be updated
  
  for(r in 1:R){
    UV_patch_vec[r, ] <- kronecker(v_current[r, ], u_current[r, ]) # form the previous iteration
  }
  
  Y_ss <- sum(apply(Y, c(1, 2), function(x) t(x) %*% x), na.rm = TRUE) # used to fast compute the residuals, removed missing 
  
  
  ff_theta_eta <- function(Q, g){
    Q_inv <- solve(Q)
    tmp <- mvrnorm(1, mu = Q_inv%*%g, Sigma = Q_inv)
    return(tmp/sqrt(sum(tmp^2)))
  }
  
  remove_attmpt <- 0  #once have attempted removal, then no more adding attempt
  R_indices <- 1
  R_num <- length(R_indices)
  
  for(s in 1:n_total){
    
    
    if(s%%20 == 0){
      print(paste(s, "-th iteration", sep = ""))
      print(paste( "Active ranks: ", R_indices, sep = ""))
      print("Tau_r")
      print(tau_current) 
    }
    
    
    if(s%%100 == 0){
      
      if(length(R_indices)>1){
        rank_remove <- R_indices[which(apply(tau_s[(s-100+1):s, R_indices], 2, mean) < threshold)]
      }else{
        rank_remove <- R_indices[mean(tau_s[(s-100+1):s, R_indices]) < threshold]
      }
      
      
      if(length(rank_remove) >= 1){
        remove_attmpt <- 1
        R_indices <- setdiff(R_indices, rank_remove)
        R_num <- length(R_indices)
        K_T_tilde <- K_T - R_num + 1
        K_F_tilde <- K_F - R_num + 1
      }else{
        if(remove_attmpt == 0){ 
          if(length(R_indices)<R){
            
            u_current[setdiff(1:R, R_indices),] <- t(null_complement(matrix(t(u_current[R_indices, ]), ncol = R_num), universe = NULL, na.allow = TRUE))[1:length(setdiff(1:R, R_indices)),]
            v_current[setdiff(1:R, R_indices),] <- t(null_complement(matrix(t(v_current[R_indices, ]), ncol = R_num), universe = NULL, na.allow = TRUE))[1:length(setdiff(1:R, R_indices)),]
            
            tau_current[setdiff(1:R, R_indices)[1]] <- 1
            R_indices <- c(R_indices, setdiff(1:R, R_indices)[1]) # add one rank
            
            R_num <- length(R_indices)
            K_T_tilde <- K_T - R_num + 1
            K_F_tilde <- K_F - R_num + 1
            
          }else{
            print("R might be too small!")
          }
        }
      }
    }
    
    
    # Block 1
    # Block 1.1: sample the fixed effects  
    # Step (a): 
    if(type_simple){ # homogeneous variance
      
      Sigma_inv <- Sigma_inv_tilde <- matrix(0, nrow = J*K, ncol = J*K)
      Sigma2_epsilon_rep <- rep(Sigma2_epsilon_current/d_ls, J)
      Sigma2_epsilon_inv_rep <- rep(Sigma2_epsilon_current/c(t(matrix(d_ls,nrow = K_T,ncol = K_F))), J)
        

      i_inv_diagonal <- 1/(Sigma2_omega_current[1] + Sigma2_epsilon_rep) 
      scaling_fac <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_inv_diagonal))
      Sigma_inv <- diag(i_inv_diagonal) - scaling_fac *  outer(i_inv_diagonal, i_inv_diagonal, "*")
      
      i_tilde_inv_diagonal <- 1/(Sigma2_omega_current[1] +  Sigma2_epsilon_inv_rep)
      scaling_tilde_fac <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_tilde_inv_diagonal))
      Sigma_inv_tilde <- diag(i_tilde_inv_diagonal)  - scaling_tilde_fac* outer(i_tilde_inv_diagonal,i_tilde_inv_diagonal, "*")
      

    }else{
      
      Sigma_inv_i_s <- Sigma_inv_i_tilde_s <- lapply(1:n, function(x) matrix(0, nrow = J*K, ncol = J*K))
      Sigma2_epsilon_rep <- rep(Sigma2_epsilon_current/d_ls, J)
      #Sigma2_epsilon_inv_rep <- rep(Sigma2_epsilon_current/d_ls, each = J)
      Sigma2_epsilon_inv_rep <- rep(Sigma2_epsilon_current/c(t(matrix(d_ls, nrow = K_T,ncol = K_F))), J)
      
      for(i in 1:n){
        
        i_inv_diagonal <- 1/(Sigma2_omega_current[i] + Sigma2_epsilon_rep) 
        scaling_fac <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_inv_diagonal))
        Sigma_inv_i_s[[i]] <-  diag(i_inv_diagonal) - scaling_fac *  outer(i_inv_diagonal, i_inv_diagonal, "*")

        i_tilde_inv_diagonal <- 1/(Sigma2_omega_current[i] +  Sigma2_epsilon_inv_rep)
        scaling_tilde_fac <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_tilde_inv_diagonal))
        Sigma_inv_i_tilde_s[[i]] <- diag(i_tilde_inv_diagonal)  - scaling_tilde_fac* outer(i_tilde_inv_diagonal,i_tilde_inv_diagonal, "*")
        
      }
    }
    
    if(J_missing){ # compute the missing covariates
      Sigma_inv_i_s_missing <- Sigma_inv_i_tilde_s_missing <- lapply(1:length(idx_missing), function(x) matrix(0, nrow = J*K, ncol = J*K))
      
      for(qq in 1:length(idx_missing)){ # go through all missing individuals 
        
        J_indices <- which(JJ[idx_missing[qq],]!=0) # index for non-missing conditions
        existing_cov_indices <- NULL
        
        for(kk in J_indices){
          existing_cov_indices  <- c(existing_cov_indices, ((kk-1)*K + 1): (kk*K))
        }
        
        Sigma2_epsilon_rep_tmp <- rep(Sigma2_epsilon_current/d_ls, length(J_indices))
        
        i_inv_diagonal_tmp <- 1/(Sigma2_omega_current[idx_missing[qq]] + Sigma2_epsilon_rep_tmp) 
        scaling_fac_tmp <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_inv_diagonal_tmp))
        Sigma_inv_tmp <- diag(i_inv_diagonal_tmp) - scaling_fac_tmp *  outer(i_inv_diagonal_tmp, i_inv_diagonal_tmp, "*")
        Sigma_inv_i_s_missing[[qq]][existing_cov_indices, existing_cov_indices] <- Sigma_inv_tmp
        
        Sigma2_epsilon_inv_rep_tmp <- rep(Sigma2_epsilon_current/c(t(matrix(d_ls,nrow = K_T,ncol = K_F))), length(J_indices))
        i_tilde_inv_diagonal_tmp <- 1/(Sigma2_omega_current[idx_missing[qq]] +  Sigma2_epsilon_inv_rep_tmp)
        scaling_tilde_fac_tmp <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_tilde_inv_diagonal_tmp))
        Sigma_inv_i_tilde_tmp <- diag(i_tilde_inv_diagonal_tmp)  - scaling_tilde_fac_tmp* outer(i_tilde_inv_diagonal_tmp,i_tilde_inv_diagonal_tmp, "*")
        Sigma_inv_i_tilde_s_missing[[qq]][existing_cov_indices, existing_cov_indices] <- Sigma_inv_i_tilde_tmp
        
      }
      
    }
    
    
    
    for(r in R_indices){
      
      if((s == 1) & (r == 1)){
        q_r_j_s <- lapply(1:n, function(x) matrix(0, nrow = J, ncol = K)) # computed for this specific r, initialization
      }
      
      R_indices_neg <- setdiff(R_indices, r)
      
      for(i in 1:n){
        for(j in 1:J){
          if(R_num > 2){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,] -  (t(X[i,]) %*% t(delta_current[j,R_indices_neg, ])) %*% UV_patch_vec[R_indices_neg, ] # n by k 
          }
          if(R_num == 2){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,] -  (t(X[i,]) %*% delta_current[j,R_indices_neg, ]) %*% UV_patch_vec[R_indices_neg, ] # n by k 
          }
          if(R_num == 1){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,] 
          }
        }
      } # q_r_j_s may contain some missing values
      
      # Step a.1 
      if(R_num == 1){
        B_r <- diag(K_T)
      }else{
        B_r <- null_complement(matrix(t(u_current[R_indices_neg, ]), ncol = R_num-1), universe = NULL, na.allow = TRUE) # need to double check for larger R, can be replaced for Gram-Schmidt later
      }
      
      lambda_tilde_r <- X%*% t(delta_current[ ,r, ]) # dim = c(n, J)
      v_tilde_r <- t(apply(lambda_tilde_r, 1, function(row) kronecker(row, v_current[r, ])))
      
      H_r <- array(0, c(n, K_T, K_T))
      w_r <- array(0, c(n, K_T))
      
      v_outer <- array(0, c(n, J*K_F, J*K_F))  # Store outer products of v_tilde_r
      
      for (l1 in 1:(J*K_F)){
        for (l2 in 1:(J*K_F)) {
          v_outer[, l1, l2] <- v_tilde_r[, l1] * v_tilde_r[, l2]
        }
      }
      
      for(l2 in 1:(J*K_F)){
        K_F_idx <- (l2 - 1) %% K_F + 1
        J_idx <- (l2 - 1) %/% K_F + 1
        q_indices <-  ((K_F_idx-1)*K_T+1): (K_F_idx*K_T)
        
        for(l1 in 1:(J*K_F)){
          block_row <- l1; block_col = l2; block_size = K_T
          row_indices <- ((block_row - 1) * block_size + 1):(block_row * block_size)
          col_indices <- ((block_col - 1) * block_size + 1):(block_col * block_size) 
          
          if(type_simple){
            inv_block_mat  <-  Sigma_inv[row_indices, col_indices] 
            H_r[idx_no_missing,,] <- H_r[idx_no_missing,,] + array(t(sapply(v_outer[idx_no_missing, l1, l2], FUN = function(a){inv_block_mat*a})), dim = c(length(idx_no_missing), K_T, K_T)) 
            if(J_missing){
              for(qq in 1:length(idx_missing)){
                if(JJ[idx_missing[qq],J_idx]!=0){ # this condition exist
                  inv_block_mat <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices]
                  H_r[idx_missing[qq],,] <- H_r[idx_missing[qq],,] + v_outer[idx_missing[qq], l1, l2]*inv_block_mat 
                }
              }
            }
            for(i in 1:n){
              if(i %in% idx_no_missing){
                  w_r[i, ] <- w_r[i, ] +  v_tilde_r[i, l1] * (inv_block_mat %*% q_r_j_s[[i]][J_idx, q_indices])
              }else{
                if(JJ[i,J_idx]!=0){ # this condition exist
                  qq <- which(idx_missing == i)
                  inv_block_mat <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices]
                  w_r[i, ] <- w_r[i, ] +  v_tilde_r[i, l1] * (inv_block_mat %*% q_r_j_s[[i]][J_idx, q_indices])
                }
              }
            }
          }else{
            for(i in 1:n){
                if(i %in% idx_no_missing){
                  inv_block_mat  <-  Sigma_inv_i_s[[i]][row_indices, col_indices]
                  H_r[i,,] <- H_r[i,,] + v_outer[i, l1, l2] * inv_block_mat
                  w_r[i, ] <- w_r[i, ] +  v_tilde_r[i, l1] * (inv_block_mat %*% q_r_j_s[[i]][J_idx, q_indices])
                }else{
                  if(JJ[i,J_idx]!=0){ # this condition exist
                    qq <- which(idx_missing == i)
                    inv_block_mat <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices]
                    H_r[i,,] <- H_r[i,,] + v_outer[i, l1, l2] * inv_block_mat
                    w_r[i, ] <- w_r[i, ] +  v_tilde_r[i, l1] * (inv_block_mat %*% q_r_j_s[[i]][J_idx, q_indices])
                  }
                }
              }
           }
        }
      }
      
      
      
      #Q_u  <- diag(rep(1/Sigma2_theta, K_T_tilde)) +
      Q_u  <-  t(B_r)%*%apply(H_r, c(2,3), sum)%*%B_r
      g_u  <- t(B_r) %*% apply(w_r, 2, sum)
      theta_r_current <- tryCatch({c(rFisherBingham(1, mu = g_u, Aplus = -Q_u/2, mtop = 10000))}, error = function(e){ff_theta_eta(Q_u, g_u)})
      u_current[r,] <- B_r%*%theta_r_current
      
      
      # Step a.2
      if(R_num == 1){
        S_r <- diag(K_F)
      }else{
        S_r <- null_complement(matrix(t(v_current[R_indices_neg, ]), ncol = R_num-1), universe = NULL, na.allow = TRUE) # can be replaced for Gram-Schmidt later
      }
      
      u_tilde_r <- t(apply(lambda_tilde_r, 1, function(row) kronecker(row, u_current[r, ]))) # dim = c(n, J*K_T)
      
  
      L_r <- array(0, c(n, K_F, K_F))
      h_r <- array(0, c(n, K_F))

      
      u_outer <- array(0, c(n, J*K_T, J*K_T))  # Store outer products of v_tilde_r
      
      for (l1 in 1:(J*K_T)){
        for (l2 in 1:(J*K_T)) {
          u_outer[, l1, l2] <- u_tilde_r[, l1] * u_tilde_r[, l2]
        }
      }
      
      for(l2 in 1:(J*K_T)){
        K_T_idx <- (l2 - 1) %% K_T + 1
        J_idx <- (l2 - 1) %/% K_T + 1
        q_indices <-  ((K_T_idx-1)*K_F+1): (K_T_idx*K_F)
        
        for(l1 in 1:(J*K_T)){
          block_row <- l1; block_col = l2; block_size = K_F
          row_indices <- ((block_row - 1) * block_size + 1):(block_row * block_size)
          col_indices <- ((block_col - 1) * block_size + 1):(block_col * block_size) 
          
          if(type_simple){
            inv_block_mat  <-  Sigma_inv_tilde[row_indices, col_indices] 
            L_r[idx_no_missing,,] <- L_r[idx_no_missing,,] + array(t(sapply(u_outer[idx_no_missing, l1, l2], FUN = function(a){inv_block_mat*a})), dim = c(length(idx_no_missing), K_F, K_F)) 
            
            if(J_missing){
              for(qq in 1:length(idx_missing)){
                if(JJ[idx_missing[qq],J_idx]!=0){ # this condition exist
                  inv_block_mat <- Sigma_inv_i_tilde_s_missing[[qq]][row_indices, col_indices]
                  L_r[idx_missing[qq],,] <- L_r[idx_missing[qq],,] + u_outer[idx_missing[qq], l1, l2]*inv_block_mat 
                }
              }
            }
            for(i in 1:n){
              if(i %in% idx_no_missing){
                h_r[i, ] <- h_r[i, ] +  u_tilde_r[i, l1] * (inv_block_mat %*% c(t(matrix(q_r_j_s[[i]][J_idx, ], nrow = K_T, ncol = K_F)))[q_indices])
              }else{
                if(JJ[i,J_idx]!=0){ # this condition exist
                  qq <- which(idx_missing == i)
                  inv_block_mat <- Sigma_inv_i_tilde_s_missing[[qq]][row_indices, col_indices]
                  h_r[i, ] <- h_r[i, ] +  u_tilde_r[i, l1] * (inv_block_mat %*% c(t(matrix(q_r_j_s[[i]][J_idx, ], nrow = K_T, ncol = K_F)))[q_indices])
                }
              }
            }
          }else{
            
            for(i in 1:n){
              if(i %in% idx_no_missing){
                inv_block_mat  <-  Sigma_inv_i_tilde_s[[i]][row_indices, col_indices]
                L_r[i,,] <- L_r[i,,] + u_outer[i, l1, l2] * inv_block_mat
                h_r[i, ] <- h_r[i, ] +  u_tilde_r[i, l1] * (inv_block_mat %*% c(t(matrix(q_r_j_s[[i]][J_idx, ], nrow = K_T, ncol = K_F)))[q_indices])
              }else{
                if(JJ[i,J_idx]!=0){ # this condition exist
                  qq <- which(idx_missing == i)
                  inv_block_mat <- Sigma_inv_i_tilde_s_missing[[qq]][row_indices, col_indices]
                  L_r[i,,] <- L_r[i,,] + u_outer[i, l1, l2] * inv_block_mat
                  h_r[i, ] <- h_r[i, ] +  u_tilde_r[i, l1] * (inv_block_mat %*% c(t(matrix(q_r_j_s[[i]][J_idx, ], nrow = K_T, ncol = K_F)))[q_indices])
                }
              }
            }
          }
        }
      }
      
      #Q_v  <- diag(rep(1/Sigma2_eta, K_F_tilde)) + 
      Q_v  <-  t(S_r)%*%apply(L_r, c(2,3), sum)%*%S_r
      g_v  <- t(S_r) %*% apply(h_r, 2, sum)
      eta_r_current <- tryCatch({c(rFisherBingham(1, mu = g_v, Aplus = -Q_v/2, mtop = 10000))}, error = function(e){ff_theta_eta(Q_v, g_v)})
      
      v_current[r,] <- S_r%*%eta_r_current
      UV_patch_vec[r, ] <- kronecker(v_current[r, ], u_current[r, ]) # does affect the partial residuals
    }
      
      

    
    # Step (b)
    for(r in R_indices){
      
      R_indices_neg <- setdiff(R_indices, r)
      
      for(i in 1:n){
        for(j in 1:J){
          if(R_num > 2){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,] -  (t(X[i,]) %*% t(delta_current[j,R_indices_neg, ])) %*% UV_patch_vec[R_indices_neg, ] # n by k 
          }
          if(R_num == 2){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,] -  (t(X[i,]) %*% delta_current[j,R_indices_neg, ]) %*% UV_patch_vec[R_indices_neg, ] # n by k 
          }
          if(R_num == 1){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,] 
          }
        }
      }
      
      Delta_r_s <- array(NA, c(n, J, p, p)) # for a specific R
      Chi_r_s <- array(NA, c(n, J, p))
      
      
      for(i in 1:n){
        
        for(j in 1:J){
          block_row <- j; block_col = j; block_size = K
          row_indices <- col_indices <- ((block_row - 1) * block_size + 1):(block_row * block_size)
          
          if(type_simple){
            if(i %in% idx_no_missing){
              inv_block_mat  <- Sigma_inv[row_indices, col_indices] 
            }else{
              qq <- which(idx_missing == i)
              inv_block_mat  <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices] 
            }
          }else{
            if(i %in% idx_no_missing){
              inv_block_mat  <- Sigma_inv_i_s[[i]][row_indices, col_indices] 
            }else{
              qq <- which(idx_missing == i)
              inv_block_mat  <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices] 
            }
          }
          
          if(p == 1){
            Delta_r_s[i,j,,] <- X[i,]^2 * (t(UV_patch_vec[r,])%*%inv_block_mat%*% UV_patch_vec[r,])
          }else{
            Delta_r_s[i,j,,] <- (X[i,]%*%t(X[i,])) * c(t(UV_patch_vec[r,])%*%inv_block_mat%*% UV_patch_vec[r,])
          }
          
          tmp_sum <- rep(0, K)
          for(jj in setdiff(1:J, j)){
            block_row <- j; block_col = jj; block_size = K
            col_indices <- ((block_col - 1) * block_size + 1):(block_col * block_size) 
            if(type_simple){
              if(i %in% idx_no_missing){
                inv_block_mat_tmp  <- Sigma_inv[row_indices, col_indices] 
              }else{
                qq <- which(idx_missing == i)
                inv_block_mat_tmp  <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices] 
              }
            }else{
              if(i %in% idx_no_missing){
                inv_block_mat_tmp  <- Sigma_inv_i_s[[i]][row_indices, col_indices] 
              }else{
                qq <- which(idx_missing == i)
                inv_block_mat_tmp  <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices] 
              }
            }
            if(JJ[i,jj]==1){
              tmp_sum <- tmp_sum + inv_block_mat_tmp %*% (q_r_j_s[[i]][jj, ] - c(t(X[i,]) %*% delta_current[jj,r, ])* UV_patch_vec[r,])
            }
          }
          Chi_r_s[i,j, ] <- X[i,] * c(t(UV_patch_vec[r,]) %*% (inv_block_mat %*% q_r_j_s[[i]][j, ] + tmp_sum))
        }
      }
      
      for(j in 1:J){
        if(p==1){ # only has one covariate 
          Q_delta <- (1/tau_current[r])*(1/Sigma2_delta) + sum(Delta_r_s[,j,,])
          g_delta <- sum(Chi_r_s[,j, ])
          delta_current[j,r,] <- rnorm(1, mean = solve(Q_delta)%*%g_delta, sd = solve(Q_delta))
        }else{
          Q_delta <-  (1/tau_current[r])*diag(1/Sigma2_delta) + apply(Delta_r_s[,j,,], c(2,3), sum)
          g_delta <- apply(Chi_r_s[which(JJ[,j]==1),j, ], 2, sum)
          delta_current[j,r,] <- mvrnorm(1, mu = solve(Q_delta)%*%g_delta, Sigma = solve(Q_delta))
        }
      }
    }
    
    # Step (b.2)
    delta_sigma_inv_delta <- matrix(0, J, R)
    
    if(p >1){
      for(jj in 1:J){
        for(rr in R_indices){
          delta_sigma_inv_delta[jj,rr] <- t(delta_current[jj,rr,])%*%diag(1/Sigma2_delta)%*%delta_current[jj,rr,]
        }
      }
    }else{
      for(jj in 1:J){
        for(rr in R_indices){
          delta_sigma_inv_delta[jj,rr] <- delta_current[jj,rr,]^2*(1/Sigma2_delta)
        }
      }
    }
    
    b_current <- rep(0, R)
    for(r in R_indices){
      
      p_laplace <- -p*J/2+ 1
      b_laplace <-  sum(delta_sigma_inv_delta[,r])
      
      if(b_current[r]==1){ 
        a_laplace <- 2/h0
      }else{
        a_laplace <- 2/h1
      }
      
      tau_current[r] <- (p_laplace-1 + sqrt((p_laplace-1)^2 + a_laplace*b_laplace))/a_laplace
      p_bernouli_tmp_1 <- dlaplace(tau_current[r], 0, h0)*pi_current
      p_bernouli_tmp_2 <- dlaplace(tau_current[r], 0, h1)*(1-pi_current)
      
      p_bernouli <- p_bernouli_tmp_1/(p_bernouli_tmp_1+ p_bernouli_tmp_2)
      b_current[r] <- rbinom(1, 1, p_bernouli)
    }
    
    tau_s[s,] <- tau_current
    
    pi_current <- rbeta(1, aa_beta + sum(b_current), bb_beta + R_num  - sum(b_current))
    
    
    alpha_current <- array(0, dim = c(K, J, n))
    
    for(r in R_indices){
      X_delta <- X %*% t(delta_current[, r, ])  # Resulting in an n x J matrix
      for (k in 1:K) {
        alpha_current[k, , ] <- alpha_current[k, , ] + t(X_delta)*UV_patch_vec[r, k]
      }
    }
    
    # Block 1.2: sample gamma (the random effects for each person)
    gamma_current <- matrix(NA, K, n) 
    for(l in 1:K){
      Sigma_1_vec <- rep(Sigma2_omega_current, each = J) + Sigma2_epsilon_current/d_ls[l]
      qq <- ((c(t(Y_tilde[,,l])) - c(alpha_current[l,,]))/Sigma_1_vec)
      idx_tmp <- which(!is.na(qq))
      g_gamma <- t(Z)[, idx_tmp]%*%qq[idx_tmp]
      Q_gamma_inv_vec <- 1/(1/Sigma2_gamma_current +   diag(t(Z[idx_tmp,]) %*%diag(1/Sigma_1_vec[idx_tmp]) %*%Z[idx_tmp,]))
      gamma_current[l, ] <- rnorm(n, mean = Q_gamma_inv_vec * g_gamma, sd = sqrt(Q_gamma_inv_vec))
    }
    
    
    # Block 1.3: sample omega  # (needs to be revised) 
    D_omega_vec <- rep(Sigma2_omega_current, each = J)
    omega_current <- array(NA, dim = c(K, J, n))
    
    for(l in 1:K){
      g_omega <- (d_ls[l]/Sigma2_epsilon_current)* ((c(t(Y_tilde[,,l])) - c(alpha_current[l,,]) -  Z%*%gamma_current[l,]))
      idx_tmp <- which(is.na(g_omega))
      if(length(idx_tmp)!=0){
        g_omega[idx_tmp] <- 0
      }
      Q_omega_inv_vec <- 1/(1/D_omega_vec + d_ls[l]/Sigma2_epsilon_current)
      omagea_tmp <- rnorm(n * J, mean = Q_omega_inv_vec * g_omega, sd = sqrt(Q_omega_inv_vec))
      if(length(idx_tmp)!=0){
        omagea_tmp[idx_tmp] <- NA
      }
      omega_current[l,,] <- array(omagea_tmp, dim = c(J, n))
    }
    
    
    res_tmp <- matrix(NA, nrow = n, ncol = J) # save the sum of squared residuals for i, j combination (except for the 1st ter, y^Ty)
    for(i in 1:n){
      for(j in 1:J){
        beta_tmp <-alpha_current[,j,i] + gamma_current[,i] + omega_current[,j, i]
        beta_tmp_star <- beta_tmp*d_ls
        res_tmp[i,j] <- -2* t(beta_tmp_star)%*%Y_tilde[i,j,] + t(beta_tmp_star)%*%beta_tmp
      }
    }
    
    # Block 2  
    # Block 2.1 
    #Sigma2_epsilon_current <- 1/rgamma(1, shape = (n*J*TT*FF)/2, rate = sum((Y - Y_pred)^2)/2) 
    Sigma2_epsilon_current <- 1/rgamma(1, shape = (J_n*TT*FF)/2, rate = (sum(res_tmp, na.rm = TRUE) + Y_ss)/2) 
    
    # Block 2.2 
    Sigma2_gamma_current <- 1/rgamma(1, aa_gamma + n*K/2,  rate = bb_gamma + sum(gamma_current^2)/2) 
    
    # Block 2.3
    if(type_simple){
      Sigma2_omega_current <- rep(1 / rgamma(1, shape =  aa_omega + J_n* K / 2, rate =  bb_omega + sum(omega_current^2, na.rm = TRUE) / 2), n)
    }else{
      Sigma2_omega_current <- 1 / rgamma(n, shape =  aa_omega + apply(JJ, 1, sum) * K / 2, rate =  bb_omega + colSums(omega_current^2, dims = 2, na.rm = TRUE) / 2)
    }
    
    # save the posterior samples 
    if(save_all){  
      Sigma2_epsilon_s[s]  <- Sigma2_epsilon_current
      Sigma2_gamma_s[s] <- Sigma2_gamma_current
      Sigma2_omega_s[s,] <-Sigma2_omega_current
      v_s[s,,] <- v_current
      u_s[s,,] <- u_current
      alpha_s[s,,,] <- alpha_current
      delta_s[s,,, ] <- delta_current
      gamma_s[s,,] <- gamma_current
      omega_s[s,,,]<- omega_current
      sample_idx <- (n_burn+1):n_save
      
    }else{
      if(s > n_burn){
        Sigma2_epsilon_s[s - n_burn]  <- Sigma2_epsilon_current
        Sigma2_gamma_s[s - n_burn] <- Sigma2_gamma_current
        Sigma2_omega_s[s - n_burn,] <-Sigma2_omega_current
        v_s[s - n_burn,,] <- v_current
        u_s[s - n_burn,,] <- u_current
        alpha_s[s - n_burn,,,] <- alpha_current
        delta_s[s - n_burn,,, ] <- delta_current
        gamma_s[s - n_burn,,] <- gamma_current
        omega_s[s - n_burn,,,] <- omega_current
        sample_idx <- 1:n_save
      }
    }
  }
  
  Sigma2_s <- list(Sigma2_epsilon_s = Sigma2_epsilon_s, Sigma2_gamma_s = Sigma2_gamma_s,  Sigma2_omega_s  =  Sigma2_omega_s)
  dat2return <- list(R_indices = R_indices, C_t  = C_t , C_f = C_f, tmp_init = tmp_init, alpha_s = alpha_s, omega_s = omega_s, gamma_s = gamma_s,  v_s = v_s, u_s = u_s, delta_s = delta_s, Sigma2_s = Sigma2_s)
  
  return(dat2return)
}

