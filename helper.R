## align order and sign of posterior samples of delta
match.order <- function(delta_s, Delta){
  
  J <- dim(Delta)[1]
  R <- dim(Delta)[2]
  
  vec <- 1:R
  R_perms <- permutations(n = R, r = R, v = 1:R)
  J_perms <- permutations(n = J, r = J, v = 1:J)
  condition_signs <- array(NA, dim = c(nrow(R_perms), nrow(J_perms), R))
  
  eval_res <- matrix(0, nrow(R_perms), nrow(J_perms))
  for(r_idx in 1:nrow(R_perms)){
    for(j_idx in 1:nrow(J_perms)){
      R_order <- R_perms[r_idx,]
      J_order <- J_perms[j_idx,]
      delta_s_tmp <- delta_s[, J_order, R_order]
      for(r in 1:R){ # for each rank, the signs for different conditions need to be adjusted together 
        min_idx <- which.min(c(sum(sweep(delta_s_tmp[,,r], 2,Delta[,r,], FUN = "-")^2), sum(sweep(-delta_s_tmp[,,r], 2,Delta[,r,], FUN = "-")^2)))
        eval_res[r_idx, j_idx] <- eval_res[r_idx, j_idx] +  min(c(sum(sweep(delta_s_tmp[,,r], 2,Delta[,r,], FUN = "-")^2), sum(sweep(-delta_s_tmp[,,r], 2,Delta[,r,], FUN = "-")^2)))
        if(min_idx == 2){
          condition_signs[r_idx, j_idx, r] <- -1 # sign needs to be adjusted 
        }else{
          condition_signs[r_idx, j_idx, r] <- 1
        }
      }
    }
  }
  
  indices <- which(eval_res== min(eval_res), arr.ind = TRUE)
  condition_sign <- condition_signs[indices[1], indices[2],]
  R_order <- R_perms[indices[1],]
  J_order <- J_perms[indices[2],]
  
  return(list(R_order = R_order, J_order = J_order,condition_sign = condition_sign ))
  
}



match.order.p <- function(delta_s, Delta){
  
  J <- dim(Delta)[1]
  R <- dim(Delta)[2]
  
  vec <- 1:R
  R_perms <- permutations(n = R, r = R, v = 1:R)
  J_perms <- permutations(n = J, r = J, v = 1:J)
  condition_signs <- array(NA, dim = c(nrow(R_perms), nrow(J_perms), R))
  
  eval_res <- matrix(0, nrow(R_perms), nrow(J_perms))
  for(r_idx in 1:nrow(R_perms)){
    for(j_idx in 1:nrow(J_perms)){
      R_order <- R_perms[r_idx,]
      J_order <- J_perms[j_idx,]
      delta_s_tmp <- delta_s[, J_order, R_order, ]
      for(r in 1:R){ # for each rank, the signs for different conditions need to be adjusted together 
        min_idx <- which.min(c(sum(sweep(delta_s_tmp[,,r, ], c(2,3),Delta[,r,], FUN = "-")^2), sum(sweep(-delta_s_tmp[,,r, ], c(2,3),Delta[,r,], FUN = "-")^2)))
        eval_res[r_idx, j_idx] <- eval_res[r_idx, j_idx] +  min(c(sum(sweep(delta_s_tmp[,,r,], c(2,3),Delta[,r,], FUN = "-")^2), sum(sweep(-delta_s_tmp[,,r,], c(2,3),Delta[,r,], FUN = "-")^2)))
        if(min_idx == 2){
          condition_signs[r_idx, j_idx, r] <- -1 # sign needs to be adjusted 
        }else{
          condition_signs[r_idx, j_idx, r] <- 1
        }
      }
    }
  }
  
  indices <- which(eval_res== min(eval_res), arr.ind = TRUE)
  condition_sign <- condition_signs[indices[1], indices[2],]
  R_order <- R_perms[indices[1],]
  J_order <- J_perms[indices[2],]
  
  return(list(R_order = R_order, J_order = J_order,condition_sign = condition_sign ))
  
}


prepocess <- function(TT, FF, K_T, K_F){
  
  # preprocess the data 
  K <- K_T * K_F
  tt <- seq(0, 1, length.out = TT)  # time points 
  ff <- seq(0, 1, length.out = FF) # frequency points
  
  # generate the basis evaluations
  bsMat_tt <- ns(tt,  knots = seq(0, 1, length.out = K_T)[2:(K_T-1)], intercept = TRUE)
  bsMat_ff <- ns(ff,  knots = seq(0, 1, length.out = K_F)[2:(K_F-1)], intercept = TRUE)
  
  O <- kronecker(bsMat_ff, bsMat_tt)
  
  C_t <- svd(bsMat_tt)$v
  C_f <- svd(bsMat_ff)$v
  
  O_tilde <- kronecker(bsMat_ff%*%C_f, bsMat_tt%*%C_t) 
  
  return(list(O = O,  bsMat_tt = bsMat_tt, bsMat_ff = bsMat_ff, O_tilde = O_tilde, C_t = C_t, C_f = C_f))
  
}


generate_Z0 <- function(n, J) {
  Z0 <- do.call(rbind, replicate(n, diag(J), simplify = FALSE))
  return(Z0)
}

generate_Z1 <- function(n, J) {
  Z1 <- matrix(0, n * J, n)
  for (i in 1:n) {
    Z1[((i-1)*J + 1):(i*J), i] <- 1
  }
  return(Z1)
}

# extract block from matrices 
extract_block <- function(M, block_row, block_col, block_size) {
  row_indices <- ((block_row - 1) * block_size + 1):(block_row * block_size)
  col_indices <- ((block_col - 1) * block_size + 1):(block_col * block_size)
  return(M[row_indices, col_indices])
}

extract_block_vec <- function(qq, block_idx, block_size){
  block_matrix <- matrix(qq, ncol = block_size, byrow = TRUE)
  return(block_matrix[block_idx, ])
  
}


init_params <- function(inits,  K_T, K_F, J, R, n, p){

  set.seed(123)
  
  if(is.null(inits$Sigma2_epsilon)){
    Sigma2_epsilon <- 0.02
  }else{
    Sigma2_epsilon <- inits$Sigma2_epsilon
  }
  
  if(is.null(inits$Sigma2_gamma)){
    Sigma2_gamma <- 0.2
  }else{
    Sigma2_gamma <- inits$Sigma2_gamma
  }
  
  if(is.null(inits$Sigma2_omega)){
    Sigma2_omega <- rep(0.1, n)
  }else{
    Sigma2_omega <- inits$Sigma2_omega
  }
  
  if(is.null(inits$u)){
    u <- t(randortho(K_T)[, 1:R])
  }else{
    u <- inits$u
  }
  
  if(is.null(inits$v)){
    v <- t(randortho(K_F)[, 1:R])
  }else{
    v <-inits$v
  }
  
  if(is.null(inits$delta)){
    delta <- array(runif(R*J*p, -0.1, 0.1) , dim = c(J, R, p))
  } else{
    delta <- inits$delta
  }
  
  return(list(Sigma2_epsilon = Sigma2_epsilon,  Sigma2_gamma = Sigma2_gamma,  Sigma2_omega = Sigma2_omega, u = u, v = v, delta = delta))
}


set_params <- function(params){
  
  if(is.null(params$Sigma2_theta)){
    Sigma2_theta <- 1
  }else{
    Sigma2_theta <- params$Sigma2_theta
  }
  
  if(is.null(params$Sigma2_eta)){
    Sigma2_eta <- 1
  }else{
    Sigma2_eta <- params$Sigma2_eta
  }
  
  if(is.null(params$Sigma2_delta)){
    Sigma2_delta <- rep(10, p)
  }else{
    Sigma2_delta  <- params$Sigma2_delta
  }
  
  if(is.null(params$aa_gamma)){ # parameter for gamma dist
    aa_gamma <- 1
  }else{
    aa_gamma <- params$aa_gamma
  }
  
  if(is.null(params$bb_gamma)){ # parameter for gamma dist
    bb_gamma <- 1
  }else{
    bb_gamma <- params$bb_gamma
  }
  
  if(is.null(params$aa_omega)){ # parameter for gamma dist
    aa_omega <- 1
  }else{
    aa_omega <- params$aa_omega
  }
  
  if(is.null(params$bb_omega)){ # parameter for gamma dist
    bb_omega <- 1
  }else{
    bb_omega <- params$bb_omega
  }
  
  if(is.null(params$aa_beta)){ # parameter for gamma dist
    aa_beta <- 1
  }else{
    aa_beta <- params$aa_beta
  }
  
  if(is.null(params$bb_beta)){ # parameter for gamma dist
    bb_beta <- 1
  }else{
    bb_beta <- params$bb_beta
  }
  
  if(is.null(params$h0)){ 
    h0 <- 0.01
  }else{
    h0 <- params$h0
  }
  
  if(is.null(params$h1)){ 
    h1 <- 1
  }else{
    h1 <- params$h1
  }

  if(is.null(params$h_laplace)){ 
    h_laplace <- 1
  }else{
    h_laplace <- params$h_laplace
  }
  
  
  return(list(Sigma2_theta = Sigma2_theta, Sigma2_eta= Sigma2_eta, Sigma2_delta = Sigma2_delta, 
              aa_gamma = aa_gamma, bb_gamma = bb_gamma, aa_omega = aa_omega, bb_omega = bb_omega,
              aa_beta = aa_beta, bb_beta = bb_beta, h0 = h0, h1 = h1, h_laplace = h_laplace))
}



