###
# We use the SICS algorithm from https://github.com/laxuntus/SICS
###

library(elasticnet)

# Function for fixing signs columnwise
fix_signs <- function(B){
  K <- ncol(B)
  p <- nrow(B)
  
  for(i in 1:K){
    my_col <- B[, i]
    nonzero_ind <- which(abs(my_col) >= 1e-4)[1]
    B[, i] <- my_col*sign(my_col[nonzero_ind])
  }
  
  B
}

# SICS is the main algorithm for sparse invariant coordinate selection
# input:
#   S1, S2 are (scatter) matrices
#   K is the number of extracted invariant coordinates
#   varnum_vec is a vector of non-zero components for each coordinate
#or a number corresponding to same number for every coordinate
# output:
#   B is the ICS loading matrix

SICS = function(S1, S2, K = ncol(S1), varnum_vec = NULL, maximiter = 500){
  p = ncol(S1)
  
  if(length(varnum_vec) == 1){
    varnum_vec <- rep(varnum_vec, K)
  }
  if(is.null(varnum_vec)){
    varnum_vec <- rep(p, K)
  }
  
  S1_inv = solve(S1)
  eig_S1 = eigen(S1, symmetric = TRUE)
  eig_S1_values = ifelse(eig_S1$values<0, 0, eig_S1$values)
  S1_sqrt = eig_S1$vectors %*% diag(sqrt(eig_S1_values)) %*% t(eig_S1$vectors)
  S1_invsqrt <- solve(S1_sqrt)
  
  eig_S2 = eigen(S2, symmetric = TRUE)
  eig_S2_values = ifelse(eig_S2$values<0, 0, eig_S2$values)
  R = eig_S2$vectors %*% diag(sqrt(eig_S2_values)) %*% t(eig_S2$vectors)
  
  A <- eigen(S1_invsqrt %*% S2 %*% S1_invsqrt)$vectors[, 1:K, drop = FALSE]
  B <- S1_invsqrt %*% A
  B <- fix_signs(B)
  
  B_temp <- matrix(0, p, K)
  
  iterat = 0
  while (sum((B_temp-B)^2) > 1e-12) {
    iterat = iterat+1
    
    B_temp = B
    
    for (i in 1:K) {
      B[, i] <- solvebeta(R, R%*%S1_invsqrt%*%A[, i], paras = c(1e-8, varnum_vec[i]), sparse = "varnum")
    }
    B <- fix_signs(B)
    
    svd_A = svd(S1_invsqrt%*%S2%*%B)
    A_0 = svd_A$u %*% t(svd_A$v)
    
    O = eigen(t(A_0)%*%S1_invsqrt%*%S2%*%S1_invsqrt%*%A_0)$vectors
    A = A_0%*%O
    
    if(iterat == maximiter){
      break
    }
  }
  
  col_norms = apply(B, 2, function(x) sqrt(sum(x^2)))
  B =  sweep(B, 2, col_norms, "/")
  
  return(list(A=A, B=B, iters=iterat))
}

###
# MTPCA
###

# Following functions SigmaMatrix and mean_vec handle computing the formulas in MTPCA

# SigmaMatrix computes the latent covariance matrix
SigmaMatrix = function(df, types, sigma_norm){
  
  formula_B2 <- function(x) {qnorm(mean(x))^2}
  formula_E2 <- function(x) {-log(2) + log(mean(x^2)) - 2*log(mean(x))}
  formula_N2 <- function(x) {mean(x^2) - mean(x)^2 - sigma_norm^2}
  formula_P2 <- function(x) {log(mean(x^2) - mean(x)) - 2*log(mean(x))}
  formula_BE <- function(b, e) {sqrt(1+qnorm(mean(b))^2)*(qnorm(mean(b))-qnorm(mean(e*b)/mean(e)))}
  formula_BN <- function(b, n) {(mean(n*b)-mean(n)*mean(b))*sqrt(1+qnorm(mean(b))^2)/dnorm(qnorm(mean(b)))}
  formula_BP <- function(b, p) {sqrt(1+qnorm(mean(b))^2)*(qnorm(mean(p*b)/mean(p))-qnorm(mean(b)))}
  formula_EE <- function(x, y) {log(mean(x*y))-log(mean(x))-log(mean(y))}
  formula_EN <- function(e, n) {(mean(n)*mean(e)-mean(n*e))/mean(e)}
  formula_EP <- function(e, p) {log(mean(p)) + log(mean(e)) - log(mean(p*e))}
  formula_NN <- function(x, y) {mean(x*y)-mean(x)*mean(y)}
  formula_NP <- function(n, p) {(mean(n*p)-mean(n)*mean(p))/mean(p)}
  formula_PP <- function(x, y) {log(mean(x*y))-log(mean(x))-log(mean(y))}
  
  formula_BB = function(x,y){
    sigma_x <- qnorm(mean(x))
    sigma_y <- qnorm(mean(y))
    BB_cov <- function(z){
      arg <- z/(sqrt(1 + sigma_x^2)*sqrt(1 + sigma_y^2))
      pmvnorm(upper = c(sigma_x, sigma_y), mean = c(0, 0), sigma = matrix(c(1, arg, arg, 1), 2, 2)) - mean(x*y)
    }
    limit_1 <- -1*abs(sigma_x*sigma_y)
    limit_2 <- abs(sigma_x*sigma_y)
    if(sign(BB_cov(limit_1)) == sign(BB_cov(limit_2))){
      if(abs(BB_cov(limit_1)) < abs(BB_cov(limit_2))){
        return(limit_1)
      } else {
        return(limit_1)
      }
    } else {
      return(uniroot(BB_cov, interval = c(-1, 1)*abs(sigma_x*sigma_y))$root)
    }
  }
  
  pair_fun <- list(
    BB = formula_BB,
    BE = formula_BE,
    BN = formula_BN,
    BP = formula_BP,
    EE = formula_EE,
    EN = formula_EN,
    EP = formula_EP,
    NN = formula_NN,
    NP = formula_NP,
    PP = formula_PP)
  diag_fun <- list(
    B = formula_B2,
    E = formula_E2,
    N = formula_N2,
    P = formula_P2)
  
  p <- ncol(df)
  result_mat <- matrix(NA, p, p)
  
  for (i in seq_len(p)) {
    for (j in seq_len(i)) {
      type_i <- types[i]
      type_j <- types[j]
      var_i <- df[,i]
      var_j <- df[,j]
      
      if (i == j) {
        fun <- diag_fun[[type_i]]
        result_mat[i, j] <- fun(var_i)
      } else {
        ty <- sort(c(type_i, type_j))
        key <- paste0(ty, collapse = "")
        fun <- pair_fun[[key]]
        if(ty[1] == type_i){result_mat[i, j] <- fun(var_i, var_j)}
        else{result_mat[i, j] <- fun(var_j, var_i)}
        result_mat[j, i] = result_mat[i, j]
      }
    }
  }
  result_mat
}

# mean_vec computes the latent mean vector
mean_vec = function(df, types){
  p = length(types)
  mu = numeric(p)
  for (j in 1:p) {
    if(types[j] == "B"){
      temp = qnorm(mean(df[,j]))
      mu[j] = temp*sqrt(1+temp^2)
    }
    if(types[j] == "E"){
      mu[j] = -0.5*log(2)-2*log(mean(df[,j]))+0.5*log(mean(df[,j]^2))
    }
    if(types[j] == "N"){
      mu[j] = mean(df[,j])
    }
    if(types[j] == "P"){
      mu[j] = 2*log(mean(df[,j]))-0.5*log(mean(df[,j]^2)-mean(df[,j])^2)
    }
  }
  return(mu)
}

# MTPCA is the main function of Mixed-type PCA
# input:
#   df is a data.frame or a n*p matrix 
#   types is a vector giving the types of the variables
#     normal: "N", exponential "E", Poisson "P" and Bernoulli "B"
#   d is the number of components to be extracted
#   varnum is the number of non-zero coefficient
#   sigma_norm is the variance parameter for normal variables
# output:
#   U is the loading matrix
#   lambda is a vector of the component variances
#   Sigma is the latent covariance matrix
#   mu is a vector of the latent means
MTPCA = function(df, types, d=NULL, varnum=ncol(df), sigma_norm=NULL){
  
  p = ncol(df)
  
  norms = which(types=="N")
  no_norms = length(norms)
  if(is.null(d)){
    if(no_norms==0) {d=p}
    else if (!is.null(sigma_norm)){d=p}
    else if (no_norms > 1) {d=no_norms-1}
    else {d=p}
  }
  
  if(is.null(sigma_norm) & d >= no_norms & no_norms != 0) {
    print("EI!")
    return()
  }
  
  if(no_norms >= 2 & is.null(sigma_norm)){
    C = cov(df[,norms])
    sigma_norm = mean(eigen(C)$values[(d+1):no_norms])
  }
  
  Sigma = SigmaMatrix(df, types, sigma_norm)
  
  mu = mean_vec(df, types)
  
  if(varnum<p) {U = SICS(diag(p), Sigma, d, varnum)$B}
  else {U = eigen(Sigma)$vectors[,1:d]}
  
  return(list(U = U, lambda = eigen(Sigma)$values[1:d], Sigma = Sigma, mu = mu))
}

###
# Estimating the component scores
###

library(mvtnorm)

# log_likeli computes the conditional log-likelihood of the component values z
# zi is the z-vector of dimension d
# xi is the corresponding x-vector of dimension p
# ind is vector of the form c("N", "N", "E", "P", "P", "P", "B") of length p
# mu is the mu-vector of length p
# Lambda_vec is the diagonal of the eigenvalue matrix Lambda (length d)
# U is orthonormal matrix of loadings, size p x d
# s2_vec is sigma^2's for the Gaussian variables, assumed to have the same order as for the Gaussian variables in xi and ind
log_likeli <- function(zi, xi, ind, mu, Lambda_vec, U, s2_vec){
  
  temp <- -0.5*sum(Lambda_vec*zi^2)
  
  index_set <- which(ind == "N")
  
  if(length(index_set) > 0){
    for(j in 1:length(index_set)){
      my_index <- index_set[j]
      my_s2 <- s2_vec[j]
      
      temp <- temp - 0.5*(1/my_s2)*(xi[my_index] - mu[my_index] - sum(c(U[my_index, ])*zi))^2
    }
  }
  
  index_set <- which(ind == "E")
  
  if(length(index_set) > 0){
    for(j in 1:length(index_set)){
      my_index <- index_set[j]
      
      temp <- temp - (xi[my_index]*exp(mu[my_index] + sum(c(U[my_index, ])*zi)) - mu[my_index] - sum(c(U[my_index, ])*zi))
    }
  }

  index_set <- which(ind == "P")
  
  if(length(index_set) > 0){
    for(j in 1:length(index_set)){
      my_index <- index_set[j]
      
      temp <- temp - (exp(mu[my_index] + sum(c(U[my_index, ])*zi)) - xi[my_index]*(mu[my_index] + sum(c(U[my_index, ])*zi)))
    }
  }
  
  index_set <- which(ind == "B")
  
  if(length(index_set) > 0){
    for(j in 1:length(index_set)){
      my_index <- index_set[j]
      
      if(xi[my_index] == 0){
        temp <- temp + log(1 - pnorm(mu[my_index] + sum(c(U[my_index, ])*zi)))
      } else {
        temp <- temp + log(pnorm(mu[my_index] + sum(c(U[my_index, ])*zi)))
      }
    }
  }
  
  temp
}

# gradient is the gradient of the log-likeli
gradient <- function(zi, xi, ind, mu, Lambda_vec, U, s2_vec){
  
  temp <- -Lambda_vec*zi

  index_set <- which(ind == "N")
  
  if(length(index_set) > 0){
    for(j in 1:length(index_set)){
      my_index <- index_set[j]
      my_s2 <- s2_vec[j]
      
      temp <- temp + (1/my_s2)*(xi[my_index] - mu[my_index] - sum(c(U[my_index, ])*zi))*U[my_index, ]
    }
  }

  index_set <- which(ind == "E")
  
  if(length(index_set) > 0){
    for(j in 1:length(index_set)){
      my_index <- index_set[j]
      
      temp <- temp - (xi[my_index]*exp(mu[my_index] + sum(c(U[my_index, ])*zi))*U[my_index, ] - U[my_index, ])
    }
  }

  index_set <- which(ind == "P")
  
  if(length(index_set) > 0){
    for(j in 1:length(index_set)){
      my_index <- index_set[j]
      
      temp <- temp - (exp(mu[my_index] + sum(c(U[my_index, ])*zi))*U[my_index, ] - xi[my_index]*U[my_index, ])
    }
  }

  index_set <- which(ind == "B")
  
  if(length(index_set) > 0){
    for(j in 1:length(index_set)){
      my_index <- index_set[j]
      
      if(xi[my_index] == 0){
        temp <- temp - (dnorm(mu[my_index] + sum(c(U[my_index, ])*zi))/(1 - pnorm(mu[my_index] + sum(c(U[my_index, ])*zi))))*U[my_index, ]
      } else {
        temp <- temp + (dnorm(mu[my_index] + sum(c(U[my_index, ])*zi))/pnorm(mu[my_index] + sum(c(U[my_index, ])*zi)))*U[my_index, ]
      }
    }
  }
  
  temp
}




# Gradient descent for a single observation xi
gradient_descent_i <- function(xi, ind, mu, Lambda_vec, U, s2_vec, maxiter = 100){
  
  d <- length(Lambda_vec)
  
  z0 <- rep(0, d)
  grad0 <- gradient(z0, xi, ind, mu, Lambda_vec, U, s2_vec)
  
  iter <- 0
  crit <- 1
  maxiter <- 100
  
  # Initial step size
  epsilon <- 0.05
  
  while(iter < maxiter & crit > 1e-6){
    
    z1 = z0 + epsilon*grad0
    
    grad1 <- gradient(z1, xi, ind, mu, Lambda_vec, U, s2_vec)
    
    # Update step size
    epsilon <- abs(sum((z1 - z0)*(grad1 - grad0)))/sum((grad1 - grad0)^2)
    
    crit <- sqrt(sum((z1 - z0)^2))
    z0 <- z1
    grad0 <- grad1
    
    iter <- iter + 1
  }
  
  c(z1, iter)
}


# Gradient descent for all observations
gradient_descent_all <- function(x, ind, mu, Lambda_vec, U, s2_vec, maxiter = 100){
  
  n <- nrow(x)
  d <- length(Lambda_vec)
  
  res <- matrix(0, n, d + 1)
  
  for(i in 1:n){
    res[i, ] <- gradient_descent_i(x[i, ], ind, mu, Lambda_vec, U, s2_vec)
  }
  
  res
}
