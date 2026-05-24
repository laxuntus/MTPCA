library(MASS)
library(ICtest)
library(ggplot2)
library(reshape2)

# coldiff calculates the Frobenius distance of two matrices,
# taking account the possibility of flipping signs of columns
coldiff = function(U_true, U_est){
  p = ncol(U_est)
  res = 0
  for (j in 1:p) {
    res = res + min(sum((U_true[,j]-U_est[,j])^2), sum((U_true[,j]+U_est[,j])^2))
  }
  return(sqrt(res))
}

# parameters used in the simulation
Ns = c(250, 1500) # sample sizes
L = diag(c(5,3)) # covariance of the latent components
M = 2000 # number of repetitions

res = array(dim = c(length(Ns), M, 7))

# The real sparse latent loading matrices

# case q = 5
U_pre = 1/sqrt(5)*matrix(c(1,1,1,0,0,-1,-1,0,
                           1,-1,0,1,-1,0,0,1), 8,2)
q = 5

# case q = 2 (choose this by removing the #-characters)
# U_pre = 1/sqrt(2)*matrix(c(1,0,-1,0,0,0,0,0,
#                  0,-1,0,1,0,0,0,0), 8,2)
# q = 2

# The simulation
for (j in 1:length(Ns)) {
  for (i in 1:M) {
    U = U_pre[sample(1:8, replace=FALSE),] # permuting columns of U
    # calculating the sigma corresponding to the indetifiability constraint
    sigma_bin = (U %*% L %*% t(U))[1,1]+0.1
    mu = c(sqrt(sigma_bin*(1+sigma_bin)), 1, 2, 3, 4, 5, 6, 7)
    Z = mvrnorm(Ns[j],c(0,0),L)
    W = mu + U %*% t(Z)
    X1 = rbinom(Ns[j],1,pnorm(W[1,]))
    X2 = rexp(Ns[j], exp(W[2,]))
    X3 = rpois(Ns[j], exp(W[3,]))
    X4 = rnorm(Ns[j], W[4,], 1)
    X5 = rnorm(Ns[j], W[5,], 1)
    X6 = rnorm(Ns[j], W[6,], 1)
    X7 = rnorm(Ns[j], W[7,], 1)
    X8 = rnorm(Ns[j], W[8,], 1)
    df1 = data.frame(X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, X6=X6, X7=X7, X8=X8)
    
    U_sparse = try(MTPCA(df1, c("B", "E", "P", "N", "N", "N", "N", "N"), 2, q)$U, silent = TRUE)
    if(length(U_sparse)==1){
      res[j,i,] = rep(NA,7)
      next
    }
    
    U_sparse2 = try(MTPCA(df1, c("B", "E", "P", "N", "N", "N", "N", "N"), 2, q+2)$U, silent = TRUE)
    if(length(U_sparse)==1){
      res[j,i,] = rep(NA,7)
      next
    }
    
    U_SICS = try(SICS(diag(8), cor(df1), 2, q)$B, silent = TRUE)
    if(length(U_SICS)==1){
      res[j,i,] = rep(NA,5)
      next
    }
    
    U_SICS2 = try(SICS(diag(8), cor(df1), 2, q+2)$B, silent = TRUE)
    if(length(U_SICS)==1){
      res[j,i,] = rep(NA,5)
      next
    }
    

    U_nonsparse = MTPCA(df1, c("B", "E", "P", "N", "N", "N", "N", "N"), 2, 8)$U
    U_pca = prcomp(scale(df1))$rotation[,1:2]
    U_rand = rorth(8)[,1:2]
    
    res[j,i,1] = coldiff(U, U_sparse)
    res[j,i,2] = coldiff(U, U_sparse2)
    res[j,i,3] = coldiff(U, U_nonsparse)
    res[j,i,4] = coldiff(U, U_pca)
    res[j,i,5] = coldiff(U, U_SICS)
    res[j,i,6] = coldiff(U, U_SICS2)
    res[j,i,7] = coldiff(U, U_rand)

    if(i%%50==0){
      print(c(j,i))
    }
  }
}

# extracting the results from the simulation
results = melt(res)[c(1,3,4)]
colnames(results) = c("n", "method", "error")
results$n = factor(results$n, labels = c("250", "1500"))
results$method = rep(c("MTPCA (s)", "MTPCA (s+)", "MTPCA (ns)", "PCA", "SPCA", "SPCA+", "rand"), each = 2*M)

ggplot(results, aes(x=method, y=error, fill=n)) +
  geom_violin(position="dodge") +
  stat_summary(fun=mean, geom="crossbar", position = position_dodge(width = 0.9)) +
  scale_x_discrete(limits=c("MTPCA (s)", "MTPCA (s+)", "MTPCA (ns)", "PCA", "SPCA", "SPCA+", "rand")) +
  theme_bw() +
  ylab("Mean squared error") +
  xlab("Method") +
  ggtitle("q = 5") +
  theme(plot.title = element_text(hjust = 0.5))
