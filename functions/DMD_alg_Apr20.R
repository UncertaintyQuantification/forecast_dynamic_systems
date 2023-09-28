library(MASS)
#### DMD implementation
### check minimize L2 norm (reality) d=6
# DMD_alg <- function(output, method="ratio_sum", threshold=0.999, 
#                     threshold_eigen_A_tilde=10^{-8}, s=1){
#   ## record rank r
#   n = ncol(output)
#   k = nrow(output)
#   
#   Y_index = seq(2, n, by=s)
#   X = output[,(Y_index-1)]
#   Y = output[,Y_index]
#   SVD_X = svd(X)
#   if(method=="ratio_sum"){
#     sum_d = sum(SVD_X$d)
#     cumsum_d = cumsum(SVD_X$d)
#     r = sum(cumsum_d/sum_d <= threshold) + 1
#     }
#   else{
#     r = sum(SVD_X$d > threshold)
#   }
#   if(r>1){
#     U = SVD_X$u[,1:r]
#     D = diag(SVD_X$d[1:r])
#     D_inv = diag(1/SVD_X$d[1:r])
#     V = SVD_X$v[,1:r]
#   }else{
#     U = as.vector(SVD_X$u[,1:r])
#     D = as.matrix((SVD_X$d[1:r]))
#     D_inv = as.matrix(1/SVD_X$d[1:r])
#     V = as.vector(SVD_X$v[,1:r])
#     
#   }
#   
#   A_tilde = t(U)%*%Y%*%V%*%D_inv
#   eigen_A_tilde = eigen(A_tilde)
#   non_zero_eigen_val_count <- sum(abs(Re(eigen_A_tilde$values)) > threshold_eigen_A_tilde)
#   
#   phi = Y %*% V %*% D_inv %*% eigen_A_tilde$vectors[,1:non_zero_eigen_val_count] %*% diag(1/eigen_A_tilde$values[1:non_zero_eigen_val_count])
#   phi_inv = ginv(phi)
#   if(r>1){
#   A_hat_DMD = phi%*% diag(eigen_A_tilde$values[1:non_zero_eigen_val_count]) %*% phi_inv
#   }else{
#     A_hat_DMD = phi%*% (eigen_A_tilde$values[1:non_zero_eigen_val_count]) %*% phi_inv
#     
#   }
#   if(is.complex(A_hat_DMD)){
#     A_hat_DMD = Re(A_hat_DMD)
#   }
#   res <- list(A_hat = A_hat_DMD,
#               eigen_value = eigen_A_tilde$values,
#               eigen_vector = phi,
#               rank = r,
#               SVD_X_d = SVD_X$d)
#   return(res)
# }

DMD_alg <- function(output, method="ratio_sum", threshold=0.999, 
                         threshold_eigen_A_tilde=10^{-8}, s=1, d=1){
  # this function can implement both (exact) DMD and higher order DMD. 
  # the default setting is DMD
  # to implement higher order DMD, simply modify the value of s and d.
  # d is the number of timesteps we stack in HODMD;  
  # s is the timestep we skip.
  n = ncol(output)
  k = nrow(output)
  output_tilde <- matrix(NA, nrow=d*k, ncol=n-d+1)
  for(i in 1:(n-d+1)){
      output_tilde[,i] = as.vector(output[,i:(i+d-1)])
  }
  Y_index = seq(2, n-d+1, by=s)
  X = output_tilde[,(Y_index-1)]
  Y = output_tilde[, Y_index]

  SVD_X = svd(X)
  if(method=="ratio_sum"){
    sum_d = sum(SVD_X$d)
    cumsum_d = cumsum(SVD_X$d)
    r = sum(cumsum_d/sum_d <= threshold) + 1
  }
  else{
    r = sum(SVD_X$d > threshold)
  }
  if(r>1){
    U = SVD_X$u[,1:r]
    D = diag(SVD_X$d[1:r])
    D_inv = diag(1/SVD_X$d[1:r])
    V = SVD_X$v[,1:r]
  }else{
    U = as.vector(SVD_X$u[,1:r])
    D = as.matrix((SVD_X$d[1:r]))
    D_inv = as.matrix(1/SVD_X$d[1:r])
    V = as.vector(SVD_X$v[,1:r])
    
  }
  
  A_tilde = t(U)%*%Y%*%V%*%D_inv
  eigen_A_tilde = eigen(A_tilde)
  non_zero_eigen_val_count <- sum(abs(Re(eigen_A_tilde$values)) > threshold_eigen_A_tilde)
  
  phi = Y %*% V %*% D_inv %*% eigen_A_tilde$vectors[,1:non_zero_eigen_val_count] %*% diag(1/eigen_A_tilde$values[1:non_zero_eigen_val_count])
  phi_inv = ginv(phi)
  if(r>1){
    A_hat_DMD = phi%*% diag(eigen_A_tilde$values[1:non_zero_eigen_val_count]) %*% phi_inv
  }else{
    A_hat_DMD = phi%*% (eigen_A_tilde$values[1:non_zero_eigen_val_count]) %*% phi_inv
    
  }
  if(is.complex(A_hat_DMD)){
    A_hat_DMD = Re(A_hat_DMD)
  }
  
  in_sample_pred <- matrix(NA, ncol=(n-d), nrow=k)
  for(i in 1:(n-d)){
    in_sample_pred[,i] <- (A_hat_DMD %*% output_tilde[,i])[((d-1)*k+1):(d*k)]
  }
  rmse_in_sample <- mean((in_sample_pred - output[,(d+1):n])^2)
  res <- list(A_hat = A_hat_DMD,
              rank = r,
              in_sample_pred = in_sample_pred,
              rmse_in_sample = rmse_in_sample,
              eigen_value = eigen_A_tilde$values,
              eigen_vector = phi,
              rank = r,
              SVD_X_d = SVD_X$d,
              Y=Y,
              X=X)
  return(res)
  
  
}

# HODMD_alg <- function(output, d, s,
#                       threshold=0.999, 
#                       threshold_eigen_A_tilde=10^{-8}){
#   # this function is used to implement higher order DMD
#   n = ncol(output)
#   k = nrow(output)
#   output_tilde <- matrix(NA, nrow=d*k, ncol=n-d+1)
#   for(i in 1:(n-d+1)){
#     output_tilde[,i] = as.vector(output[,i:(i+d-1)])
#   }
#   DMD_fit <- DMD_alg(output_tilde, s=s, threshold=threshold,threshold_eigen_A_tilde=threshold_eigen_A_tilde)
#   A_hat <- DMD_fit$A_hat
#   
#   in_sample_pred <- matrix(NA, ncol=(n-d), nrow=k)
#   for(i in 1:(n-d)){
#     in_sample_pred[,i] <- (A_hat %*% output_tilde[,i])[((d-1)*k+1):(d*k)]
#   }
#   rmse_in_sample <- mean((in_sample_pred - output[,(d+1):n])^2)
#   res <- list(A_hat = A_hat,
#               r = DMD_fit$rank,
#               in_sample_pred = in_sample_pred,
#               rmse_in_sample = rmse_in_sample)
#   return(res)
# }
# 
# 
