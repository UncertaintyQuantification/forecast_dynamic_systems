### code for Example1: Lorenz 96
library(deSolve) 
library(pracma)
library(plot3D)
library(RobustGaSP)
library(lhs)
library(lattice)
library(MASS)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(gridExtra)


source('functions/DMD_alg_Apr20.R')

################# simulation of Lorenz 96 #################
Lorenz96_model<-function(time, state, parameters){
  # parameters = c(F, the sd of a white noise)
  x <- state
  dx <- rep(0, k)
  dx[1] = (x[2] - x[k-1])*x[k] - x[1] + parameters
  dx[2] = (x[3] - x[k])*x[1] - x[2] + parameters
  dx[k] = (x[1] - x[k-2])*x[k-1] - x[k] + parameters
  for(i in 3:(k-1))
    dx[i] = (x[i+1] - x[i-2])*x[i-1] - x[i] + parameters
  list(dx)
}

pred_rk4_ppgp_local<-function(model,state_cur,testing_input_local,h,mean_pred=T){
  ### predictive function for PP-GP
  if(mean_pred){
    k1=predict(model, testing_input_local)$mean 
    testing_input_local=matrix((state_cur+k1*h/2)[loc_neighbor_index],k,4)
    k2=predict(model, testing_input_local)$mean
    testing_input_local=matrix((state_cur+k2*h/2)[loc_neighbor_index],k,4)
    k3=predict(model, testing_input_local)$mean 
    testing_input_local=matrix((state_cur+k3*h)[loc_neighbor_index],k,4)
    k4=predict(model, testing_input_local)$mean 
  }else{
    ##simulate
    k1=simulate(model,testing_input_local,num_sample=1)
    testing_input_local=matrix((state_cur+k1*h/2)[loc_neighbor_index],k,4)
    k2=simulate(model,testing_input_local,num_sample=1)
    testing_input_local=matrix((state_cur+k2*h/2)[loc_neighbor_index],k,4)
    k3=simulate(model,testing_input_local,num_sample=1)
    testing_input_local=matrix((state_cur+k3*h)[loc_neighbor_index],k,4)
    k4=simulate(model,testing_input_local,num_sample=1)
  }
  return(state_cur+(k1+2*k2+2*k3+k4)*h/6 )
}

k = 40 # dimension of snapshot
N = 1*10^3 # total time steps (training + testing)

set.seed(123)
p0 = (rWishart(1, df=k, Sigma=diag(k)))[,,1] # the initial sates x_0 is drawn from N(0, p0), where p0 is drawn from a Wishart distribution
L_p0 = t(chol(p0))
x0 <- as.vector(L_p0%*%matrix(rnorm(k),k,1))
t = 0 
h = 0.01 # step-size for euler/rk4,

times <- seq(from=0,by=h, length.out=N+1)
t <- times[-1]
out_ode <- matrix(0, k, N+1)
out_ode[,1] <- x0
tau = 1

for(j in 1:N){
  ##for rk4 it uses smoothing obs to solve
  x_new <- ode(y = out_ode[,j], times = c((j-1)*h,  j*h), func = Lorenz96_model, parms = 8, method = "rk4")
  
  x_new <- t(as.matrix(x_new[,-1]))
  out_ode[,j+1] <- x_new[,2]
}

X_true <- out_ode[,-1] 

###################### PP-GP ######################
## build local index
location_index=as.numeric(1:k)
loc_neighbor_index=matrix(NA,k,4) ###the true one has 4 neighbor
for(i in 1:k){
  if(i==1){
    loc_neighbor_index[i,]=c(k-1,k,1,2)
  }else if(i==2)(
    loc_neighbor_index[i,]=c(k,1,2,3)
    
  )else if(i==k){
    loc_neighbor_index[i,]=c(k-2,k-1,k,1)
  }else{
    loc_neighbor_index[i,]=c(i-2,i-1,i,i+1)
    
  }
}

##form input 
n=500
input=matrix(NA,n,4)
output=rep(NA,n)
n_training = 100
n_testing = N - n_training
f_all=matrix(NA,k,N)
ground_truth=t(X_true[,((n_training+1):N)])


for(i in 1:N){
  f_all[,i]=Lorenz96_model(i, state= X_true[,i], parameters=8)[[1]]
}
f_training=f_all[,1:n_training] # k times n_training
set.seed(0)
index_training=sample(1:(n_training*k),n) # a vector of length n(=500)
output=f_training[index_training] 
input=matrix(NA,n,4)
count=1
for(i in index_training){
  loc_i=i%%40
  t_i=floor((i-1)/40)+1
  if(loc_i==0){loc_i=40}
  input[count,]=X_true[loc_neighbor_index[loc_i,],t_i]
  count=count+1
}

splom((input),(input),pch=20,cex=.5)

## build emulator
m_f=rgasp(design=input,response=(output))

## out-of-sample forecast, 900 steps
state_cur=X_true[,n_training] 
testing_input_cur=matrix(state_cur[loc_neighbor_index],k,4)

pred_ppgp_local_record=matrix(NA,N-n_training,k)
for(i in 1: (N-n_training)){
  print(i)
  pred_ppgp_local_record[i,]=pred_rk4_ppgp_local(model=m_f,state_cur=state_cur,testing_input_local=testing_input_cur,h=h,mean_pred=T)
  state_cur=pred_ppgp_local_record[i,]
  testing_input_cur=matrix(state_cur[loc_neighbor_index],k,4)
  
}


## sample M chains
M=100
pred_ppgp_sample_local_record=array(NA,dim=c(N-n_training,k,M))
for(i_M in 1:M){
  print(i_M)
  state_cur=X_true[,n_training]
  testing_input_cur=matrix(state_cur[loc_neighbor_index],k,4)
  
  for(i in 1: (N-n_training)){
    pred_ppgp_sample_local_record[i,,i_M]=pred_rk4_ppgp_local(model=m_f,state_cur=state_cur,testing_input_local=testing_input_cur,h=h,mean_pred=F)
    state_cur=pred_ppgp_sample_local_record[i,,i_M]
    testing_input_cur=matrix(state_cur[loc_neighbor_index],k,4)
    
  }
}

mean_sample=matrix(NA,N-n_training,k)
median_sample=matrix(NA,N-n_training,k)
lower95_sample=matrix(NA,N-n_training,k)
upper95_sample=matrix(NA,N-n_training,k)
sd_sample = matrix(NA, N-n_training, k)

## draw samples
for(i in 1: (N-n_training)){
  for(j in 1:k){
    mean_sample[i,j]=mean(pred_ppgp_sample_local_record[i,j,])
    lower_median_upper=quantile(pred_ppgp_sample_local_record[i,j,],c(0.025,0.5,0.975))
    sd_sample[i,j] = sd(pred_ppgp_sample_local_record[i,j,])
    lower95_sample[i,j]=lower_median_upper[1]
    median_sample[i,j]=lower_median_upper[2]
    upper95_sample[i,j]=lower_median_upper[3]
    
  }
}


###################### AR(1) ######################
pred.fgasp.exp = matrix(0, k, n_testing)
lower95.fgasp.exp=matrix(0, k, n_testing)
upper95.fgasp.exp=matrix(0, k, n_testing)
input_fgasp = as.numeric(1:n_training) 
testing_input=as.numeric((n_training+1):N)
for(grid_index in 1:k){
  output_fgasp = X_true[grid_index, 1:n_training]
  fgasp.exp = fgasp(input_fgasp, output_fgasp, kernel_type = "exp",
                    have_noise = T)
  est_all = optim(c(log(1),log(.02)),log_lik,object=fgasp.exp,
                  method="L-BFGS-B",
                  control = list(fnscale=-1))
  fgasp.pred=predict(param=est_all$par, object=fgasp.exp,testing_input=testing_input)
  pred.fgasp.exp[grid_index,]=fgasp.pred@mean
  lower95.fgasp.exp[grid_index,]=fgasp.pred@mean-1.96*sqrt(fgasp.pred@var)
  upper95.fgasp.exp[grid_index,]=fgasp.pred@mean+1.96*sqrt(fgasp.pred@var)
}



###################### DMD ######################
## fit DMD
DMD_fit <- DMD_alg(X_true[,1:n_training], threshold = 0.99) 
A_hat_dmd <- DMD_fit$A_hat

## in-sample forecast
in_sample_pred <- A_hat_dmd %*% X_true[,1:(n_training-1)]
Sigma_dmd <- (X_true[,1] %*% t(X_true[,1]) + (X_true[,2:n_training]-in_sample_pred) %*% t(X_true[,2:n_training]-in_sample_pred))/n_training

## out-of-sample forecast
out_sample_pred_dmd <- matrix(NA, k, n_testing)
cur_state = X_true[,n_training]
for(i in 1:n_testing){
  out_sample_pred_dmd[, i] <- A_hat_dmd %*% cur_state
  cur_state <- out_sample_pred_dmd[,i]
}

## variance of out-of-sample, lower&upper bound
pred_variance_dmd <- matrix(NA,nrow=k, ncol=n_testing)
cur_cov_dmd_r <- Sigma_dmd
for(i in 1:n_testing){
  pred_variance_dmd[,i] = diag(cur_cov_dmd_r)
  cur_cov_dmd_r <- A_hat_dmd %*% cur_cov_dmd_r %*% t(A_hat_dmd) + Sigma_dmd
}

ci_upper_dmd <- out_sample_pred_dmd + 1.96*sqrt(pred_variance_dmd)
ci_lower_dmd <- out_sample_pred_dmd - 1.96*sqrt(pred_variance_dmd)


###################### HODMD ######################
d=6
s=3
hodmd_fit <- DMD_alg(X_true[,1:n_training], d=d, s=s, threshold = 0.99)
A_hat_hodmd <- hodmd_fit$A_hat

## in-sample prediction
in_sample_pred_hodmd <- matrix(NA, nrow=k, ncol=(n_training-d))
cur_state = matrix(c(X_true[,1:d]), ncol=1)
for(i in 1:ncol(in_sample_pred_hodmd)){
  in_sample_pred_hodmd[,i] = (A_hat_hodmd %*% cur_state)[(k*(d-1)+1):(k*d)]
  cur_state = matrix(c(X_true[,(i+1):(i+d)]), ncol=1)
}

X_hodmd = hodmd_fit$X
Y_hodmd = hodmd_fit$Y
Sigma_hodmd <- (X_hodmd[,1] %*% t(X_hodmd[,1]) + (Y_hodmd - A_hat_hodmd %*% X_hodmd)%*%t(Y_hodmd - A_hat_hodmd %*% X_hodmd))/(ncol(X_hodmd)+1)

## out-of-sample prediction
out_sample_pred_hodmd <- matrix(NA, nrow=k, ncol=n_testing)
cur_state <- matrix((X_true[, (n_training-d+1):n_training]), ncol=1)
for(i in 1:n_testing){
  out_sample_pred_hodmd[,i] = (A_hat_hodmd %*% cur_state)[(k*(d-1)+1):(k*d)]
  cur_state <- matrix(c(cur_state[-(1:k)], out_sample_pred_hodmd[,i]), ncol=1)
}

## variance of out-of-sample prediction, compute lower/upper bound
pred_variance_hodmd <- matrix(NA,nrow=k, ncol=n_testing)
cur_cov_hodmd = Sigma_hodmd
for(i in 1:n_testing){
  pred_variance_hodmd[,i] = diag(cur_cov_hodmd)[(k*(d-1)+1):(k*d)]
  cur_cov_hodmd <- A_hat_hodmd %*% cur_cov_hodmd %*% t(A_hat_hodmd) + Sigma_hodmd
}

ci_upper_hodmd <- out_sample_pred_hodmd + 1.96*sqrt(pred_variance_hodmd)
ci_lower_hodmd <- out_sample_pred_hodmd - 1.96*sqrt(pred_variance_hodmd)

############## Table 2 ##############
### 500-step forecast
RMSE_500step = c(sqrt(mean( ((pred.fgasp.exp[,1:500])-t(ground_truth[1:500,]))^2)),
                 sqrt(mean( ((out_sample_pred_dmd[,1:500])-t(ground_truth[1:500,]))^2)),
                 sqrt(mean( ((out_sample_pred_hodmd[,1:500])-t(ground_truth[1:500,]))^2)),
                 sqrt(mean( (pred_ppgp_local_record[1:(500),]-ground_truth[1:(500),])^2)))
P_500step = c(mean((t(ground_truth[1:500,]) > lower95.fgasp.exp[,1:500]) & (t(ground_truth[1:500,]) < upper95.fgasp.exp[,1:500])),
              mean((t(ground_truth[1:500,]) > ci_lower_dmd[,1:500]) & (t(ground_truth[1:500,]) < ci_upper_dmd[,1:500])),
              mean((t(ground_truth[1:500,]) > ci_lower_hodmd[,1:500]) & (t(ground_truth[1:500,]) < ci_upper_hodmd[,1:500])),
              length(which(ground_truth[1:(500),]<upper95_sample[1:(500),] &ground_truth[1:(500),]>lower95_sample[1:(500),]))/length(ground_truth[1:(500),]))
L_500step = c(mean((upper95.fgasp.exp-lower95.fgasp.exp)[,1:500]),
              mean((ci_upper_dmd-ci_lower_dmd)[,1:500]),
              mean((ci_upper_hodmd-ci_lower_hodmd)[,1:500]),
              mean(abs(upper95_sample-lower95_sample)[1:(500),]))
measurement_500step = data.frame(model = c("AR1", "DMD", "HODMD", "PP-GP"),
                                 RMSE = RMSE_500step,
                                 P = P_500step,
                                 L = L_500step)
### 900-step forecast
RMSE_900step = c(sqrt(mean( ((pred.fgasp.exp)-t(ground_truth))^2)),
                 sqrt(mean( ((out_sample_pred_dmd)-t(ground_truth))^2)),
                 sqrt(mean( ((out_sample_pred_hodmd)-t(ground_truth))^2)),
                 sqrt(mean( (pred_ppgp_local_record-ground_truth)^2)))
P_900step = c(mean((t(ground_truth) > lower95.fgasp.exp) & (t(ground_truth) < upper95.fgasp.exp)),
              mean((t(ground_truth) > ci_lower_dmd) & (t(ground_truth) < ci_upper_dmd)),
              mean((t(ground_truth) > ci_lower_hodmd) & (t(ground_truth) < ci_upper_hodmd)),
              length(which(ground_truth<upper95_sample &ground_truth>lower95_sample))/length(ground_truth))
L_900step = c(mean((upper95.fgasp.exp-lower95.fgasp.exp)),
              mean((ci_upper_dmd-ci_lower_dmd)),
              mean((ci_upper_hodmd-ci_lower_hodmd)),
              mean(abs(upper95_sample-lower95_sample)))
measurement_900step = data.frame(model = c("AR1", "DMD", "HODMD", "PP-GP"),
                                 RMSE = RMSE_900step,
                                 P = P_900step,
                                 L = L_900step)


############## Fig. 1  ##############
N0=n_training
plot_index_set = c(10,20,30,40)
plot_list = list()
for(plot_index in plot_index_set){
  ppgp_pred_df = data.frame(time = (N0+1):N,
                            point_est = mean_sample[, plot_index], # or pred_ppgp_local_record?
                            lower = lower95_sample[, plot_index],
                            upper = upper95_sample[, plot_index]) 
  ar1_pred_df = data.frame(time = (N0+1):N,
                           point_est = pred.fgasp.exp[plot_index,])
  dmd_pred_df = data.frame(time = (N0+1):N,
                           point_est = out_sample_pred_dmd[plot_index,])
  
  hodmd_pred_df = data.frame(time = (N0+1):N,
                             point_est = out_sample_pred_hodmd[plot_index,])
  
  point_df = data.frame(time = rep(1:N, 5),
                        method = factor(rep(c("Truth","AR(1)", "DMD", "HODMD", "PP-GP"), each=N),levels=c("Truth","AR(1)","DMD", "HODMD", "PP-GP")),
                        value = c(X_true[plot_index,],
                                  c(rep(NA, N0),pred.fgasp.exp[plot_index,]),
                                  c(rep(NA, N0),out_sample_pred_dmd[plot_index,]), 
                                  c(rep(NA, N0),out_sample_pred_hodmd[plot_index,]),
                                  c(rep(NA, N0),mean_sample[, plot_index])))
  
  plot_list[[as.character(plot_index)]] = ggplot() +
    geom_ribbon(data=ppgp_pred_df, aes(x=time, ymin=lower, ymax=upper), 
                fill="#00009c", alpha = 0.2) +
    geom_path(data=point_df, aes(time, value, group=method,
                                 col=method, linetype=method), linewidth=1.1) +
    scale_linetype_manual(values = c("solid","dashed","dashed", "dashed", "dashed")) +
    
    scale_color_manual(values=c("black","#E48586", "#B28264", "#FACC5F", "#00009c")) + 
    geom_vline(xintercept=N0, linetype="dotted") + 
    geom_vline(xintercept=lyapunov_exp[plot_index/10]*100,
               color="red", linetype="dotted") + 
    xlab('step') +
    ylab("state") +
    theme_classic() + 
    theme(legend.position="bottom", 
          legend.text = element_text(size = 15),
          legend.title= element_blank(),
          axis.title = element_text(size = 14), 
          axis.text = element_text(size=14)
    )
}
shared_legend <- get_legend(plot_list[[as.character(plot_index_set[1])]])
p1 = plot_grid(plot_list[[as.character(plot_index_set[1])]] + theme(legend.position = "none"),
               plot_list[[as.character(plot_index_set[2])]] + theme(legend.position = "none"),
               plot_list[[as.character(plot_index_set[3])]] + theme(legend.position = "none"),
               plot_list[[as.character(plot_index_set[4])]] + theme(legend.position = "none"),
               nrow=2, ncol=2, align = 'vh')
plot_grid(shared_legend, p1, ncol=1, 
          rel_heights = c(0.2, 1))


