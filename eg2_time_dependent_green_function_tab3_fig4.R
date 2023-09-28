#install.packages("BiocManager")
#BiocManager::install("rhdf5")

library(MASS)
library(rhdf5)
library(plot3D)
library(ggplot2)
library(cowplot)
library(RobustGaSP)
library(gridExtra)

source('functions/DMD_alg_Apr20.R')

###################### load data ######################

h5f = H5Fopen("data/td_tf6600.00.3_COHSEX.h5")
output_all_r <- rep(NA, 6601)
output_all_i <- rep(NA, 6601)

for(i in 1:4){
  for(j in 1:4){
    output_all_r = rbind(output_all_r, h5f$data$time_domain$lessG_all$r[,i,j,1])
    output_all_i = rbind(output_all_i, h5f$data$time_domain$lessG_all$i[,i,j,1])
    
  }
}
output_all_r = output_all_r[2:17,]
output_all_i = output_all_i[2:17,]

############## size of training/testing data ##############
n_training=2500
starting=n_training-1000
ending=n_training
n_testing = dim(output_all_r)[2]-n_training-1
n_testing1000=1000
n_testing2000=2000
k = nrow(output_all_r)


###################### PP-GP ######################
set.seed(0)
index_training=sample(starting:n_training,size=800) 
num_training=length(index_training)
output_training_r=matrix(NA,num_training,16)
input_training_r=matrix(NA,num_training,16)

# format input and output
count=0
for(t in index_training){
  count=count+1
  input_training_r[count,]=as.numeric(output_all_r[,t-1])
  output_training_r[count,]=as.numeric(output_all_r[,t])
}

# build emulator
m_ppgasp_r=ppgasp(design=input_training_r,
                  response=output_training_r,
                  nugget.est=T,isotropic=T,num_initial_values=2)
m_ppgasp_r@beta_hat 
m_ppgasp_r@nugget 

# out-of-sample forecast 
m_predict_record_r=matrix(NA,n_testing,16)
input_cur_r=output_all_r[,ending]
for(t_testing in 1:n_testing){
  m_predict_record_r[t_testing,]=predict(m_ppgasp_r,testing_input=matrix(input_cur_r,1,16))$mean
  
  input_cur_r= m_predict_record_r[t_testing,]
}

# sample M chains 
M=100
m_predict_record_sample_r=array(NA,dim=c(n_testing,16,M))
rt_sample_all_mat_r=matrix(rt(16*M*n_testing,df=m_ppgasp_r@num_obs-1),M*n_testing,16)

index_sample=1
for(i_M in 1:M){
  print(i_M)
  input_cur_r=output_all_r[,ending]
  input_cur_i=output_all_i[,ending]
  for(t_testing in 1:n_testing){
    m_pred_r=predict(m_ppgasp_r,testing_input=matrix(input_cur_r,1,16))
    m_predict_record_sample_r[t_testing,,i_M]=m_pred_r$mean+m_pred_r$sd*matrix(rt_sample_all_mat_r[index_sample,],1,16) ##mean and sd, it is a t in ppgp 
    input_cur_r = m_predict_record_sample_r[t_testing,,i_M]
    index_sample=index_sample+1
  }
}

m_predict_record_median_r=matrix(NA,n_testing,16)
m_predict_record_lower95_r=matrix(NA,n_testing,16)
m_predict_record_upper95_r=matrix(NA,n_testing,16)


for(t_testing in 1:n_testing){
  for(j in 1:16){
    median_lower_upper_r=quantile(m_predict_record_sample_r[t_testing,j,],c(0.025,0.5,0.975))
    m_predict_record_lower95_r[t_testing,j]=median_lower_upper_r[1]
    m_predict_record_median_r[t_testing,j]=median_lower_upper_r[2]
    m_predict_record_upper95_r[t_testing,j]=median_lower_upper_r[3]
    
    
  }
}


###################### AR(1) ######################
pred.fgasp.exp = matrix(0, k, n_testing2000)
upper95.fgasp.exp = matrix(0, k, n_testing2000)
lower95.fgasp.exp = matrix(0, k, n_testing2000)
input_fgasp = as.numeric(1:n_training)
testing_input = as.numeric((n_training+1):(n_training+n_testing2000))
for(grid_index in 1:k){
  output_fgasp = output_all_r[grid_index, 1:n_training]
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
# fit DMD
dmd_fit_r <- DMD_alg(output_all_r[,1:n_training]) 

A_hat_dmd_r <- dmd_fit_r$A_hat

# in-sample forecast
in_sample_pred_r <- A_hat_dmd_r %*% output_all_r[,1:(n_training-1)]
Sigma_dmd_r <- (output_all_r[,1] %*% t(output_all_r[,1]) + 
                  (output_all_r[,2:n_training]-in_sample_pred_r) %*% t(output_all_r[,2:n_training]-in_sample_pred_r))/n_training


# out-of-sample forecast
out_sample_pred_dmd_r <- matrix(NA, nrow=k, ncol=n_testing)
cur_state = output_all_r[,n_training]

# out-of-sample variance, lower&upper bound
for(i in 1:n_testing){
  out_sample_pred_dmd_r[,i] = A_hat_dmd_r %*% cur_state
  cur_state = out_sample_pred_dmd_r[,i]
}


pred_variance_dmd_r <- matrix(NA,nrow=k, ncol=n_testing)
cur_cov_dmd_r <- Sigma_dmd_r

for(i in 1:n_testing){
  pred_variance_dmd_r[,i] = diag(cur_cov_dmd_r)
  cur_cov_dmd_r <- A_hat_dmd_r %*% cur_cov_dmd_r %*% t(A_hat_dmd_r) + Sigma_dmd_r
}

ci_upper_dmd <- out_sample_pred_dmd_r + 1.96*sqrt(pred_variance_dmd_r)
ci_lower_dmd <- out_sample_pred_dmd_r - 1.96*sqrt(pred_variance_dmd_r)


###################### HODMD ######################
d=6
s=3
# fit HODMD
hodmd_fit_r <- DMD_alg(output_all_r[,1:n_training], d=d, s=s, threshold = 0.9999)
plot(hodmd_fit_r$SVD_X_d, type="l")
A_hat_hodmd_r <- hodmd_fit_r$A_hat

# in-sample forecast
in_sample_pred_hodmd <- matrix(NA, nrow=k, ncol=(n_training-d))
cur_state = matrix(c(output_all_r[,1:d]), ncol=1)
for(i in 1:ncol(in_sample_pred_hodmd)){
  in_sample_pred_hodmd[,i] = (A_hat_hodmd_r %*% cur_state)[(k*(d-1)+1):(k*d)]
  cur_state = matrix(c(output_all_r[,(i+1):(i+d)]), ncol=1)
}


# out-of-sample forecast
out_sample_pred_hodmd <- matrix(NA, nrow=k, ncol=n_testing)
cur_state <- matrix((output_all_r[, (n_training-d+1):n_training]), ncol=1)
for(i in 1:n_testing){
  out_sample_pred_hodmd[,i] = (A_hat_hodmd_r %*% cur_state)[(k*(d-1)+1):(k*d)]
  cur_state <- matrix(c(cur_state[-(1:k)], out_sample_pred_hodmd[,i]), ncol=1)
}


X_hodmd = hodmd_fit_r$X
Y_hodmd = hodmd_fit_r$Y
Sigma_hodmd <- (X_hodmd[,1] %*% t(X_hodmd[,1]) + (Y_hodmd - A_hat_hodmd_r %*% X_hodmd)%*%t(Y_hodmd - A_hat_hodmd_r %*% X_hodmd))/(ncol(X_hodmd)+1)


# variance of out-of-sample, lower&upper bound
pred_variance_hodmd <- matrix(NA,nrow=k, ncol=n_testing)
cur_cov_hodmd <- Sigma_hodmd

for(i in 1:n_testing){
  pred_variance_hodmd[,i] = diag(cur_cov_hodmd)[(k*(d-1)+1):(k*d)]
  cur_cov_hodmd <- A_hat_hodmd_r %*% cur_cov_hodmd %*% t(A_hat_hodmd_r) + Sigma_hodmd
}

ci_upper_hodmd <- out_sample_pred_hodmd + 1.96*sqrt(pred_variance_hodmd)
ci_lower_hodmd <- out_sample_pred_hodmd - 1.96*sqrt(pred_variance_hodmd)


##################### Table 3 #####################
### 1000-step forecast
RMSE_1000step = c(sqrt(mean( (pred.fgasp.exp[,1:n_testing1000]-output_all_r[,(n_training+1):(n_training+n_testing1000)])^2)),
                  sqrt(mean( ((out_sample_pred_dmd_r[,1:n_testing1000])-output_all_r[,(n_training+1):(n_training+n_testing1000)])^2)),
                  sqrt(mean( ((out_sample_pred_hodmd[,1:n_testing1000])-output_all_r[,(n_training+1):(n_training+n_testing1000)])^2)),
                  sqrt(mean( (t(m_predict_record_r[1:n_testing1000,])-output_all_r[,(n_training+1):(n_training+n_testing1000)])^2)))

P_1000step = c(length(which(upper95.fgasp.exp[,1:n_testing1000]>output_all_r[,(ending+1):(ending+n_testing1000)]&lower95.fgasp.exp[,1:n_testing1000]<output_all_r[,(ending+1):(ending+n_testing1000)]))/length(output_all_r[,(ending+1):(ending+n_testing1000)]),
               mean((output_all_r[,(n_training+1):(n_training+n_testing1000)] > ci_lower_dmd[,1:n_testing1000]) & (output_all_r[,(n_training+1):(n_training+n_testing1000)] < ci_upper_dmd[,1:n_testing1000])),
               mean((output_all_r[,(n_training+1):(n_training+n_testing1000)] > ci_lower_hodmd[,1:n_testing1000]) & (output_all_r[,(n_training+1):(n_training+n_testing1000)] < ci_upper_hodmd[,1:n_testing1000])),
               length(which(m_predict_record_upper95_r[1:n_testing1000,]>t(output_all_r[,(ending+1):(ending+n_testing1000)])&m_predict_record_lower95_r[1:n_testing1000,]<t(output_all_r[,(ending+1):(ending+n_testing1000)])))/length(output_all_r[,(ending+1):(ending+n_testing1000)]))

L_1000step = c(mean(abs(upper95.fgasp.exp[,1:n_testing1000]-lower95.fgasp.exp[,1:n_testing1000])),
               mean(ci_upper_dmd[,1:n_testing1000]-ci_lower_dmd[,1:n_testing1000]),
               mean(ci_upper_hodmd[,1:n_testing1000] - ci_lower_hodmd[,1:n_testing1000]) ,
               mean(abs(m_predict_record_upper95_r[1:n_testing1000,]-m_predict_record_lower95_r[1:n_testing1000,])))

measurement_1000step = data.frame(model = c("AR1", "DMD", "HODMD", "PP-GP"),
                                  RMSE = RMSE_1000step,
                                  P = P_1000step,
                                  L = L_1000step)
### 2000-step forecast
RMSE_2000step = c(sqrt(mean( (pred.fgasp.exp[,1:n_testing2000]-output_all_r[,(n_training+1):(n_training+n_testing2000)])^2)),
                  sqrt(mean( ((out_sample_pred_dmd_r[,1:n_testing2000])-output_all_r[,(n_training+1):(n_training+n_testing2000)])^2)),
                  sqrt(mean( ((out_sample_pred_hodmd[,1:n_testing2000])-output_all_r[,(n_training+1):(n_training+n_testing2000)])^2)),
                  sqrt(mean( (t(m_predict_record_r[1:n_testing2000,])-output_all_r[,(n_training+1):(n_training+n_testing2000)])^2)))

P_2000step = c(length(which(upper95.fgasp.exp[,1:n_testing2000]>output_all_r[,(ending+1):(ending+n_testing2000)]&lower95.fgasp.exp[,1:n_testing2000]<output_all_r[,(ending+1):(ending+n_testing2000)]))/length(output_all_r[,(ending+1):(ending+n_testing2000)]),
               mean((output_all_r[,(n_training+1):(n_training+n_testing2000)] > ci_lower_dmd[,1:n_testing2000]) & (output_all_r[,(n_training+1):(n_training+n_testing2000)] < ci_upper_dmd[,1:n_testing2000])),
               mean((output_all_r[,(n_training+1):(n_training+n_testing2000)] > ci_lower_hodmd[,1:n_testing2000]) & (output_all_r[,(n_training+1):(n_training+n_testing2000)] < ci_upper_hodmd[,1:n_testing2000])) ,
               length(which(m_predict_record_upper95_r[1:n_testing2000,]>t(output_all_r[,(ending+1):(ending+n_testing2000)])&m_predict_record_lower95_r[1:n_testing2000,]<t(output_all_r[,(ending+1):(ending+n_testing2000)])))/length(output_all_r[,(ending+1):(ending+n_testing2000)]))

L_2000step = c(mean(abs(upper95.fgasp.exp[,1:n_testing2000]-lower95.fgasp.exp[,1:n_testing2000])),
               mean(ci_upper_dmd[,1:n_testing2000]-ci_lower_dmd[,1:n_testing2000]),
               mean(ci_upper_hodmd[,1:n_testing2000] - ci_lower_hodmd[,1:n_testing2000]),
               mean(abs(m_predict_record_upper95_r[1:n_testing2000,]-m_predict_record_lower95_r[1:n_testing2000,])))

measurement_2000step = data.frame(model = c("AR1", "DMD", "HODMD", "PP-GP"),
                                  RMSE = RMSE_2000step,
                                  P = P_2000step,
                                  L = L_2000step)

##################### Fig.4 #####################
scientific_10 <- function(x) {
  x = signif(x,1)
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

plot_list = list()
plot_index_set = c(4,8,12,16)
num_plot=1000
num_skip=2000
for(i in 1:length(plot_index_set)){
  plot_index = plot_index_set[i]
  truth = output_all_r[plot_index,(num_skip+1):(n_training+num_plot)]
  ar1_pred = c(rep(NA, (n_training-num_skip)), pred.fgasp.exp[plot_index, 1:num_plot])
  dmd_pred = c(rep(NA, (n_training-num_skip)), out_sample_pred_dmd_r[plot_index,1:num_plot])
  hodmd_pred = c(rep(NA, (n_training-num_skip)), out_sample_pred_hodmd[plot_index,1:num_plot])
  ppgp_pred = c(rep(NA, (n_training-num_skip)), m_predict_record_r[1:num_plot, plot_index])
  
  plot1_df <- data.frame(step=rep((num_skip+1):(n_training+num_plot), 5),
                         method = factor(rep(c("Truth","AR(1)","DMD", "HODMD", "PP-GP"), each=(n_training+num_plot-num_skip)),levels=c("Truth","AR(1)","DMD", "HODMD", "PP-GP")),
                         pred = c(truth, ar1_pred, dmd_pred, hodmd_pred, ppgp_pred))
  ppgp_ci <- data.frame(step=(n_training+1):(n_training+num_plot),
                        lower = m_predict_record_lower95_r[1:num_plot,plot_index],
                        upper = m_predict_record_upper95_r[1:num_plot,plot_index])
  
  plot_list[[as.character(plot_index)]] =   ggplot() +
    geom_ribbon(data=ppgp_ci, aes(x=step, ymin=lower, ymax=upper), 
                fill="#00009c", alpha = 0.2) + 
    geom_path(data=plot1_df, aes(step, pred, group=method,
                                 col=method, linetype=method), linewidth=0.5) +
    scale_y_continuous(breaks=seq(min(na.omit(plot1_df$pred)), max(na.omit(plot1_df$pred)), length.out=4),
                       labels = scientific_10) + 
    scale_linetype_manual(values = c("solid","dashed","dashed", "dashed", "dashed")) +
    scale_color_manual(values=c("black","#916DB3", "#800020", "#FACC5F", "#4169E1")) + 
    geom_vline(xintercept=n_training, linetype="dotted") + 
    xlab('step') +
    ylab("state") +
    theme_classic() + 
    theme(legend.position="bottom", 
          legend.text = element_text(size = 13),
          legend.title= element_blank(),
          axis.title = element_text(size = 10), 
          axis.text = element_text(size=10)
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




