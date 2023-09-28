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


################# Fig.3 #################
sample_ind=1
pred_ppgp_sample_local_record_sd = matrix(NA, 40, 900)
state_cur = X_true[,100]
testing_input_cur = matrix(state_cur[loc_neighbor_index],k,4)
for(i in 1:900){
  pred_ppgp_sample_local_record_sd[,i] = predict(m_f, testing_input_cur)$sd
  state_cur = pred_ppgp_sample_local_record[i,,sample_ind]
  testing_input_cur = matrix(state_cur[loc_neighbor_index],k,4)
}

se_1_df <- data.frame(step=101:1000,
                      se = pred_ppgp_sample_local_record_sd[1,],
                      abs_pred_error = abs(mean_sample[,1]-X_true[1,101:1000])) %>%
  mutate(cum_abs_pred_error = cumsum(abs_pred_error),
         mean_abs_pred_error = cum_abs_pred_error/seq(1,900,1))

p1_se= ggplot(se_1_df, aes(x=step)) +
  geom_line(aes(y=se)) +
  geom_line(aes(y=mean_abs_pred_error/450), col="#00009c") + 
  scale_y_continuous(
    name = "standard deviation",
    limits = c(0,0.0018),
    sec.axis = sec_axis( trans=~.*450, name="cummulative mean abs error")
  ) + 
  theme_classic() +
  theme(axis.title.y.right = element_text(colour = "#00009c"),
        plot.margin = unit(c(0.5,0.5,0.5,0),"cm")) 

se_21_df <- data.frame(step=101:1000,
                       se = pred_ppgp_sample_local_record_sd[21,],
                       abs_pred_error = abs(mean_sample[,21]-X_true[21,101:1000])) %>%
  mutate(cum_abs_pred_error = cumsum(abs_pred_error),
         log_cum_abs_pred_error = log(cum_abs_pred_error),
         mean_abs_pred_error = cum_abs_pred_error/seq(1,900,1))

p21_se = ggplot(se_21_df, aes(x=step)) +
  geom_line(aes(y=se)) +
  geom_line(aes(y=mean_abs_pred_error/200),col="#00009c") + 
  scale_y_continuous(
    name = "standard deviation",
    limits = c(0,0.004),
    sec.axis = sec_axis( trans=~.*200, name="cummulative mean abs error")
  ) + 
  theme_classic() +
  theme(axis.title.y.right = element_text(colour = "#00009c"),
        plot.margin = unit(c(0.5,0,0.5,0.5),"cm"))
(p21_se)

grid.arrange(p1_se, p21_se, nrow=1)
