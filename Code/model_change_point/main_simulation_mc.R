######################################################################################
###   Dynamic Latent Space Network Models with Attractor and Edge Persistence    #####
######################################################################################
require("igraph")
require("MCMCpack")
require("inline")
require("RcppArmadillo")
require("vegan")
require("ROCR")

source("functions_adaptive2.R")
source('getData_reduced3_change.R')
source("functions_get_dic_cpp2.R") # c.get_loglik() loaded

# ===============================================
# Argument parser for command line input
# ===============================================
# args = commandArgs(trailingOnly = TRUE)
# SEED = as.numeric(args[1]) 
# n = as.numeric(args[2])
# TT = as.numeric(args[3])
# p = as.numeric(args[4])
# tau = as.numeric(args[5]); sigma = as.numeric(args[6])
# alpha = as.numeric(args[7]); delta = as.numeric(args[8])
# gamma1_1 = as.numeric(args[9]); # within-group attraction parameters
# gamma1_2 = as.numeric(args[10]); # within-group attraction parameters
# gamma2 = as.numeric(args[11]) # between-group attraction parameters
# Ns = as.numeric(args[12]) 
# tuneZs_val = as.numeric(args[13]) 
# sigma_fixed = as.numeric(args[14])
# sigma_fixed_val = as.numeric(args[15])
# gamma1_1_change = as.numeric(args[16]); # within-group attraction parameters, changed value in the middle
# gamma1_2_change = as.numeric(args[17]); # within-group attraction parameters, changed value in the middle
# gamma2_change = as.numeric(args[18]); 
# kk = as.numeric(args[19])
# change_points_chosen = as.numeric(args[20])

SEED = 123
n = 100; TT = 10; p = 2
Ns = 20000
tuneA_val = 2 ; tuneD_val = 2; tuneZs_val = 3
sigma_fixed = 1; sigma_fixed_val = 1
tau = 1; sigma = 1
alpha = 1; delta = 3
gamma1_1 = 0.6;        # within-group attraction parameters
gamma1_2 = 0.6;        # within-group attraction parameters
gamma2 = -0.2          # between-group attraction parameters
gamma1_1_change = 0.8; # within-group attraction parameters, changed value after the change-point
gamma1_2_change = 0.2; # within-group attraction parameters, changed value after the change-point
gamma2_change = -0.5; 
kk = 6                             # true change-point
# change_points_chosen = c(3, 7)     # specified change-point
change_points_chosen = 5           # specified change-point

############################################################
### MCMC settings
############################################################
tuneA_val = 2 
tuneD_val = 2
burnin = round(Ns/5)
save_image_yes = 0
save_image = 'None'

############################################################
### Get Synthetic Data
############################################################
# -------------------
# Parameters needed for generating data set
# -------------------
pi = rep(c(1, 2), each=n/2) # memberships, known
# one change point at kk
gamma1 = matrix(0, TT, 2)
gamma1[,1] = c(rep(gamma1_1, kk-1), rep(gamma1_1_change, TT+1-kk))
gamma1[,2] = c(rep(gamma1_2, kk-1), rep(gamma1_2_change, TT+1-kk))
gamma2_vec = c(rep(gamma2, kk-1), rep(gamma2_change, TT+1-kk))

# # e.g. multiple change points for TT=10
# gamma1 = matrix(0, 10, 2)
# gamma1[,1] = c(0.7, 0.7, 0.7, 0.7, 0.7, 0.5, 0.5, 0.5, 0.5, 0.5)
# gamma1[,2] = c(0.6, 0.6, 0.6, 0.6, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2)
# gamma2_vec = c(0.2, 0.2, 0.2, 0.2, -0.1, -0.1, -0.1, -0.5, -0.5, -0.5)

cat('gamma1_1:', gamma1[,1], '\n')
cat('gamma1_2:', gamma1[,2], '\n')
cat('gamma2_vec:', gamma2_vec, '\n')

# -------------------
# Generate dataset
# -------------------
# Y: Observed network sequence, dimension NxNxTT
cat('SEED:', SEED, '\n')
set.seed(SEED)
temp = getData_change_3(n, TT, p, tau, sigma, gamma1, gamma2_vec, alpha, delta, pi)
Y = temp$Y
Z_true = temp$X # Ground truth of latent positions, dimension nxpxTT, e.g. Z_true[,,t] for latent positions at time t

# -------------------
# Divide time-series into segments by change points
# -------------------
Y_segment_set = list()
if(length(change_points_chosen)==0){
  cat('No change point is given. Fit on whole Y...')
  Y_segment_set[[1]] = Y
}
if(length(change_points_chosen)==1){
  cat('change_points_chosen:', change_points_chosen, '\n')
  cat('One change point is given. Fit on\n')
  Y_segment_set[[1]] = Y[,,1:(change_points_chosen-1)]
  Y_segment_set[[2]] = Y[,,(change_points_chosen-1):TT] # for the second segment, sample conditional on time change_points_chosen-1
  cat(1:(change_points_chosen-1),'\n')
  cat((change_points_chosen-1):TT, '\n')
}
if(length(change_points_chosen)>1){
  cat('change_points_chosen:', change_points_chosen, '\n')
  cat('Multiple change points are given. Fit on\n')
  len_cp = length(change_points_chosen)
  Y_segment_set[[1]] = Y[,,1:(change_points_chosen[1]-1)]
  cat(1:(change_points_chosen[1]-1), '\n')
  for(i in 2:length(change_points_chosen)){
    Y_segment_set[[i]] = Y[,,(change_points_chosen[i-1]-1):(change_points_chosen[i]-1)]
    cat((change_points_chosen[i-1]-1):(change_points_chosen[i]-1),'\n')
  }
  Y_segment_set[[len_cp+1]] = Y[,,(change_points_chosen[len_cp]-1):TT]
  cat((change_points_chosen[len_cp]-1):TT, '\n')
}

############################################################
### Fit models for each Y in Y_segment_set
############################################################
DICs = numeric(length(Y_segment_set))
for(ii in 1:length(Y_segment_set)){
  Y = Y_segment_set[[ii]]
  subwindow = ii
  cat('--------------------------------\n')
  cat('subwindow: ', ii, '\n')
  cat('Dim of Y in subwindow:', dim(Y), '\n')
  n = dim(Y)[1]
  TT = dim(Y)[3]
  ### Get MCMC Samples and posterior analyze
  source('runMCMC.R')
  DICs[ii] = DIC_temp
  cat('\n\n')
  # save image
  if(save_image_yes){
    out_name = paste('homedirec/sw_', subwindow, '_', save_image, sep='')
    save.image(out_name)
    cat('Image was saved in', out_name, '\n')
  }
  cat('\n\n')
}
cat(DICs, '\n')
cat('Sum of DICs is: ', sum(DICs), '\n')



