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
# source('getData_reduced3_change.R')
source("functions_get_dic_cpp2.R") # c.get_loglik() loaded

# ===============================================
# Argument parser for command line input
# ===============================================
# args = commandArgs(trailingOnly = TRUE)
# filename = args[1]
# sigma_fixed = as.numeric(args[2])
# sigma_fixed_val = as.numeric(args[3])
# Ns = as.numeric(args[4])
# save_image = args[5]
# save_image_yes = as.numeric(args[6])
# SEED = as.numeric(args[7]) 
# tuneA_val = as.numeric(args[8])
# tuneD_val = as.numeric(args[9])
# tuneZs_val = as.numeric(args[10])
# cp_index = as.numeric(args[11])
# combn_cp2_index = as.numeric(args[12])

filename = '../../Data/nets_twitter.RData'
sigma_fixed = 1
sigma_fixed_val = 1
Ns = 20000
# Ns = 50000
save_image = 'None'
save_image_yes = 0
SEED = 13281
tuneA_val = 2
tuneD_val = 2
tuneZs_val = 4

###########
# Load data
############
load(filename)
cat(filename, '\n')
pi = ifelse(pi=='R', 2, 1)
n = dim(Y)[1]
TT = dim(Y)[3]

# set change points
change_points_chosen = 6
# change_points_chosen = c(5, 10)
# cp_index = 99; combn_cp2_index=6
# if(cp_index==100){
#   change_points_chosen = combn(3:10, 2)[,combn_cp2_index]
# }else{
#   change_points_chosen = c(cp_index:TT)
# }
# if(cp_index==99){
#   change_points_chosen = combn_cp2_index
# }

# cat densities over time
all = apply(Y, 3, FUN = function(Yt){
  mean(Yt[lower.tri(Yt)]) 
})
within_group1 = apply(Y[pi==1, pi==1,], 3, FUN = function(Yt){
  mean(Yt[lower.tri(Yt)])
})
within_group2 = apply(Y[pi==2, pi==2,], 3, FUN = function(Yt){
  mean(Yt[lower.tri(Yt)])
})
between_groups = apply(Y[pi==1, pi==2,], 3, mean)

cat('Densities over time:\n')
cat('all: ', all, '\n')
cat('within_group1: ', within_group1, '\n')
cat('within_group2: ', within_group2, '\n')
cat('between_groups: ', between_groups, '\n')
cat('dimension: ', dim(Y), '\n')


set.seed(SEED)
cat('SEED,', SEED, '\n')

# -------------------
# Parameters for MCMC 
# -------------------
burnin = round(Ns*(3/10))
p = 2   

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
to_plot_yearly_evolution_data_1cp = matrix(0, length(Y_segment_set), 21)
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
  to_plot_yearly_evolution_data_1cp[ii,] = c(ii, param1, param2, param3, param4, param5)
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



