n = dim(Y)[1]
TT = dim(Y)[3]
# -------------------
# Set tuning parameters for MCMC 
# -------------------
tuneZs = matrix(rep(tuneZs_val, n*TT), c(n, TT))
tuneA = tuneA_val
tuneD = tuneD_val

#-------------------------------------
## Priors and Initialization for MCMC  
#-------------------------------------
Alpha = numeric(Ns)   
Delta = numeric(Ns)
w = matrix(0, n, Ns) # first/secondrow for gamma1_1/gamma1_2, third row for gamma2, each column contains samples for one iteration
t2 = numeric(Ns)
s2 = numeric(Ns)
X = list()  
AccRate = numeric(4+n*TT)
# ---------1. initialize latent positions X (GMDS, Sarkar and Moore, 2005) --------------- # 
if(subwindow==1){
  dissim = array(0, dim=dim(Y))
  for(tt in 1:TT) dissim[,,tt] = shortest.paths(graph=graph.adjacency(Y[,,tt]), mode="all")
  dissim[which(dissim==Inf, arr.ind=TRUE)] = max(c(dissim[which(dissim!=Inf, arr.ind=TRUE)]))
  X[[1]] = array(0, c(p, TT, n)) # latent positions at iteration 1, each layer for each agent, each columns are positions at each time
  X[[1]][,1,] = t(cmdscale(d=dissim[,,1], k=p)) # centered
  temp.lambda = 1
  H = matrix(-1/n ,n, n)+diag(n)
  for(tt in 2:TT){
    temp = 1/(1+temp.lambda)*H%*%(-1/2*dissim[,,tt]^2)%*%H +
      temp.lambda/(1+temp.lambda)*t(X[[1]][,tt-1,])%*%X[[1]][,tt-1,]
    temp = eigen(temp)
    X[[1]][,tt,] = t(temp$vectors[,1:p]%*%diag(temp$values[1:p])^(1/2))
    X[[1]][,tt,] = t(vegan::procrustes(X=t(X[[1]][,tt-1,]), Y=t(X[[1]][,tt,]), scale=FALSE)$Yrot)
  }
}else{
  # initialize using the estimated values from last sub-window
  X[[1]] = array(0, c(p, TT, n)) # latent positions at iteration 1, each layer for each agent, each columns are positions at each time
  for (t in 1:TT){
    for (i in 1:n){
      X[[1]][,t,i] = EX[,tt,i]
    }
  }
}

# ---------2. initialize alpha, delta, gamma1_1, gamma1_2, gamma2, t2, s2 --------------- # 
Alpha[1] = 0
Delta[1] = 0
w[1,1] = 0.5
w[2,1] = 0.5
w[3,1] = -0.5
t2[1] = sum(X[[1]][,1,]*X[[1]][,1,])/(n*p) # t2 is redundant in the fitting period for subwindow>1
s2[1] = 0.0001  # redundant since s2 is fixed. 

# ---------3. set priors --------------- # 
# alpha, delta (flat)
nuAlpha = Alpha[1]
xiAlpha = 100
nuDelta = Delta[1]
xiDelta = 100
# gamma1, gamma2 (informative)
nuGamma1 = 0.5
xiGamma1 = 10
nuGamma2 = -0.5
xiGamma2 = 10
# t2: Tau^2 ~ IG()
shapeT2 = 2.05
scaleT2 = (shapeT2-1)*t2[1]
# s2: redundant since s2 is fixed
shapeS2 = 9
scaleS2 = 0.0001

#-------------------------------------
## Run MCMC
#-------------------------------------
cat('MCMC settings: \n')
cat('      n, TT, Ns, burnin, tuneA, tuneD, tuneZs[2,1] = ', c(n, TT, Ns, burnin, tuneA, tuneD, tuneZs[2,1]), '\n')
cat('sigma_fixed', sigma_fixed, '\n')
cat('sigma_fixed value', sigma_fixed_val, '\n')
cat('Init for alpha, delta, gamma1_1, gamma1_2, gamma2, t2, s2 = ', Alpha[1], Delta[1], w[1,1], w[2,1], w[3,1], t2[1], s2[1], '\n')
cat('Priors for alpha: nuAlpha, xiAlpha: ', nuAlpha, xiAlpha, '\n')
cat('Priors for Delta: nuDelta, xiDelta: ', nuDelta, xiDelta, '\n')
cat('Priors for gamma1s: nuGamma1, xiGamma1: ', nuGamma1, xiGamma1, '\n')
cat('Priors for gamma2: nuGamma2, xiGamma2: ', nuGamma2, xiGamma2, '\n')
cat('Priors for tau2: IG(shapeT2, scaleT2): ', shapeT2, scaleT2, '\n')
# cat('Priors for sigma2: IG(shapeS2, scaleS2): ', shapeS2, scaleS2, '\n')

pb = txtProgressBar(min=2, max=Ns, style=3)
system.time({
  for(it in 2:Ns){
    # ------------ 1. Draw Samples for X, alpha, delta-------------------------- #
    RN = rnorm(n*TT*p) # random numbers for the update of X_ti, t=1,...,TT, i=1,...,n
    RNAD = rnorm(2)    # random numbers for the update of alpha, delta
    Draws1 = c.update1(X[[it-1]], c(n,p,TT,1), tuneZs, Y,
                       Alpha[it-1], Delta[it-1], tuneA, tuneD, w[,it-1],
                       t2[it-1], s2[it-1], xiAlpha, xiDelta, nuAlpha,
                       nuDelta, CAUCHY=0, RN, RNAD, pi, it, subwindow)
    X[[it]] = Draws1[[1]] # pxTTxn, Xit, i.e. Zti 
    Alpha[it] = Draws1[[2]]
    Delta[it] = Draws1[[3]]
    AccRate = AccRate + Draws1[[4]]
    tuneA = Draws1[[5]]
    tuneD = Draws1[[6]]
    tuneZs = Draws1[[7]]
    
    if(it==burnin){ #using latent positions at burnin-th iteration(centered) as reference
      Xit0 <- t(X[[it]][,1,])
      for(tt in 2:TT) Xit0 <- rbind(Xit0, t(X[[it]][,tt,])) # Xito (n*TT X p)
    }
    if(it>burnin){
      XitCentered <- t(X[[it]][,1,]) # nXp
      for(tt in 2:TT) XitCentered <- rbind(XitCentered,t(X[[it]][,tt,]))
      procr <- vegan::procrustes(X=Xit0,Y=XitCentered,scale=FALSE)$Yrot
      for(tt in 1:TT){
        X[[it]][,tt,] <- t(procr[((tt-1)*n+1):(tt*n),]) # X[[it]] p X TT X n
      }
    }
    if(it < Ns) X[[it+1]] = X[[it]]
    
    # ------------ 2. Draw Samples for tau^2, sigma^2-------------------------- #
    Draws2 <- c.t2s2Parms(X[[it]], c(n,p,TT,1), shapeT2,
                          shapeS2, scaleT2, scaleS2, Y, w[,it-1], pi)
    if(subwindow==1){
      t2[it] <- rinvgamma(1, shape=Draws2[[1]], scale=Draws2[[2]])
      s2[it] <- rinvgamma(1, shape=Draws2[[3]], scale=Draws2[[4]])
    }else{
      t2[it] = t2[it-1]
      # Don not run for subwindow>1, sigma not fixed
    }
    if (sigma_fixed){
      s2[it] = sigma_fixed_val^2
    }
    Ai_1s = Draws2[[5]] # saved intermediate results, of dimension p x TT x n
    Ai_2s = Draws2[[6]]
    
    # ------------ 3. Draw Samples for gamma_w, gamma_b -------------------------- #
    # Draw samples for gamma1's
    ati = X[[it]][,-1,] - X[[it]][,-TT,] - w[3,it-1]*Ai_2s[,-1,] # px(TT-1)*n
    bti = Ai_1s[,-1,]
    # gamma1_1
    if(length(dim(ati))==3){
      sum_ab = sum(ati[,,pi==1]*bti[,,pi==1])
      sum_bb = sum(bti[,,pi==1]^2)
    }else{ # when TT=2
      sum_ab = sum(ati[,pi==1]*bti[,pi==1])
      sum_bb = sum(bti[,pi==1]^2)
    }
    mu_g1 = (sum_ab/s2[it] + nuGamma1/xiGamma1)/(sum_bb/s2[it] + 1/xiGamma1)
    s_g1 = sqrt(1/((sum_bb/s2[it]) + (1/xiGamma1)))
    w[1,it] = rnorm(1, mu_g1, s_g1)
    # gamma1_2
    if(length(dim(ati))==3){
      sum_ab = sum(ati[,,pi==2]*bti[,,pi==2])
      sum_bb = sum(bti[,,pi==2]^2)
    }else{
      sum_ab = sum(ati[,pi==2]*bti[,pi==2])
      sum_bb = sum(bti[,pi==2]^2)
    }
    mu_g1 = (sum_ab/s2[it] + nuGamma1/xiGamma1)/(sum_bb/s2[it] + 1/xiGamma1)
    s_g1 = sqrt(1/((sum_bb/s2[it]) + (1/xiGamma1)))
    w[2,it] = rnorm(1, mu_g1, s_g1)
    
    # Draw samples for gamma2
    ati = X[[it]][,-1,] - X[[it]][,-TT,] - rep(w[pi[1:n],it], each=p*(TT-1))*Ai_1s[,-1,] # px(TT-1)*n
    bti = Ai_2s[,-1,]
    sum_ab = sum(ati*bti)
    sum_bb = sum(bti^2)
    mu_g2 = (sum_ab/s2[it] + nuGamma2/xiGamma2)/(sum_bb/s2[it] + 1/xiGamma2)
    s_g2 = sqrt(1/((sum_bb/s2[it]) + (1/xiGamma2)))
    w[3,it] = rnorm(1, mu_g2, s_g2)
    
    setTxtProgressBar(pb,it)
  }
  close(pb)
})
# Acceptance rate
cat('Acceptance rate: alpha, delta: ', AccRate[1:2]/(it-1), '\n')
cat('Acceptance rate: latent positions X_ti (summary) : ', summary(AccRate[-c(1:4)]/(it-1)), '\n')

#-------------------------------------
## Analyze MCMC samples
#-------------------------------------
# ---------------------
# Posterior means: 
# ---------------------
E_alpha = mean(Alpha[-c(1:burnin)])
E_delta = mean(Delta[-c(1:burnin)])
E_gamma1_1 = mean(w[1,-c(1:burnin)])
E_gamma1_2 = mean(w[2,-c(1:burnin)])
E_gamma2 = mean(w[3,-c(1:burnin)])
E_t2 = mean(t2[-c(1:burnin)])
E_s2 = mean(s2[-c(1:burnin)])
cat('subwindow', subwindow, 'Posterior Mean (alpha, delta, gamma1_1, gamma1_2, gamma2, t, s): ', c(E_alpha, E_delta, E_gamma1_1, E_gamma1_2, E_gamma2, sqrt(E_t2), sqrt(E_s2)), "\n")
cat('subwindow', subwindow, 'Posterior Mean (alpha, delta, gamma1_1, gamma1_2, gamma2, t2, s2): ', c(E_alpha, E_delta, E_gamma1_1, E_gamma1_2, E_gamma2, E_t2, E_s2), "\n")  

# ---------------------
# Posterior distribution
# ---------------------
summary_table = NULL
param1 = c(mean(Alpha[-c(1:burnin)]), sd(Alpha[-c(1:burnin)]), quantile(Alpha[-c(1:burnin)], c(0.025, 0.975)))
param2 = c(mean(Delta[-c(1:burnin)]), sd(Delta[-c(1:burnin)]), quantile(Delta[-c(1:burnin)], c(0.025, 0.975)))
param3 = c(mean(w[1, -c(1:burnin)]), sd(w[1, -c(1:burnin)]), quantile(w[1, -c(1:burnin)], c(0.025, 0.975)))
param4 = c(mean(w[2, -c(1:burnin)]), sd(w[2, -c(1:burnin)]), quantile(w[2, -c(1:burnin)], c(0.025, 0.975)))
param5 = c(mean(w[3, -c(1:burnin)]), sd(w[3, -c(1:burnin)]), quantile(w[3, -c(1:burnin)], c(0.025, 0.975)))
cat('subwindow', subwindow, 'Alpha: mean, sd, CI_left, CI_right: ', param1, '\n')
cat('subwindow', subwindow, 'Delta: mean, sd, CI_left, CI_right: ', param2, '\n')
cat('subwindow', subwindow, 'gamma1_1: mean, sd, CI_left, CI_right: ', param3, '\n')
cat('subwindow', subwindow, 'gamma1_2: mean, sd, CI_left, CI_right: ', param4, '\n')
cat('subwindow', subwindow, 'gamma2: mean, sd, CI_left, CI_right: ', param5, '\n')

# ---------------------
# AUC
# ---------------------
N1 = Ns-burnin
EX = array(0,dim=c(p,TT,n))
for(it1 in c(1:N1)){
  EX = EX + X[[burnin+it1]]
}
EX = EX/N1
predY = array(0,dim=dim(Y))
for(tt in 1:TT){
  for(i in 1:n){
    for(j in c(1:n)[-i]){
      dijt = sqrt(t(EX[,tt,i]-EX[,tt,j])%*%(EX[,tt,i]-EX[,tt,j]))
      if (tt == 1){
        expon = E_alpha - dijt
      }else{
        expon = E_alpha + E_delta*Y[i,j,tt-1] - dijt
      }
      predY[i,j,tt] = 1/(1+exp(-expon))
    }
  }
}
ROC_pred = numeric(0)
ROC_Y = numeric(0)
if(subwindow==1){
  for(tt in 1:TT){
    ROC_pred = c(ROC_pred, c(predY[,,tt][upper.tri(predY[,,tt])]))
    ROC_Y = c(ROC_Y, c(Y[,,tt][upper.tri(Y[,,tt])]))
  }
}else{
  for(tt in 2:TT){ #for fitting with subwindow>1, taking off the first one, tt from 2
    ROC_pred = c(ROC_pred, c(predY[,,tt][upper.tri(predY[,,tt])]))
    ROC_Y = c(ROC_Y, c(Y[,,tt][upper.tri(Y[,,tt])]))
  }
}
pred = prediction(predictions=ROC_pred, labels=ROC_Y)
AUCPMean = slot(performance(pred, "auc"),"y.values")[[1]]  
# print(AUCPMean)
cat('subwindow', subwindow, 'AUC = ', AUCPMean, '\n')

# ---------------------
# AUC over time
# ---------------------
predY = array(0,dim=dim(Y))
for(tt in 1:TT){
  for(i in 1:n){
    for(j in c(1:n)[-i]){
      dijt = sqrt(t(EX[,tt,i]-EX[,tt,j])%*%(EX[,tt,i]-EX[,tt,j]))
      if (tt == 1){
        expon = E_alpha - dijt
      }else{
        expon = E_alpha + E_delta*Y[i,j,tt-1] - dijt
      }
      predY[i,j,tt] = 1/(1+exp(-expon))
    }
  }
}
AUC_t = numeric(TT)
for(tt in 1:TT){
  ROC_pred = c(predY[,,tt][upper.tri(predY[,,tt])])
  ROC_Y = c(Y[,,tt][upper.tri(Y[,,tt])])
  pred = prediction(predictions=ROC_pred, labels=ROC_Y)
  AUC_t[tt] = slot(performance(pred, "auc"),"y.values")[[1]]
}
cat('subwindow', subwindow, 'auc over time: ', AUC_t,'\n')

# ---------------------
# Calculate DIC
# --------------------
Ds = numeric(Ns-burnin)
for(it in 1:(Ns-burnin)){
  Ds[it] = c.get_loglik(Y, Alpha[it+burnin], Delta[it+burnin], X[[it+burnin]], c(n,p,TT), subwindow)[[1]]
}
D_bar = -2*mean(Ds)
D = -2*c.get_loglik(Y, E_alpha, E_delta, EX, c(n,p,TT), subwindow)[[1]]
cat('subwindow', subwindow, 'burnin:', burnin, ',', burnin/Ns, '\n')
cat('subwindow', subwindow, 'Expected deviance D_bar:', D_bar, '\n')
cat('subwindow', subwindow, 'Fitted deviance D(theta_hat):', D, '\n')
cat('subwindow', subwindow, 'Effective params p_D (=D_bar-D(theta_hat)):', D_bar-D, '\n')
cat('subwindow', subwindow, 'DIC (=D_bar+p_D):', 2*D_bar-D, '\n')
DIC_temp = 2*D_bar-D

# # save image
# if(save_image_yes){
#   out_name = paste('homedirec/sw_', subwindow, '_', save_image, sep='')
#   save.image(out_name)
#   cat('Image was saved in', out_name, '\n')
# }