############################################################
###   Get Synthetic Data from Model with no change point ###
############################################################
sigmoid = function(x){
  1/(1 + exp(-x))
}
getA = function(neighbors, X, i, t, p){
  if(sum(neighbors)==0){
    A = 0
  }else if(sum(neighbors)==1){
    A = X[neighbors,,t-1] - X[i,,t-1]
  }else{
    if (p>1){
      avg = apply(X[neighbors,,t-1], 2, mean) # average of latent position of neighbors
    }else{
      avg = mean(X[neighbors,,t-1]) # average of latent position of neighbors
    }
    A = avg - X[i,,t-1]
  }
  return(A)
}
# Get synthetic data Y, which is a series of network observations over time, 
# Args:
# n: number of vertices in network
# TT: number of time points
# p: dimension of latent space
# tau: sd for the initial latent position, eg. X_it ~ N(0, tau2 I_p) (actor i at time t)
# sigma: sd for the transition distribution, eg. X_it+1 ~ N(X_it + gamma1*A1() + gamma2*A2(), sigma2 I_p)
# gamma1: within-group attraction parameter, a vector
# gamma2: between-group attraction parameter
# alpha: baseline level
# delta: edge persistence
# pi: membership of each node i
# Return: 
# X: nxpxT, X[2,,t] represents the latent position (\in R^p) for actor 2 at time t 
# Y: nxnxT, Y[,,t] represents the adjacency matrix at t
# P: nxnxT, P[i,j,t] represents link prob. between i, j at t
getData2 = function(n, TT, p, tau, sigma, gamma1, gamma2, alpha, delta, pi){
  X = array(0, dim=c(n,p,TT))
  Y = array(0, dim=c(n,n,TT))
  P = array(0, dim=c(n,n,TT))
  # initialize latent position at time 1
  X[,,1] = matrix(rnorm(n*p, 0, tau), n, p)
  # initialize random network at time 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      s = sqrt(sum((X[i,,1] - X[j,,1])^2)) # euclidean distance between i, j at t
      pp = sigmoid(alpha - s)
      P[i,j,1] = pp
      P[j,i,1] = P[i,j,1]
      Y[i,j,1] = sample(c(1, 0), 1, prob=c(pp, 1-pp))
      Y[j,i,1] = Y[i,j,1]
    }
  }
  for(t in 2:TT){
    # update latent positions X[,,t]
    X[,,t] = X[,,t-1] + matrix(rnorm(n*p, 0, sigma), n, p)
    for(i in 1:n){
      # add attraction effect
      neighbors1 = (Y[,,t-1][i, ] == 1 & pi==pi[i]) # neighbors for actor i with same membership 
      neighbors2 = (Y[,,t-1][i, ] == 1 & pi!=pi[i]) # neighbors for actor i with different membership
      A1 = getA(neighbors1, X, i, t, p)
      A2 = getA(neighbors2, X, i, t, p)
      X[i,,t] = X[i,,t] + gamma1[pi[i]]*A1 + gamma2*A2
    }
    # update observed network Y[,,t]
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        s = sqrt(sum((X[i,,t] - X[j,,t])^2)) # euclidean distance between i, j at t
        pp = sigmoid(alpha + delta*Y[i,j,t-1] - s)
        P[i,j,t] = pp
        P[j,i,t] = P[i,j,t]
        Y[i,j,t] = sample(c(1, 0), 1, prob=c(pp, 1-pp))
        Y[j,i,t] = Y[i,j,t]
      }
    }
  }
  return(list(X=X, Y=Y, P=P))
}

# # ----------
# # example
# # ----------
# n = 20
# TT = 10
# p = 2
# pi = rep(c(1, 2), each=(n/2))
# tau = 1
# sigma = 0.01
# alpha = 1
# delta = 3
# gamma1 = c(0.7, 0.4)
# gamma1 = c(0.7, -0.2)
# gamma2 = -0.5
# set.seed(1122)
# temp11 = getData2(n, TT, p, tau, sigma, gamma1, gamma2, alpha, delta, pi)
# Y = temp11$Y
# Z = temp11$X
# P = temp11$P

