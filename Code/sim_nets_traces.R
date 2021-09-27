library(igraph)
library(ggplot2)

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
      avg = mean(X[neighbors,,t-1]) # average for 1-d
    }
    A = avg - X[i,,t-1]
  }
  return(A)
}
# get synthetic data Y, which is a series of network observations over time, 
# of dimension nxnxT
# X: nxpxT, X[2,,t] represents the latent position (\in R^p) for actor 2 at time t 
# p: dimension of latent space
# tau: sd for the initial latent position, eg. X_it ~ N(0, tau2 I_p) (actor i from group 1 at time t)
# nu: initial latent position for group 2: adding nu to make the two groups start at distinctly different mean location 
# sigma: sd for the transition distribution, eg. X_it+1 ~ N(X_it + gamma*A(), sigma2 I_p)
# gamma: attraction parameter
# alpha: baseline level
# delta: edge persistance
# return: Y: nxpxT, X[,,t] represents the adjacency matrix at t
getData = function(n, TT, p, tau, sigma, gamma1, gamma2, alpha, delta, pi, nu){
  X = array(0, dim=c(n,p,TT))
  Y = array(0, dim=c(n,n,TT))
  P = array(0, dim=c(n,n,TT))
  # initialize latent position at time 1 by memberships
  # X[,,1] = matrix(rnorm(n*p, 0, tau), n, p)
  X[pi==1,,1] = matrix(rnorm(sum(pi==1)*p, 0, tau), sum(pi==1), p)
  X[pi==2,,1] = matrix(rnorm(sum(pi==2)*p, nu, tau), sum(pi==2), p)
  
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
      neighbors2 = (Y[,,t-1][i, ] == 1 & pi!=pi[i]) # neighbors for actor i
      A1 = getA(neighbors1, X, i, t, p)
      A2 = getA(neighbors2, X, i, t, p)
      #A0 = getA(neighbors1|neighbors2, X, i, t, p)
      # X[i,,t] = X[i,,t] + gammas[i]*A
      X[i,,t] = X[i,,t] + gamma1*A1 + gamma2*A2 
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


# plots
draw = function(data, ylim=c(-20, 20)){
  Z = data$X
  n = length(Z[,,1])
  title= bquote(atop(gamma^w==.(gamma1)~','~gamma^b==.(gamma2)~','
                     ~alpha==.(alpha)~','~delta==.(delta)~','
                     ~sigma==.(sigma)~','~tau==.(tau)))
  plot(1:TT, Z[1,,], type='l', ylim=ylim, ylab=expression('Z'[t][','][i]), 
       xlab='Time', col='lightsteelblue', lwd=2, cex.lab=1.2)
  title(title, line=1,  cex.main=1)
  for(i in 2:n){
    if(pi[i]==1){
      lines(1:TT, Z[i,,], type='l', col='lightsteelblue', lwd=2)
    }else{
      lines(1:TT, Z[i,,], type='l', col='indianred', lwd=2)
      cat('ddd')
    }
  }
}

draw_link_prob = function(data, timespot, pi, n){
  Node<-matrix (1:n, n, 1)
  Link = matrix (NA, n*(n-1)/2, 3) # all possible edges
  colnames(Link)[3]<-'weight' 
  colnames(Node)[1]<-'order'
  k=1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      Link[k, 1] = i
      Link[k, 2] = j
      Link[k, 3] = data$P[i, j, timespot]
      k=k+1
    }
  }
  net <- graph.data.frame(Link, Node, directed=F) 
  l <- layout.circle(net) 
  vertex.color = c('lightsteelblue', 'indianred ')[pi]
  E(net)$width <- 1+E(net)$weight/12
  plot(net, 
       edge.color=c ('black', 'blue', 'cyan')[((E(net)$weight<=2/3)+(E(net)$weight<=1/3) + 1)],
       layout=l, xlab=paste('t=', timespot), cex=2,
       vertex.color=vertex.color, vertex.size=15, vertex.label=NA) 
}

# Generate figures for flocking ===============================
n = 10
TT = 50
p = 1
pi = rep(c(1, 2), each=(n/2))
tau = 1
sigma = 1

nu = 10
alpha = 1
delta = 2
gamma1 = 0.3
gamma2 = 0.5

set.seed(1)
temp11 = getData(n, TT, p, tau, sigma, gamma1, gamma2, alpha, delta, pi, nu)
data = temp11

# traces
par(mfrow=c(1,1), cex.lab=2)
par(mar=c(5.1, 5.1, 4.1, 2.1))
draw(temp11, ylim=c(-4, 16)) #flocking

# nets
par(mfrow=c(2,3), cex.lab=2)
par(mar=c(5,1,1,0)+.1)
draw_link_prob(temp11, 1, pi, n)
draw_link_prob(temp11, 20, pi, n)
plot(NULL)
draw_link_prob(temp11, 30, pi, n)
draw_link_prob(temp11, 50, pi, n)
plot(NULL)
plot(legend(0,0.45,legend=c('p>2/3','1/3<p<2/3',' p<1/3'), col=c('black', 'blue','cyan'),
            lty=1,cex=2, title='Line'))

# Generate figures for polarization ===============================
n = 10
TT = 50
p = 1
pi = rep(c(1, 2), each=(n/2))
tau = 1
sigma = 1

nu = 0
tau = 1
alpha = 2
delta = 3 
gamma1 = 0.7
gamma2 = -0.1

set.seed(1) 
temp11 = getData(n, TT, p, tau, sigma, gamma1, gamma2, alpha, delta, pi, nu)
data = temp11

# traces
par(mfrow=c(1,1), cex.lab=2)
par(mar=c(5.1, 5.1, 4.1, 2.1))
draw(temp11, ylim=c(-10, 10)) # polarization

# nets
par(mfrow=c(2,3), cex.lab=2)
par(mar=c(5,1,1,0)+.1)
draw_link_prob(temp11, 1, pi, n)
draw_link_prob(temp11, 20, pi, n)
plot(NULL)
draw_link_prob(temp11, 30, pi, n)
draw_link_prob(temp11, 50, pi, n)
plot(NULL)
plot(legend(0,0.45,legend=c('p>2/3','1/3<p<2/3',' p<1/3'), col=c('black', 'blue','cyan'),
            lty=1,cex=2, title='Line'))
