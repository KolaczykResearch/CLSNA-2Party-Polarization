# ----------------------------------------------
# Posterior mean latent positions showing temporal dynamics
# ----------------------------------------------
par(mfrow=c(2,5))
xl = c(1.15*min(EX[1,,]),1.15*max(EX[1,,]))
yl = c(1.15*min(EX[2,,]),1.15*max(EX[2,,]))
lims = range(c(xl,yl))
xl = c(1.15*min(Z_true[,1,]),1.15*max(Z_true[,1,]))
yl = c(1.15*min(Z_true[,2,]),1.15*max(Z_true[,2,]))
lims = range(c(xl,yl))
for(t in 1:TT){
  plot(t(EX[,t,pi==1]),xlim=lims,ylim=lims,xlab="",ylab="",
       pch=16,cex=0.75, main=paste('t =', t))
  points(t(EX[,t,pi==2]),xlim=lims,ylim=lims,xlab="",ylab="",
         pch=16,cex=0.75,col='red')
}

## ----------------
# trace plots of parameter samples
## ---------------
par(mfrow=c(2,3))
plot(Alpha[burnin:Ns],type="l",ylab=expression(alpha),
     xlab="",cex.lab=1,cex.axis=1)
plot(Delta[burnin:Ns],type="l",ylab=expression(delta),
     xlab="",cex.lab=1,cex.axis=1)
plot(w[1, burnin:Ns],type="l",ylab=expression(gamma^w[1]),
     xlab="",cex.lab=1,cex.axis=1)
plot(w[2, burnin:Ns],type="l",ylab=expression(gamma^w[2]),
     xlab="",cex.lab=1,cex.axis=1)
plot(w[3, burnin:Ns],type="l",ylab=expression(gamma^b),
     xlab="",cex.lab=1,cex.axis=1)
plot(t2[burnin:Ns],type="l",ylab=expression(tau^2),
     xlab="",cex.lab=1,cex.axis=1)

## ----------------
# trace plots of latent positions
## ---------------
x_trace=numeric(Ns)
for(i in 1:Ns){
  agent=4
  t = 10
  x_trace[i] = X[[i]][1,t,agent]
}
plot(x_trace[burnin:Ns],type="l",ylab=expression(x[it]),
     xlab="",cex.lab=1,cex.axis=1)






