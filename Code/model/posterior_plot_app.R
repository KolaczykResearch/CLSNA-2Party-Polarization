distance = function(x, y){
  return(sqrt(sum((x-y)^2)))
}

# ---------------------------
# -Figure 5 for twitter 
# ---------------------------
# load('twitter_run_yearly_2010_2020.RData') 
# load image obtained by running main_application.R

#### latent positions ----------------------------- 
par(mfrow=c(2,3))
xl = c(1.15*min(EX[1,,]),1.15*max(EX[1,,]))
yl = c(1.15*min(EX[2,,]),1.15*max(EX[2,,]))
lims = range(c(xl,yl))
for(t in c(1,3,5,7,9,11)){
  plot(t(EX[,t,pi==1]),xlim=lims,ylim=lims,xlab="",ylab="",
       pch=16,cex=0.75, main=paste(t+2009), col='dodgerblue')
  # points(mean(t(EX[,t,pi==1])[,1]), mean(t(EX[,t,pi==1])[,2]), pch=19)
  points(t(EX[,t,pi==2]),xlim=lims,ylim=lims,xlab="",ylab="",
         pch=16,cex=0.75,col='firebrick')
  # points(mean(t(EX[,t,pi==2])[,1]), mean(t(EX[,t,pi==2])[,2]), pch=19)
  # cat(sqrt(sum((apply(t(EX[,t,pi==1]),2,mean)-apply(t(EX[,t,pi==2]),2,mean))^2)), '\n')
}

#### MEAN distance- mean pairwise distance ----------------------------
## mean distances - the average of distances between all pairs of objects, 
# where each pair is made up of one object from each group. 
mean_distance_1 = numeric(TT)
mean_distance_2 = numeric(TT)
between_distance = numeric(TT)
cat('between-party distance:\n')
for(t in 1:TT){
  group1 = t(EX[,t,pi==1])
  temp = 0
  ct = 0
  for(i in 1:(nrow(group1)-1)){
    for(j in 2:nrow(group1)){
      temp = temp + distance(group1[i,], group1[j,])
      ct = ct+1
    }
  }
  mean_distance_1[t] = temp/ct
}
for(t in 1:TT){
  group2 = t(EX[,t,pi==2])
  temp = 0
  ct = 0
  for(i in 1:(nrow(group2)-1)){
    for(j in 2:nrow(group2)){
      temp = temp + distance(group2[i,], group2[j,])
      ct = ct+1
    }
  }
  mean_distance_2[t] = temp/ct
}
cat('within-party distance:\n')
for(t in 1:TT){
  group1 = t(EX[,t,pi==1])
  group2 = t(EX[,t,pi==2])
  temp = 0
  for(i in 1:nrow(group1)){
    for(j in 1:nrow(group2)){
      temp = temp + distance(group1[i,], group2[j,])
    }
  }
  between_distance[t]=temp/(nrow(group1)*nrow(group2))
}

par(mfrow=c(1,1))
plot(2010:2020, mean_distance_1, type='b', pch=16, col='dodgerblue', lwd=2, 
     xlab='Year', ylab='Mean distance', ylim=c(5,20))
lines(2010:2020, mean_distance_2, col='firebrick', lwd=2, type='b', pch=16)
lines(2010:2020, between_distance, lwd=2, type='b', pch=16,)
legend('topright', c('within Democrats', 'within Republicans', 'between two groups'), 
       col=c('dodgerblue', 'firebrick', 'black'), lty=c(1,1,1), lwd=c(2,2,2), pch=c(16,16,16), cex=0.9)

# ---------------------------
# -Figure 7 for twitter 
# ---------------------------
# load('reddit_political_year.RData') 
# load image obtained by running main_application.R

#### latent positions ----------------------------- 
# Posterior mean latent positions showing temporal dynamics
par(mfrow=c(1,5))
xl = c(1.15*min(EX[1,,]),1.15*max(EX[1,,]))
yl = c(1.15*min(EX[2,,]),1.15*max(EX[2,,]))
lims = range(c(xl,yl))
for(t in 1:TT){
  plot(t(EX[,t,pi==1]),xlim=lims,ylim=lims,xlab="",ylab="",
       pch=16,cex=0.75, main=paste(t+2014, t+2015, sep='-'), col='dodgerblue')
  points(t(EX[,t,pi==2]),xlim=lims,ylim=lims,xlab="",ylab="",
         pch=16,cex=0.75,col='firebrick')
  points(mean(t(EX[,t,pi==1])[,1]), mean(t(EX[,t,pi==1])[,2]), pch=20, col='green')
  points(mean(t(EX[,t,pi==2])[,1]), mean(t(EX[,t,pi==2])[,2]), pch=20, col='green')
}

#### MEAN distance- mean pairwise distance ----------------------------
## mean distances - the average of distances between all pairs of objects, 
# where each pair is made up of one object from each group. 
mean_distance_1 = numeric(TT)
mean_distance_2 = numeric(TT)
between_distance = numeric(TT)
cat('between-party distance:\n')
for(t in 1:TT){
  group1 = t(EX[,t,pi==1])
  temp = 0
  ct = 0
  for(i in 1:(nrow(group1)-1)){
    for(j in 2:nrow(group1)){
      temp = temp + distance(group1[i,], group1[j,])
      ct = ct+1
    }
  }
  mean_distance_1[t] = temp/ct
}
for(t in 1:TT){
  group2 = t(EX[,t,pi==2])
  temp = 0
  ct = 0
  for(i in 1:(nrow(group2)-1)){
    for(j in 2:nrow(group2)){
      temp = temp + distance(group2[i,], group2[j,])
      ct = ct+1
    }
  }
  mean_distance_2[t] = temp/ct
}
cat('within-party distance:\n')
for(t in 1:TT){
  group1 = t(EX[,t,pi==1])
  group2 = t(EX[,t,pi==2])
  temp = 0
  for(i in 1:nrow(group1)){
    for(j in 1:nrow(group2)){
      temp = temp + distance(group1[i,], group2[j,])
    }
  }
  between_distance[t]=temp/(nrow(group1)*nrow(group2))
}

par(mfrow=c(1,1))
plot(1:5, mean_distance_1, type='b', pch=16, col='dodgerblue', lwd=2, 
     xlab='Year', ylab='Mean distance', ylim=c(0,4.5), xaxt = "n")
lines(1:5, mean_distance_2, col='firebrick', lwd=2, type='b', pch=16)
lines(1:5, between_distance, lwd=2, type='b', pch=16,)
legend('bottomright', c('within Democrats', 'within Republicans', 'between two groups'), 
       col=c('dodgerblue', 'firebrick', 'black'), lty=c(1,1,1), lwd=c(2,2,2), pch=c(16,16,16), cex=0.9)
axis(1, at=1:5, labels=c('2015-2016', '2016-2017', '2017-2018', '2018-2019', '2019-2020'))


