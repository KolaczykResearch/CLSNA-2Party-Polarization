library(igraph)
library(ggplot2)

# ------------------------------------------------------------------------------------------------------------
# # twitter density plot 
# -------------------------------------------------------------------------------------------------------------
load('../Data/nets_twitter.RData')
TT = dim(Y)[3]
all = apply(Y, 3, FUN = function(Yt){
  mean(Yt[lower.tri(Yt)]) 
})
within_group1 = apply(Y[pi=='D', pi=='D',], 3, FUN = function(Yt){
  mean(Yt[lower.tri(Yt)])
})
within_group2 = apply(Y[pi=='R', pi=='R',], 3, FUN = function(Yt){
  mean(Yt[lower.tri(Yt)])
})
between_groups = apply(Y[pi=='D', pi=='R',], 3, mean)
density = c(all, within_group1, within_group2, between_groups)
type = rep(c('Among all members', 'Among Democrats', 'Among Republicans ', 'Between two groups'), each=TT)
t = rep(2010:2020, times=4)
t = as.factor(t)
gg_data = data.frame(Density=density, Type=type, t=t)
p1 = ggplot(gg_data, aes(x=t, y=Density, group=Type, color=Type)) + 
  geom_line() + ylim(0, 0.8) + theme_bw() + 
  scale_color_manual(values=c("black", "dodgerblue", "firebrick", "darkgreen")) + 
  xlab('Year') + theme(legend.position=c(0.2, 0.8)) + 
  theme(legend.text=element_text(size=11))
p1

# ------------------------------------------------------------------------------------------------------------
# # reddit density plot 
# -------------------------------------------------------------------------------------------------------------
load('../Data/nets_reddit.RData')
TT = dim(Y)[3]
all = apply(Y, 3, FUN = function(Yt){
  mean(Yt[lower.tri(Yt)]) 
})
within_group1 = apply(Y[pi=='D', pi=='D',], 3, FUN = function(Yt){
  mean(Yt[lower.tri(Yt)])
})
within_group2 = apply(Y[pi=='R', pi=='R',], 3, FUN = function(Yt){
  mean(Yt[lower.tri(Yt)])
})
between_groups = apply(Y[pi=='D', pi=='R',], 3, mean)
density = c(all, within_group1, within_group2, between_groups)
type = rep(c('Among all members', 'Among Democrats', 'Among Republicans ', 'Between two groups'), each=TT)
t = rep(1:5, times=4)
t = as.factor(t)
gg_data = data.frame(Density=density, Type=type, t=t)
p1 = ggplot(gg_data, aes(x=t, y=Density, group=Type, color=Type)) + 
  geom_line() + ylim(0, 1) + theme_bw() + 
  scale_color_manual(values=c("black", "dodgerblue", "firebrick", "darkgreen")) + 
  xlab('One-year period')  + theme(legend.position=c(0.2, 0.2)) + 
  scale_x_discrete(labels=c("1" = "2015-2016", "2" = "2016-2017",
                            "3" = "2017-2018", "4" = "2018-2019", "5" = "2019-2020")) + 
  theme(legend.text=element_text(size=11))
p1
