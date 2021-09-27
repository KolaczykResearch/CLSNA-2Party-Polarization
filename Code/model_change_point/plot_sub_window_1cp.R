library('ggplot2')

####################
# TWITTER
###################
data = to_plot_yearly_evolution_data_1cp # obtained by running main_application_mc.R
p_total = 2 # number of time periods

# sort data by the first column (year)
data = data[order(data[,1]), ]
colnames(data) = c('Year', 'alpha', 'alpha_se', 'alpha_l', 'alpha_r', 
                  'delta', 'delta_se', 'delta_l', 'delta_r',
                  'gamma11', 'gamma11_se', 'gamma11_l', 'gamma11_r',
                  'gamma12', 'gamma12_se', 'gamma12_l', 'gamma12_r',
                  'gamma2', 'gamma2_se', 'gamma2_l', 'gamma2_r')
data

# -------------
# plot for delta
# -------------
# bar plots
data=data.frame(data)
pd = position_dodge(0.01) # move them .01 to the left and right
data$Year=as.factor(data$Year)
p1 = ggplot(data, aes(x=Year, y=delta)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin=delta_l, ymax=delta_r), colour="black", width=.2, position=pd) +
  ylab(bquote('Edge persistence')) +
  ggtitle(bquote('Evolution of edge persistence'~delta)) +
  scale_x_discrete(labels=c('1' = "2010-2014", '2' = "2015-2020")) + 
  theme_bw() 
p1

# -----------------------------------------
# Plot for gamma^w_1, gamma^w_2, gamma^w_3
# -----------------------------------------
mean = c(data[,'gamma11'], data[,'gamma12'], data[,'gamma2'])
ci_l = c(data[,'gamma11_l'], data[,'gamma12_l'], data[,'gamma2_l'])
ci_r = c(data[,'gamma11_r'], data[,'gamma12_r'], data[,'gamma2_r'])
Year = rep(1:p_total, 3)
Party = rep(c('Within Democrats', 'Within Republicans', 'Between the two'), each=p_total)

df=data.frame(mean, ci_l, ci_r, Year, Party)
df$Year = as.factor(df$Year)
df$Party = factor(df$Party, levels=c('Within Democrats', 'Within Republicans', 'Between the two'))

p3 = ggplot(df, aes(x=Year, y=mean, color=Party, fill=Party)) +  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=ci_l, ymax=ci_r), colour="black", width=.2, position=position_dodge(.9)) +
  ylab('Attraction/repulsion') +
  ggtitle(bquote('Evolution of attraction/repulsion, '~gamma[1]^w~', '~gamma[2]^w~', '~gamma^b)) +
  scale_x_discrete(labels=c('1' = "2010-2014", '2' = "2015-2020")) + 
  theme_bw() + 
  scale_colour_manual(values=c('Within Democrats'='#619cff', 'Within Republicans'='#F8766D', 'Between the two'='#E69F00')) + 
  scale_fill_manual(values=c('Within Democrats'='#619cff', 'Within Republicans'='#F8766D', 'Between the two'='#E69F00')) + 
  theme(legend.position='bottom', legend.title=element_blank(), legend.text=element_text(size=11))
p3


####################
# REDDIT
###################
data = to_plot_yearly_evolution_data_1cp # obtained by running main_application_mc.R
p_total = 2

# sort data by the first column (year)
data=data.frame(data)
data = data[order(data[,1]), ]
colnames(data) = c('Year', 'alpha', 'alpha_se', 'alpha_l', 'alpha_r', 
                   'delta', 'delta_se', 'delta_l', 'delta_r',
                   'gamma11', 'gamma11_se', 'gamma11_l', 'gamma11_r',
                   'gamma12', 'gamma12_se', 'gamma12_l', 'gamma12_r',
                   'gamma2', 'gamma2_se', 'gamma2_l', 'gamma2_r')
data

# -------------
# plot for delta
# -------------
# bar plots
pd = position_dodge(0.01) # move them .01 to the left and right
data$Year=as.factor(data$Year)
p1 = ggplot(data, aes(x=Year, y=delta)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin=delta_l, ymax=delta_r), colour="black", width=.2, position=pd) +
  ylab(bquote('Edge persistence')) +
  ggtitle(bquote('Evolution of edge persistence'~delta)) +
  scale_x_discrete(labels=c('1' = "April 2015 - March 2018", '2' = "April 2018- March 2020")) + 
  theme_bw() 
p1

# -----------------------------------------
# Plot for gamma^w_1, gamma^w_2, gamma^w_3
# -----------------------------------------
mean = c(data[,'gamma11'], data[,'gamma12'], data[,'gamma2'])
ci_l = c(data[,'gamma11_l'], data[,'gamma12_l'], data[,'gamma2_l'])
ci_r = c(data[,'gamma11_r'], data[,'gamma12_r'], data[,'gamma2_r'])
Year = rep(1:p_total, 3)
Party = rep(c('Within Democrats', 'Within Republicans', 'Between the two'), each=p_total)

df=data.frame(mean, ci_l, ci_r, Year, Party)
df$Year = as.factor(df$Year)
df$Party = factor(df$Party, levels=c('Within Democrats', 'Within Republicans', 'Between the two'))

# pd = position_dodge(0.1) # move them .01 to the left and right
p3 = ggplot(df, aes(x=Year, y=mean, color=Party, fill=Party)) +  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=ci_l, ymax=ci_r), colour="black", width=.2, position=position_dodge(.9)) +
  ylab('Attraction/repulsion') +
  ggtitle(bquote('Evolution of attraction/repulsion, '~gamma[1]^w~', '~gamma[2]^w~', '~gamma^b)) +
  scale_x_discrete(labels=c('1' = "April 2015 - March 2018", '2' = "April 2018- March 2020")) + 
  theme_bw() + 
  scale_colour_manual(values=c('Within Democrats'='#619cff', 'Within Republicans'='#F8766D', 'Between the two'='#E69F00')) + 
  scale_fill_manual(values=c('Within Democrats'='#619cff', 'Within Republicans'='#F8766D', 'Between the two'='#E69F00')) + 
  theme(legend.position='bottom', legend.title=element_blank(), legend.text=element_text(size=11))
p3



