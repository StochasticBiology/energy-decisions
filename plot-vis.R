library(ggplot2)
library(gridExtra)

# time series plots

df.0 = read.csv("time-series-0.csv"); df.0$expt = 0
df.1 = read.csv("time-series-1.csv"); df.1$expt = 1
df.2 = read.csv("time-series-2.csv"); df.2$expt = 2
df = rbind(df.0, df.1, df.2)
df.long = melt(df, id.vars=c("t", "expt"))
ggplot(df.long[df.long$variable %in% c("pro_1", "prooff_1", "rna_1", "p_1", "pp_1"),], 
       aes(x=t, y=value, color=factor(expt))) + geom_line() + 
  scale_x_log10() +
  facet_wrap(~variable, scales="free")

# stochastic simulation runs

df.0 = read.csv("gillespie-series-0.csv"); df.0$expt = 0
df.1 = read.csv("gillespie-series-1.csv"); df.1$expt = 1
df.2 = read.csv("gillespie-series-2.csv"); df.2$expt = 2
df.3 = read.csv("gillespie-series-3.csv"); df.3$expt = 3
df = rbind(df.0, df.1, df.2, df.3)
df.long = melt(df, id.vars=c("t", "expt"))
ggplot(df.long[df.long$variable %in% c("p_1", "p_2"),], 
       aes(x=t, y=value, color=factor(variable))) + geom_line() + 
  facet_wrap(~expt, scales="free", ncol=1)

sdf.0 = read.csv("gillespie-switches-0.csv"); sdf.0$expt = 0
sdf.1 = read.csv("gillespie-switches-1.csv"); sdf.1$expt = 1
sdf.2 = read.csv("gillespie-switches-2.csv"); sdf.2$expt = 2
sdf.3 = read.csv("gillespie-switches-3.csv"); sdf.3$expt = 3
sdf = rbind(sdf.0, sdf.1, sdf.2, sdf.3)

ggplot(sdf, aes(x=log10(dt), fill=factor(expt))) + geom_histogram(aes(y = ..density..), position="dodge")

ggplot(sdf, aes(x = factor(expt), y=log10(dt))) + geom_boxplot()

# attractor basins from parameter scan

df = read.csv("param-scan.csv")

g.1 = ggplot(df[df$n==1,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param~scale) + coord_cartesian(xlim=c(0,30),ylim=c(0,30))
g.2 = ggplot(df[df$n==2,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param~scale) + coord_cartesian(xlim=c(0,30),ylim=c(0,30))
g.3 = ggplot(df[df$n==3,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param~scale) + coord_cartesian(xlim=c(0,30),ylim=c(0,30))
g.4 = ggplot(df[df$n==4,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param~scale) + coord_cartesian(xlim=c(0,30),ylim=c(0,30))

grid.arrange(g.1, g.2, g.3, g.4, nrow=2)
png("set-2023.png", width=1000, height=2000)
grid.arrange(g.1, g.2, g.3, g.4, nrow=2)
dev.off()

## focussed plots

df = read.csv("zoom-scan.csv")

zg.1 = ggplot(df[df$param==0,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(ATP~gamma3) + coord_cartesian(xlim=c(0,30),ylim=c(0,30))
zg.2 = ggplot(df[df$param==1,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(ATP~cd2) + coord_cartesian(xlim=c(0,30),ylim=c(0,30))

grid.arrange(zg.1, zg.2, nrow=1)
