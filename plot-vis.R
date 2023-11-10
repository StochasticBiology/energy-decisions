library(ggplot2)
library(ggpubr)
library(reshape)

sf = 3

# time series plots

df.0 = read.csv("time-series-0.csv"); df.0$expt = 0
df.1 = read.csv("time-series-1.csv"); df.1$expt = 1
df.2 = read.csv("time-series-2.csv"); df.2$expt = 2
df = rbind(df.0, df.1, df.2)
df.long = melt(df, id.vars=c("t", "expt"))
df.long$Experiment = "Original, ATP = 1"
df.long$Experiment[df.long$expt==1] = "Scaled cd2, ATP = 2"
df.long$Experiment[df.long$expt==2] = "Scaled gamma3, ATP = 2"
ggplot(df.long[df.long$variable %in% c("pro_1", "prooff_1", "rna_1", "p_1", "pp_1"),], 
       aes(x=t, y=value, color=Experiment)) + geom_line() + 
  scale_x_log10() + facet_wrap(~variable, scales="free") + theme_classic()

png("time-series.png", width=600*sf, height=300*sf, res=72*sf)
ggplot(df.long[df.long$variable %in% c("pro_1", "prooff_1", "rna_1", "p_1", "pp_1"),], 
       aes(x=t, y=value, color=Experiment)) + geom_line() + 
  scale_x_log10() + facet_wrap(~variable, scales="free") + theme_classic()
dev.off()

# stochastic simulation runs

df.0 = read.csv("gillespie-series-0.csv"); df.0$expt = 0
df.1 = read.csv("gillespie-series-1.csv"); df.1$expt = 1
df.2 = read.csv("gillespie-series-2.csv"); df.2$expt = 2
df.3 = read.csv("gillespie-series-3.csv"); df.3$expt = 3
df = rbind(df.0, df.1, df.2, df.3)
df.long = melt(df, id.vars=c("t", "expt"))
df.long$Experiment = "Original"
df.long$Experiment[df.long$expt==1] = "Scaled, ATP = 1"
df.long$Experiment[df.long$expt==2] = "Scaled, ATP = 0.25"
df.long$Experiment[df.long$expt==3] = "Scaled, ATP = 2"
gill.traces = ggplot(df.long[df.long$variable %in% c("p_1", "p_2"),], 
       aes(x=t, y=value, color=factor(variable))) + geom_line() + labs(color="Protein") +
  ylab("Count") +
  facet_wrap(~Experiment, scales="free", ncol=1)

ggplot(df.long[df.long$variable %in% c("p_1", "p_2"),], 
       aes(x=t, y=value, color=factor(variable))) + geom_line() + labs(color="Protein") +
  ylab("Count") +
  facet_wrap(~Experiment, scales="free", ncol=1) + theme_classic()
dev.off()

df$Experiment = "Original"
df$Experiment[df$expt==1] = "Scaled, ATP = 1"
df$Experiment[df$expt==2] = "Scaled, ATP = 0.25"
df$Experiment[df$expt==3] = "Scaled, ATP = 2"
gill.space = ggplot(df, aes(x=p_1, y=p_2, color=Experiment)) + geom_path(alpha=0.2) +
  facet_wrap(~Experiment, nrow=2) + theme_classic()  +
  theme(legend.position = "none")

sdf.0 = read.csv("gillespie-switches-0.csv"); sdf.0$expt = 0
sdf.1 = read.csv("gillespie-switches-1.csv"); sdf.1$expt = 1
sdf.2 = read.csv("gillespie-switches-2.csv"); sdf.2$expt = 2
sdf.3 = read.csv("gillespie-switches-3.csv"); sdf.3$expt = 3
sdf = rbind(sdf.0, sdf.1, sdf.2, sdf.3)

sdf$Experiment = "Original"
sdf$Experiment[sdf$expt==1] = "Scaled, ATP = 1"
sdf$Experiment[sdf$expt==2] = "Scaled, ATP = 0.25"
sdf$Experiment[sdf$expt==3] = "Scaled, ATP = 2"
gill.switches = ggplot(sdf, aes(x=log10(dt), fill=Experiment)) + 
  geom_histogram(aes(y = ..density..) , alpha=0.4,position="identity") + 
  theme_light() + facet_wrap(~Experiment, nrow=4) + xlab("log10(switch time)") +
  theme(legend.position = "none")

png("stoch.png", width=800*sf, height=600*sf, res=72*sf)
ggarrange(gill.traces, 
          ggarrange(gill.space, gill.switches, labels=c("B", "C"), nrow=1),
          labels=c("A", ""), nrow=2)
dev.off()
#ggplot(sdf, aes(x = factor(expt), y=log10(dt))) + geom_boxplot()

# attractor basins from parameter scan

df = read.csv("param-scan.csv")

var.names = c("gamma2", "gamma3", "gamma4", "cd1", "cd2", "ca1", "ca2", "lambda2", "lambda3")
scale.names = c("x0.01", "x0.1", "x1", "x10", "x100")
df$param.name = var.names[df$param+1]
df$scale.name = scale.names[df$scale+1]

g.1 = ggplot(df[df$n==1,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~scale.name) + coord_cartesian(xlim=c(0,30),ylim=c(0,30))
g.2 = ggplot(df[df$n==2,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~scale.name) + coord_cartesian(xlim=c(0,30),ylim=c(0,30))
g.3 = ggplot(df[df$n==3,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~scale.name) + coord_cartesian(xlim=c(0,30),ylim=c(0,30))
g.4 = ggplot(df[df$n==4,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~scale.name) + coord_cartesian(xlim=c(0,30),ylim=c(0,30))

ggarrange(g.1, g.2, g.3, g.4, labels=c("A", "B", "C", "D"), nrow=2, ncol=2)
png("param-scan.png", width=800*sf, height=1200*sf, res=72*sf)
ggarrange(g.1, g.2, g.3, g.4, labels=c("A", "B", "C", "D"), nrow=2, ncol=2)
dev.off()

png("param-scan-2.png", width=600*sf, height=600*sf, res=72*sf)
ggarrange(g.2)
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

ggarrange(zg.1, zg.2, nrow=1, ncol=2)
png("zoom-scan.png", width=800*sf, height=400*sf, res=72*sf)
ggarrange(zg.1, zg.2, labels=c("A", "B"), nrow=1, ncol=2)
dev.off()

# attractor basins from parameter scan with expanded coords

df = read.csv("param-scan-300.csv")

var.names = c("gamma2", "gamma3", "gamma4", "cd1", "cd2", "ca1", "ca2", "lambda2", "lambda3")
scale.names = c("x0.01", "x0.1", "x1", "x10", "x100")
df$param.name = var.names[df$param+1]
df$scale.name = scale.names[df$scale+1]

g.1 = ggplot(df[df$n==1,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~scale.name) + coord_cartesian(xlim=c(0,300),ylim=c(0,300))
g.2 = ggplot(df[df$n==2,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~scale.name) + coord_cartesian(xlim=c(0,300),ylim=c(0,300))
g.3 = ggplot(df[df$n==3,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~scale.name) + coord_cartesian(xlim=c(0,300),ylim=c(0,300))
g.4 = ggplot(df[df$n==4,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~scale.name) + coord_cartesian(xlim=c(0,300),ylim=c(0,300))

ggarrange(g.1, g.2, g.3, g.4, nrow=2, ncol=2)
png("param-scan-300.png", width=800*sf, height=1200*sf, res=72*sf)
ggarrange(g.1, g.2, g.3, g.4, labels=c("A", "B", "C", "D"), nrow=2, ncol=2)
dev.off()

