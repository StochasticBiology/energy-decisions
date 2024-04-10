library(ggplot2)
library(ggpubr)

theme_set(theme_light())
sf = 2

# labels for the various factor parameters involved
scale.labels = c("x0.01", "x0.1", "x1", "x10", "x100")
scale.labels = factor(scale.labels, levels=scale.labels)
param.labels = c("gamma2", "gamma3", "gamma4", "cd1", "cd2", "ca1", "ca2", "lambda2", "lambda3")
param.labels = factor(param.labels, levels=param.labels)
method.labels = c("RK4", "Euler")
method.labels = factor(method.labels, levels=method.labels)
IC.labels = c("0,0", "10,0", "0,10", "10,10")
IC.labels = factor(IC.labels, levels=IC.labels)

########### Experiment 0 -- compare Euler and RK4 methods
### manuscript: SI Fig XXX

df.0 = read.csv("grn-sim-0.csv")

g.0.1 = ggplot(df.0, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() + theme(legend.position="none") +
  facet_grid(~method.labels[euler+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

df.0t = read.csv("grn-sim-0-t.csv")

tx = 0
ty = 9

g.0.2 = ggplot() + geom_line(data=df.0t[df.0t$euler == 0 & df.0t$ip1==tx & df.0t$ip2==ty,], aes(x=t,y=p1), color="blue") +
  geom_line(data=df.0t[df.0t$euler == 1 & df.0t$ip1==tx & df.0t$ip2==ty,], aes(x=t,y=p1), color="blue", size=4, alpha =0.2) 

png("expt-0.png", width=600*sf, height=300*sf, res=72*sf)
print( ggarrange(g.0.1, g.0.2, widths = c(2,1), labels = c("A", "B")) )
dev.off()

############ Experiment 1 -- parameter scans for default model
### manuscript: Fig 2 (n = 2), Supp Fig XXX (all n)

df.1.1 = read.csv("grn-sim-1.1.csv")
g.1.1 = ggplot(df.1.1, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(param.labels[param+1] ~ scale.labels[scale+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

df.1.2 = read.csv("grn-sim-1.2.csv")
g.1.2 = ggplot(df.1.2, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(param.labels[param+1] ~ scale.labels[scale+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

df.1.3 = read.csv("grn-sim-1.3.csv")
g.1.3 = ggplot(df.1.3, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(param.labels[param+1] ~ scale.labels[scale+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))
 
df.1.4 = read.csv("grn-sim-1.4.csv")
g.1.4 = ggplot(df.1.4, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(param.labels[param+1] ~ scale.labels[scale+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

png("expt-1.1.png", width=600*sf, height=600*sf, res=72*sf)
print( g.1.2 )
dev.off()

png("expt-1.2.png", width=1200*sf, height=1200*sf, res=72*sf)
print( ggarrange(g.1.1, g.1.2, g.1.3, g.1.4, nrow=2, ncol=2, labels=c("A", "B", "C", "D")) )
dev.off()

############ Experiment 2 -- parameter scans for on-DNA dimerisation
### manuscript: Supp Fig XXX 

df.2.2 = read.csv("grn-sim-2.2.csv")
g.2.2 = ggplot(df.2.2, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(param.labels[param+1] ~ scale.labels[scale+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

png("expt-2.png", width=600*sf, height=600*sf, res=72*sf)
print( g.2.2 )
dev.off()

############ Experiment 3 -- ATP x gamma3 influence
### manuscript: Fig 3A

df.3 = read.csv("grn-sim-3.csv")
g.3 = ggplot(df.3, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(ATP ~ gamma3) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

############ Experiment 4 -- ATP x gamma3 influence
### manuscript: Fig 3B

df.4 = read.csv("grn-sim-4.csv")
g.4 = ggplot(df.4, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(ATP ~ cd2) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

png("expt-4.png", width=800*sf, height=400*sf, res=72*sf)
print( ggarrange(g.3, g.4, labels=c("A", "B")) )
dev.off()

############ Experiment 5 -- different initial conditions
### manuscript: Supp Fig XXXX

df.5 = read.csv("grn-sim-5.csv")
g.5 = ggplot(df.5, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(ATP ~ IC.labels[ICs+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

png("expt-5.png", width=600*sf, height=300*sf, res=72*sf)
print( g.5 )
dev.off()

############ Experiment 6 -- different modes of ATP influence
### manuscript: Supp Fig XXXX

df.6a = read.csv("grn-sim-6a.csv")
df.6b = read.csv("grn-sim-6b.csv")
df.6c = read.csv("grn-sim-6c.csv")

g.6a = ggplot(df.6a, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(ATP ~ gamma3) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))
g.6b = ggplot(df.6b, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(ATP ~ gamma3) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))
g.6c = ggplot(df.6c, aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  theme(legend.position="none") + facet_grid(ATP ~ gamma3) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

png("expt-6.png", width=800*sf, height=300*sf, res=72*sf)
print( ggarrange(g.6a, g.6b, g.6c, labels=c("A", "B", "C")) )
dev.off()

############ Experiment 7 -- bifurcation plots
### manuscript: XXX

df.7 = read.csv("grn-sim-7.csv")

sub.7.1 = df.7[df.7$gamma3 == max(df.7$gamma3),]
sub.7.2 = df.7[df.7$gamma3 != max(df.7$gamma3),]
g.7.1 = ggplot(sub.7.1, aes(x=ATP, y=p1)) + geom_point()
g.7.2 = ggplot(sub.7.2, aes(x=ATP, y=p1)) + geom_point()

png("expt-7.png", width=600*sf, height=300*sf, res=72*sf)
print( ggarrange(g.7.1, g.7.2, labels=c("A", "B")) )
dev.off()

############ Experiment 8 -- matrices of pairwise changes
### manuscript: Supp Fig XXXX

df.8 = read.csv("grn-sim-8.csv", header=FALSE)
colnames(df.8) <- colnames(df.7)

g.8.1.1 = ggplot(df.8[df.8$scale==0,], aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  facet_grid(param.labels[param+1] ~ param.labels[param2+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))
g.8.1.2 = ggplot(df.8[df.8$scale==1,], aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  facet_grid(param.labels[param+1] ~ param.labels[param2+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))
g.8.1.3 = ggplot(df.8[df.8$scale==2,], aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  facet_grid(param.labels[param+1] ~ param.labels[param2+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))
g.8.1.4 = ggplot(df.8[df.8$scale==3,], aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  facet_grid(param.labels[param+1] ~ param.labels[param2+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

g.8.2.1 = ggplot(df.8[df.8$scale==4,], aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  facet_grid(param.labels[param+1] ~ param.labels[param2+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))
g.8.2.2 = ggplot(df.8[df.8$scale==5,], aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  facet_grid(param.labels[param+1] ~ param.labels[param2+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))
g.8.2.3 = ggplot(df.8[df.8$scale==6,], aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  facet_grid(param.labels[param+1] ~ param.labels[param2+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))
g.8.2.4 = ggplot(df.8[df.8$scale==7,], aes(x=ip1,y=ip2,xend=p1,yend=p2,color=p1)) + geom_segment() +
  facet_grid(param.labels[param+1] ~ param.labels[param2+1]) + coord_cartesian(xlim = c(0,30), ylim = c(0,30))

png("expt-8.1.png", width=1200*sf, height=900*sf, res=72*sf)
print( ggarrange(g.8.1.1, g.8.1.2, g.8.1.3, g.8.1.4, labels=c("A", "B", "C", "D")) )
dev.off()

png("expt-8.2.png", width=1200*sf, height=900*sf, res=72*sf)
print( ggarrange(g.8.2.1, g.8.2.2, g.8.2.3, g.8.2.4, labels=c("A", "B", "C", "D")) )
dev.off()


