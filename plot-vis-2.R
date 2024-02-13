library(reshape)
library(ggplot2)
library(ggpubr)

param.names = c("gamma2", "gamma3", "gamma4", "cd1", "cd2", "ca1", "ca2", "lambda2", "lambda3")
scale.names = c("x0.01", "x0.1", "x1", "x10", "x100")

# read in time series from gillespie simulation set
df = data.frame()
for(param in 0:8) {
  for(scale in 0:4) {
    tdf = read.csv(paste0("gillespie-scan-series-", 10*param+scale, ".csv")); tdf$scale = scale.names[scale+1]; tdf$param = param.names[param+1]
    df = rbind(df, tdf)
  }
}

# reshape data to be labelled by changed parameter and scale of change
df.long = melt(df, id.vars=c("t", "scale", "param"))

# plot of traces as time series
gill.traces = ggplot(df.long[df.long$variable %in% c("p_1", "p_2"),], 
                     aes(x=t, y=value, color=factor(variable))) + geom_line() + labs(color="Protein") +
  ylab("Count") +
  facet_grid(param~scale, scales="free")

# trellis plot of behaviours in p1, p2 space
gill.space = ggplot(df, aes(x=log10(p_1+1), y=log10(p_2+1), color=param)) + geom_path(alpha=0.2) +
  facet_grid(param~scale) + theme_classic()  +
  theme(legend.position = "none")

# read in data on switching time from gillespie simulations
df = data.frame()
for(param in 0:8) {
  for(scale in 0:4) {
    tdf = read.csv(paste0("gillespie-scan-switches-", 10*param+scale, ".csv")); 
    if(nrow(tdf) > 0) {
      tdf$scale = scale.names[scale+1]; tdf$param = param.names[param+1]
    df = rbind(df, tdf)
    }
  }
}

# plot histograms of switching times in same format as previous plot
gill.switches = ggplot(df, aes(x=log10(dt+1), fill=param)) + 
  geom_histogram(aes(y = ..density..) , alpha=0.4,position="identity") + 
  theme_light() + facet_grid(param~scale) + xlab("log10(switch time)") +
  theme(legend.position = "none")

# output both to a file
sf = 2
png("gillespie-scan-out.png", width=900*sf, height=900*sf, res=72*sf)
print(ggarrange(gill.space, gill.switches, labels=c("A", "B")))
dev.off()

# read in data from pairwise parameter changes and set header line and param names
df = read.csv("param-matrix-scan.csv", header=FALSE)
colnames(df) = c("n", "param", "param2", "s1", "s2", "P.gamma2", "P.gamma3", "P.gamma4", "P.cd1", "P.cd2", "P.ca1", "P.ca2", "P.lambda2", "P.lambda3", "ip1", "ip2", "p1", "p2")
df$param.name = factor(param.names[df$param+1], levels=param.names)
df$param2.name = factor(param.names[df$param2+1], levels=param.names)

# produce plots for different scales of change
gn.1.1 = ggplot(df[df$n==2 & df$s1==0.5 & df$s2 == 0.5,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~param2.name) 
gn.1.2 = ggplot(df[df$n==2 & df$s1==0.5 & df$s2 == 2,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~param2.name) 
gn.2.1 = ggplot(df[df$n==2 & df$s1==2 & df$s2 == 0.5,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~param2.name) 
gn.2.2 = ggplot(df[df$n==2 & df$s1==2 & df$s2 == 2,], aes(x=ip1,y=ip2,xend=p1,yend=p2, color=factor(round(p2, digits=1)))) + 
  geom_segment(size=0.5) + 
  theme(legend.position = "none") + 
  facet_grid(param.name~param2.name)

# put together in file output
sf = 2
png("param-matrix-scan.png", width=1000*sf, height=1000*sf, res=72*sf)
ggarrange(gn.1.1, gn.1.2, gn.2.1, gn.2.2, labels = c("A. x0.5,x0.5", "B. x0.5,x2", "C. x2,x0.5", "D. x2,x2"))
dev.off()
