library(ggplot2)
library(ggpubr)

# read output from simulation code
df = read.csv("hill.csv")

# get concentration points at which promoter is 50% active in each case
scale.1 = df[which(df[,2]>0.5)[1],1]
scale.2 = df[which(df[,3]>0.5)[1],1]
scale.4 = df[which(df[,4]>0.5)[1],1]

# unscale profile plot
g.0 = ggplot(df) + geom_line(aes(x=x0,y=p11), color="red") + 
  geom_line(aes(x=x0,y=p12), color="green") + 
  geom_line(aes(x=x0,y=p14), color="blue")  +
  xlab("x0") + ylab("p1") + theme_light()

# scaled profile plot
g.1 = ggplot(df) + geom_line(aes(x=x0/scale.1,y=p11), color="red") + 
  geom_line(aes(x=x0/scale.2,y=p12), color="green") + 
  geom_line(aes(x=x0/scale.4,y=p14), color="blue")  + xlim(0,2.5) +
  xlab("Scaled x0") + ylab("p1") + theme_light()
  
png("plot-hill.png", width=600*sf, height=200*sf, res=72*sf)
ggarrange(g.0, g.1, labels = c("A", "B"), nrow=1, ncol=2)
dev.off()
