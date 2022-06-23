library(ggplot2)
library(ggrepel)
library(lhs)

I = 10
my.distance = "manhattan"
my.power = 1
n = 25
p = 3
n.bad = 1

storage = matrix(NA,I,p)
sd.store = matrix(NA,I,p)

system.time({
for(i in 1:I)
{
  sd.vec = c(runif(p-n.bad,.5,1),runif(n.bad,0,.5))
  sd.store[i,] = sd.vec
  X = matrix(rnorm(n*p),n,p)
  X = scale(X)
  X = scale(X,scale=1/sd.vec)
  
  Z.mds = MDS(X,scale=FALSE,HD.distance=my.distance,m.power=my.power)
  my.plot.seq = seq(-5,5,by=.5)
  fit = MDS.G.Biplot(X=X,Z=Z.mds,HD.distance=my.distance,m.power=my.power,plot.seq=my.plot.seq)
  LD.axis.array = fit$LD.axis.array
  stress = fit$stress
  storage[i,] = apply(LD.axis.array,3,function(c) sqrt(sum((c[nrow(c),]-c[1,])^2)))
  print(i)
}
})

storage.df = as.data.frame(storage)
storage.df = stack(storage.df)
colnames(storage.df) = c("Stress","Attribute")
p1 = ggplot(aes(y = Stress, x=factor(Attribute)), data = storage.df) + geom_boxplot() +
  labs(x="Attribute",y="Average Stress") +
  scale_x_discrete(labels=c("1","2","3")) +
  theme(axis.text=element_text(size=27),
        axis.title=element_text(size=27),panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1)) 
print(p1)
path = "~/Dropbox/Leman Research/General Dynamics/Bi-Plots/Manuscript/"
file.name = paste("bad_axis_sim",".png",sep="")
png(file = paste(path,file.name,sep=""), width = 1024, height = 768, units = "px")
print(p1)
dev.off()


observations = 1:n
attributes.all = 1:p
my.colors = c("red","blue","black")
names(my.colors) = attributes.all
my.colors2 = c(my.colors,"white")
names(my.colors2) = c(attributes.all,"")

# Compile LD axis array into a single matrix
attributes = attributes.all
LD.axis = cbind(c(apply(LD.axis.array,3,function(x)x[,1])),c(apply(LD.axis.array,3,function(x)x[,2])))
attribute.vec = rep(attributes,dim(LD.axis.array)[1])[order(rep(1:p,dim(LD.axis.array)[1]))]

# Define data frame with LD projection (Z.mds.df)
Z.mds.df = as.data.frame(Z.mds); colnames(Z.mds.df) = c("Z1","Z2")

# Define data frame with LD axes
LD.axis.df = as.data.frame(LD.axis)
LD.axis.df$attribute = attribute.vec
LD.axis.df$seq = rep(my.plot.seq,p); colnames(LD.axis.df) = c("b1","b2","attribute","seq")

# Define data frame with coordinates where I want tick marks
tick.seq = seq(min(my.plot.seq),max(my.plot.seq),by=2); tick.seq = tick.seq[tick.seq!=0]
LD.ticks = cbind(c(apply(LD.axis.array[which(my.plot.seq%in%tick.seq),,],3,function(x)x[,1])),
                 c(apply(LD.axis.array[which(my.plot.seq%in%tick.seq),,],3,function(x)x[,2])))
LD.ticks.df = as.data.frame(LD.ticks)
LD.ticks.df$seq = rep(tick.seq,p)
LD.ticks.df$attribute = rep(attributes,length(tick.seq))[order(rep(1:p,length(tick.seq)))]
colnames(LD.ticks.df) = c("b1","b2","seq","attribute")

# Define data frame to tell me where to put axis labels
endpoint.ind = seq(length(my.plot.seq),p*length(my.plot.seq),by=length(my.plot.seq))
Z.copy = Z.mds.df
Z.copy$label = ""
LD.axis.copy = LD.axis.df[,1:2]
LD.axis.copy$label = LD.axis.df$attribute
LD.axis.copy$label[-endpoint.ind] = ""
LD.axis.copy$label[endpoint.ind] = attributes
colnames(LD.axis.copy) = c("Z1","Z2","label")
label.df = as.data.frame(rbind(Z.copy,LD.axis.copy))

mdsplot <- ggplot() + labs(title=NULL,x="",y="") +
  geom_path(data=LD.axis.df,aes(x=b1,y=b2,color=factor(attribute)),linetype="dashed") +
  geom_text(data=LD.ticks.df,aes(x=b1,y=b2,label=seq,color=factor(attribute)),cex=5) +
  geom_label_repel(data=label.df,aes(x=Z1,y=Z2,
                                     label=label,group=label,color=label),size=8,fontface="bold") +
  scale_fill_manual(breaks=factor(attributes),values=my.colors2) +
  scale_color_manual(breaks=factor(attributes),values=my.colors2) +
  geom_point(data=Z.mds.df,aes(x=Z1,y=Z2),color="black",cex=3) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(),legend.position="none",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1))
print(mdsplot)

# path = "~/Dropbox/Leman Research/General Dynamics/Bi-Plots/Manuscript/"
# file.name = paste("bad_axis_sim_example",".png",sep="")
# png(file = paste(path,file.name,sep=""), width = 1024, height = 768, units = "px")
# print(mdsplot)
# dev.off()