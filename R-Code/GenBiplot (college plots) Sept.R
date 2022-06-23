library(ggplot2)
library(ggrepel)
library(lhs)
library(directlabels)

my.distance = "manhattan"
my.power = 1

raw.data = read.csv("~/Dropbox/Leman Research/Bi-plots/College_Data.csv")
data = raw.data[,colnames(raw.data)%in%c("Fac_to_Student","Tot_Enrollment","Percent_Grad_Student",
                                         "ACT_75","Percent_Admit",
                                         "Graduation_Rate","Percent_Male","Avg_Net_Cost")]
observations = raw.data[,1]
attributes.all = c("Stud/Fac","Enroll","GradStud",
                   "ACT","Admit","GradRate","Male",
                   "AvgCost")
X = scale(data)
n = nrow(X)
p = ncol(X)

Z.mds = MDS(X,HD.distance=my.distance,m.power=my.power)
my.plot.seq = seq(-5,5,by=.1)
fit = MDS.G.Biplot(X=X,Z=Z.mds,HD.distance=my.distance,m.power=my.power,plot.seq=my.plot.seq)
LD.axis.array = fit$LD.axis.array

my.colors = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')
names(my.colors) = attributes.all
my.colors2 = c(my.colors,"white")
names(my.colors2) = c(attributes.all,"")

K = 8 # How many axes to show?
ind = which(rank(rowMeans(fit$stress.mat))<(K+1))

# Compile LD axis array into a single matrix
attributes = attributes.all[ind]
LD.axis = cbind(c(apply(LD.axis.array[,,ind],3,function(x)x[,1])),c(apply(LD.axis.array[,,ind],3,function(x)x[,2])))
attribute.vec = rep(attributes,length(my.plot.seq))[order(rep(1:K,length(my.plot.seq)))]

# Define data frame with LD projection (Z.mds.df)
Z.mds.df = as.data.frame(Z.mds); colnames(Z.mds.df) = c("Z1","Z2")

# Define data frame with LD axes
LD.axis.df = as.data.frame(LD.axis)
LD.axis.df$attribute = attribute.vec
LD.axis.df$seq = rep(my.plot.seq,K); colnames(LD.axis.df) = c("b1","b2","attribute","seq")

# Define data frame with coordinates where I want tick marks
# tick.seq = seq(min(my.plot.seq),max(my.plot.seq),by=2); tick.seq = tick.seq[tick.seq!=0]
tick.seq = seq(min(my.plot.seq),max(my.plot.seq),by=2); tick.seq = tick.seq[tick.seq!=0]
LD.ticks = cbind(c(apply(LD.axis.array[which(my.plot.seq%in%tick.seq),,ind],3,function(x)x[,1])),
                 c(apply(LD.axis.array[which(my.plot.seq%in%tick.seq),,ind],3,function(x)x[,2])))
LD.ticks.df = as.data.frame(LD.ticks)
LD.ticks.df$seq = rep(tick.seq,K)
LD.ticks.df$attribute = rep(attributes,length(tick.seq))[order(rep(1:K,length(tick.seq)))]
colnames(LD.ticks.df) = c("b1","b2","seq","attribute")

# Define data frame to tell me where to put axis labels
endpoint.ind = seq(length(my.plot.seq),K*length(my.plot.seq),by=length(my.plot.seq))
Z.copy = Z.mds.df
Z.copy$label = ""
LD.axis.copy = LD.axis.df[,1:2]
LD.axis.copy$label = LD.axis.df$attribute
LD.axis.copy$label[-endpoint.ind] = ""
LD.axis.copy$label[endpoint.ind] = attributes
colnames(LD.axis.copy) = c("Z1","Z2","label")
label.df = as.data.frame(rbind(Z.copy,LD.axis.copy))

#Z.ps = 5
#Z.ts = 8
#b.ps = 8
#b.ts = 12

Z.ps = 2
Z.ts = 3
b.ps = 3
b.ts = 5

mdsplot <- ggplot() + labs(title=NULL,x="",y="") +
  geom_path(data=LD.axis.df,aes(x=b1,y=b2,color=factor(attribute)),linetype="dashed") +
  geom_text(data=LD.ticks.df,aes(x=b1,y=b2,label=seq,color=factor(attribute)),size=b.ps) +
  geom_label_repel(data=label.df,aes(x=Z1,y=Z2,
                   label=label,group=label,color=label),size=b.ts,fontface="bold") +
  scale_fill_manual(breaks=factor(attributes),values=my.colors2) +
  scale_color_manual(breaks=factor(attributes),values=my.colors2) +
  geom_point(data=Z.mds.df,aes(x=Z1,y=Z2),color="black",size=Z.ps) +
  geom_text_repel(data=Z.mds.df,aes(x=Z1,y=Z2,label=observations),color="black",size=Z.ts) + 
  scale_x_continuous(expand = c(.2,.2)) + 
   theme(axis.title.x=element_blank(),axis.text.x=element_blank(), 
        axis.ticks=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(),legend.position="none",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1))

# path = "~/Dropbox/Leman Research/General Dynamics/Bi-Plots/Manuscript/"
# if(my.distance=="sqrt.manhattan"){name="sqrtmanhattan"}else{name=my.distance}
# file.name = paste("GenBiplot_",name,".png",sep="")
# 
# png(file = paste(path,file.name,sep=""), width = 1024, height = 768, units = "px")
# print(mdsplot)
# dev.off()
print(mdsplot)
sum((as.matrix(dist(X,"manhattan"))-as.matrix(dist(Z.mds)))^2)

# head(LD.axis.df)
# 
par(mfrow=c(3,3))
for(i in 1:p)
{
  plot(LD.axis.array[,,i],main=attributes.all[i])
}

stress.mj = rep(NA,p)
for(i in 1:p)
{
 stress.mj[i] = sum((dist(X[,-i],"manhattan")-dist(Z.mds))^2) 
}
stress.all = sum((dist(X,"manhattan")-dist(Z.mds))^2)
