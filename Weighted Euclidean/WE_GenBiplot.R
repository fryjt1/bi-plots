library(ggplot2)
library(ggrepel)
library(lhs)
library(directlabels)
library(smacof)

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
HD.dist = as.matrix(dist(X,my.distance))

out = optim(par=rep(.5,p),fn=W.Euclidean.stress,X.mat=X,HD.dist=HD.dist,method="L-BFGS-B",
            lower=rep(0,p),upper=rep(10,p))
w.opt = out$par
D.curl = W.Euclidean(w.opt,X)

Z.mds = Rotation(General.MDS.SA(D=D.curl,bounds=c(-9,9)))
my.plot.seq = seq(-3,3,by=.5)
# fit = MDS.G.Biplot(X=X,Z=Z.mds,HD.distance=my.distance,m.power=my.power,plot.seq=my.plot.seq)
fit = MDS.G.Biplot.WE(X=X,Z=Z.mds,w=w.opt,plot.seq=my.plot.seq)
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

Z.ps = 2
Z.ts = 2
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
  scale_x_continuous(expand = c(.2,.2))  +
theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
     axis.ticks=element_blank(),axis.title.y=element_blank(),
     axis.text.y=element_blank(),legend.position="none",
     panel.background = element_blank(),
     panel.border = element_rect(colour = "black",fill=NA, size=1))

print(mdsplot)

# 
# path = "~/Dropbox/Leman Research/General Dynamics/Bi-Plots/Manuscript/"
# if(my.distance=="sqrt.manhattan"){name="sqrtmanhattan"}else{name=my.distance}
# file.name = paste("GenBiplot_",name,".png",sep="")
# 
# png(file = paste(path,file.name,sep=""), width = 1024, height = 768, units = "px")
# print(mdsplot)
# dev.off()
# print(mdsplot)
# 
# V = eigen(t(X)%*%X)$vectors[,1:2]
# lens = apply(V,1,function(c) sqrt(sum(c^2)))
# maxstress = apply(fit$stress.mat,1,max)
# 
# cbind(attributes,apply(LD.axis.array,3,function(c) sum(diff(c)^2)))
# sum(diff(LD.axis.array[,,1])^2)
# plot(lens,apply(LD.axis.array,3,function(c) sum(diff(c)^2)))
# 
# plot(NA,NA,xlim=range(LD.axis.array),ylim=range(LD.axis.array))
# cols = rainbow(p)
# for(i in 1:p)
# {
#   text(LD.axis.array[,1,i],LD.axis.array[,2,i],labels=round(fit$stress.mat[i,],0),cex=.75,col=cols[i])
# }
