library(ggplot2)
library(ggrepel)
library(lhs)
library(directlabels)

my.distance = "cosine"
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

attributes = attributes.all
LD.axis = cbind(c(apply(LD.axis.array,3,function(x)x[,1])),c(apply(LD.axis.array,3,function(x)x[,2])))
attribute.vec = attributes

Z.mds.df = as.data.frame(Z.mds); colnames(Z.mds.df) = c("Z1","Z2")
LD.axis.df = as.data.frame(LD.axis)
LD.axis.df$attribute = attribute.vec; colnames(LD.axis.df) = c("b1","b2","attribute")

# Define data frame to tell me where to put axis labels
Z.copy = Z.mds.df
Z.copy$label = ""
LD.axis.copy = LD.axis.df[,1:2]
LD.axis.copy$label = LD.axis.df$attribute
colnames(LD.axis.copy) = c("Z1","Z2","label")
label.df = as.data.frame(rbind(Z.copy,LD.axis.copy))

my.colors = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')
names(my.colors) = attributes.all
my.colors2 = c(my.colors,"white")
names(my.colors2) = c(attributes.all,"")

Z.ps = 5
Z.ts = 8
b.ps = 8
b.ts = 12

mdsplot <- ggplot() + labs(title=NULL,x="",y="") +
  geom_point(data=Z.mds.df,aes(x=Z1,y=Z2),color="black",size=Z.ps) +
  geom_point(data=LD.axis.df,aes(x=b1,y=b2,color=factor(attributes)),size=b.ps) +
  geom_text_repel(data=Z.mds.df,aes(x=Z1,y=Z2,label=observations),color="black",size=Z.ts) + 
  geom_label_repel(data=LD.axis.df,aes(x=b1,y=b2,color=factor(attribute),label=attribute),size=b.ts) +
  scale_fill_manual(breaks=factor(attributes),values=my.colors2) +
  scale_color_manual(breaks=factor(attributes),values=my.colors2) +
  scale_x_continuous(expand = c(.2, .2)) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), 
        axis.ticks=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(),legend.position="none",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1))

# path = "~/Dropbox/Leman Research/General Dynamics/Bi-Plots/Manuscript/"
# if(my.distance=="sqrt.manhattan"){name="SqrtManhattan"}else{name=my.distance}
# file.name = paste("GenBiplot_",name,".png",sep="")
# 
# png(file = paste(path,file.name,sep=""), width = 1024, height = 768, units = "px")
# print(mdsplot)
# dev.off()
print(mdsplot)
