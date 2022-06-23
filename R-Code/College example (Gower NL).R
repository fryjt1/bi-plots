library(ggplot2)
library(ggrepel)

my.distance = "Euclidean"
my.plot.seq = seq(-5,5,by=1)

file.path = "~/Dropbox/Leman Research/Bi-plots/College_Data.csv"
raw.data = read.csv(file.path)
raw.col.names = c("Fac_to_Student","Tot_Enrollment","Percent_Grad_Student",
                  "ACT_75","Percent_Admit",
                  "Graduation_Rate","Percent_Male","Avg_Net_Cost")
n = nrow(data)
p = ncol(data)
data = raw.data[,colnames(raw.data)%in%raw.col.names]
observations = raw.data[,1]
attributes.all = c("Student-to-Faculty","Enrollment","Percent Graduate Students",
                   "ACT Score","Percent Admitted","Graduation Rate","Percent Male",
                   "Average Cost")

gower_nl_results <- gower_nl_fnc(filepath.and.name = file.path,raw.var.names = raw.col.names, plot.var.names = attributes.all,
                                 high.D.dist.metric = my.distance, plot.seq = my.plot.seq)

n = nrow(data)
p = ncol(data)

Z.mds = gower_nl_results$LD.pts
LD.axis.array = gower_nl_results$axis.pts

attributes = attributes.all
LD.axis = cbind(c(apply(LD.axis.array,3,function(x)x[,1])),c(apply(LD.axis.array,3,function(x)x[,2])))
attribute.vec = rep(attributes,dim(LD.axis.array)[1])[order(rep(1:p,dim(LD.axis.array)[1]))]

Z.mds.df = as.data.frame(Z.mds); colnames(Z.mds.df) = c("Z1","Z2")
LD.axis.df = as.data.frame(LD.axis)
LD.axis.df$attribute = attribute.vec
LD.axis.df$seq = rep(my.plot.seq,p); colnames(LD.axis.df) = c("b1","b2","attribute","seq")

tick.seq = seq(min(my.plot.seq),max(my.plot.seq),by=1); tick.seq = tick.seq[tick.seq!=0]
LD.ticks = cbind(c(apply(LD.axis.array[which(my.plot.seq%in%tick.seq),,],3,function(x)x[,1])),
                 c(apply(LD.axis.array[which(my.plot.seq%in%tick.seq),,],3,function(x)x[,2])))
LD.ticks.df = as.data.frame(LD.ticks)
LD.ticks.df$seq = rep(tick.seq,p)
LD.ticks.df$attribute = rep(attributes,length(tick.seq))[order(rep(1:p,length(tick.seq)))]
colnames(LD.ticks.df) = c("b1","b2","seq","attribute")

endpoint.ind = seq(length(my.plot.seq),p*length(my.plot.seq),by=length(my.plot.seq))
endpoint.df = as.data.frame(LD.axis[endpoint.ind,])
endpoint.df$attribute = attributes; colnames(endpoint.df) = c("b1","b2","attribute")

mdsplot <- ggplot() + labs(title=NULL,x="",y="") +
  geom_point(data=Z.mds.df,aes(x=Z1,y=Z2),color="black") +
  geom_text_repel(data=Z.mds.df,aes(x=Z1,y=Z2,label=observations),color="black",cex=2.5) + 
  geom_path(data=LD.axis.df,aes(x=b1,y=b2,color=factor(attribute)),linetype="dotted") +
  geom_text_repel(data=endpoint.df,aes(x=b1,y=b2,color=factor(attribute),label=attribute),cex=3) +
  geom_text(data=LD.ticks.df,aes(x=b1,y=b2,label=seq,color=factor(attribute)),cex=2.5) +
  scale_x_continuous(expand = c(.2, .2)) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), 
        axis.ticks=element_blank()) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(), 
        axis.ticks=element_blank()) +
  theme(legend.position="none") 

path = "~/Dropbox/Leman Research/General Dynamics/Bi-Plots/Manuscript/"
file.name = "NL_Euclidean.png"

png(file = paste(path,file.name,sep=""), width = 400, height = 400, units = "px")
print(mdsplot)
dev.off()
