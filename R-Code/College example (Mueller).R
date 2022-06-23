library(ggplot2)
library(ggrepel)

my.distance = "Cosine"
my.plot.seq = seq(-4,4,by=1)

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

mueller_results <- mueller_fnc(filepath.and.name = file.path, raw.var.names = raw.col.names,plot.var.names = attributes.all,
                                high.D.obs.dist.metric = my.distance, #other option is Manhattan
                                mds.method = c("Euclidean"))
n = nrow(data)
p = ncol(data)

Z.mds = mueller_results$LD.pts
LD.axis.array = mueller_results$var.pts

attributes = attributes.all
LD.axis = cbind(c(apply(LD.axis.array,3,function(x)x[,1])),c(apply(LD.axis.array,3,function(x)x[,2])))
attribute.vec = attributes

Z.mds.df = as.data.frame(Z.mds); colnames(Z.mds.df) = c("Z1","Z2")
LD.axis.df = as.data.frame(LD.axis)
LD.axis.df$attribute = attribute.vec; colnames(LD.axis.df) = c("b1","b2","attribute")

mdsplot <- ggplot() + labs(title=NULL,x="",y="") +
  geom_point(data=Z.mds.df,aes(x=Z1,y=Z2),color="black") +
  geom_point(data=LD.axis.df,aes(x=b1,y=b2,color=factor(attributes))) +
  geom_text_repel(data=Z.mds.df,aes(x=Z1,y=Z2,label=observations),color="black",cex=2.5) + 
  geom_text_repel(data=LD.axis.df,aes(x=b1,y=b2,color=factor(attribute),label=attribute),cex=3) +
  scale_x_continuous(expand = c(.2, .2)) +
  # theme(axis.title.x=element_blank(),axis.text.x=element_blank(), 
  #       axis.ticks=element_blank()) +
  # theme(axis.title.y=element_blank(),axis.text.y=element_blank(), 
  #       axis.ticks=element_blank()) +
  theme(legend.position="none") 

path = "~/Dropbox/Leman Research/General Dynamics/Bi-Plots/Manuscript/"
file.name = "Mueller_Cosine.png"

png(file = paste(path,file.name,sep=""), width = 400, height = 400, units = "px")
print(mdsplot)
dev.off()
