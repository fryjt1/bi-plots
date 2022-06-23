library(ggplot2)
library(ggrepel)
library(lhs)

raw.data = read.csv("~/Dropbox/Leman Research/Bi-plots/College_Data.csv")
data = raw.data[,colnames(raw.data)%in%c("Fac_to_Student","Tot_Enrollment","Percent_Grad_Student",
                                         "ACT_75","Percent_Admit",
                                         "Graduation_Rate","Percent_Male","Avg_Net_Cost")]
observations = raw.data[,1]
attributes.all = c("Faculty-to-Student","Enrollment","Percent Graduate Students",
                   "ACT Score","Percent Admitted","Graduation Rate","Percent Male",
                   "Average Cost")
X = scale(data)
n = nrow(X)
p = ncol(X)
my.power = 2
Z.mds = Minkowski.MDS.SA(X,my.power=my.power)
Z.mds = scale(Rotation(Z.mds),scale=FALSE)
my.plot.seq = seq(-4,4,by=1); my.plot.seq = my.plot.seq[my.plot.seq!=0]
LD.axis.array = MDS.G.Biplot(X=X,Z=Z.mds,HD.distance=c("euclidean"),m.power=2,plot.seq=my.plot.seq)

LD.axis.array.all = array(NA,c(length(plot.seq),2,p))
stress.mat = matrix(NA,nrow=p,ncol=length(plot.seq))

for(i in 1:p)
{
  startx = c(0,0)
  for(j in 1:length(plot.seq))
  {
    starts = rbind(12*randomLHS(n=10,k=2)-6,startx)
    sols = matrix(NA,nrow=nrow(starts),ncol=2)
    mins = rep(NA,nrow(starts))
    
    for(k in 1:nrow(starts))
    {
      a = rep(0,p)
      a[i] = plot.seq[j]
      opt = optim(par=starts[k,],fn=Minkowski.Stress,method="L-BFGS-B",
                  lower=rep(-10,2),upper=rep(10,2),X=X,Z=Z.mds,a=a,my.power=my.power)
      sols[k,] = opt$par
      mins[k] = opt$value
    }
    LD.axis.array.all[j,,i] = sols[which.min(mins),]
    startx = sols[which.min(mins),]
    stress.mat[i,j] = mins[which.min(mins)]
  }
  print(i)
}

attributes = attributes.all
LD.axis = cbind(c(apply(LD.axis.array,3,function(x)x[,1])),c(apply(LD.axis.array,3,function(x)x[,2])))
attribute.vec = rep(attributes,dim(LD.axis.array)[1])[order(rep(1:p,dim(LD.axis.array)[1]))]

Z.mds.df = as.data.frame(Z.mds); colnames(Z.mds.df) = c("Z1","Z2")
LD.axis.df = as.data.frame(LD.axis)
LD.axis.df$attribute = attribute.vec; colnames(LD.axis.df) = c("b1","b2","attribute")

endpoint.ind = seq(length(my.plot.seq),p*length(my.plot.seq),by=length(my.plot.seq))
endpoint.df = LD.axis[endpoint.ind,]

mdsplot <- ggplot() + labs(title=paste("MDS (Minkowski, p =",my.power,") Projection"),x="",y="") +
  geom_point(data=Z.mds.df,aes(x=Z1,y=Z2),color="black") +
  geom_text_repel(data=Z.mds.df,aes(x=Z1,y=Z2,label=observations),color="black",cex=2.5) + 
  geom_path(data=LD.axis,aes(x=b1,y=b2,color=factor(attribute)),linetype="dotted") +
  geom_text_repel(data=endpoint.df,aes(x=b1,y=b2,color=factor(attribute),label=attribute),cex=3) +
  geom_text(data=LD.ticks,aes(x=V1,y=V2,label=mark,color=factor(attribute)),cex=2.5) +
  scale_x_continuous(expand = c(.2, .2)) +
  theme(legend.position="none")


print(mdsplot)
